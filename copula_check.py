#!/usr/bin/env python3
"""
copula_check.py -- reproducible (re)analysis for
"Copula-based analysis of space-time dependence in serial homicide".

INPUT  : a CSV file with one row per LINKED offence and columns
             label,date,lat,lon
         (date as YYYY-MM-DD; decimal degrees). Row order is irrelevant.
USAGE  : python3 copula_check.py series.csv [--B 1000] [--seed 1] [--zero_km 0.1]
         python3 copula_check.py pairs.csv --pairs      (CSV with columns T,D)
         python3 copula_check.py --params mu_T sigma_T xi_T mu_D sigma_D theta [--name X]
         (the second form only prints the operational table for given fitted parameters)

OUTPUT : T_i (days) and D_i (km, haversine, R = 6371 km); descriptive statistics;
         Kendall/Spearman tests; GPD & Weibull fits for T, exponential & log-normal for D
         (AIC and BIC with the SAME parameter counts, so BIC - AIC = p(ln n - 2));
         IFM and CML fits of the Gaussian, Student-t, Clayton and Gumbel copulas;
         rank-based parametric-bootstrap Cramer-von Mises test (Genest & Remillard, 2008);
         Gumbel diagnostics; empirical joint-exceedance counts; operational table.
"""
import argparse, csv, sys
from datetime import date
import numpy as np
from scipy import stats, optimize, special

R_EARTH = 6371.0
EPS = 1e-10

# ----------------------------------------------------------------------------- data
def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    h = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(h))

def load_series(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            y, m, d = map(int, r["date"].split("-"))
            rows.append((date(y, m, d), float(r["lat"]), float(r["lon"]), r.get("label", "")))
    rows.sort(key=lambda x: x[0])
    T = np.array([(rows[i + 1][0] - rows[i][0]).days for i in range(len(rows) - 1)], float)
    D = np.array([haversine(rows[i][1], rows[i][2], rows[i + 1][1], rows[i + 1][2])
                  for i in range(len(rows) - 1)])
    return rows, T, D

def load_pairs(path):
    """Alternative input: CSV with columns T,D (already computed intervals and displacements)."""
    T, D = [], []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            T.append(float(r["T"])); D.append(float(r["D"]))
    n = len(T)
    rows = [(None, None, None, f"pair{i+1}") for i in range(n + 1)]
    return rows, np.array(T), np.array(D)

# ----------------------------------------------------------------------------- copulas
def _clip(u):
    return np.clip(u, EPS, 1 - EPS)

def logc_gauss(u, v, rho):
    x, y = stats.norm.ppf(_clip(u)), stats.norm.ppf(_clip(v))
    r2 = 1 - rho ** 2
    return -0.5 * np.log(r2) - (rho ** 2 * (x ** 2 + y ** 2) - 2 * rho * x * y) / (2 * r2)

def logc_t(u, v, rho, nu):
    x, y = stats.t.ppf(_clip(u), nu), stats.t.ppf(_clip(v), nu)
    r2 = 1 - rho ** 2
    lj = (special.gammaln((nu + 2) / 2) - special.gammaln(nu / 2) - np.log(nu * np.pi)
          - 0.5 * np.log(r2) - (nu + 2) / 2 * np.log1p((x ** 2 - 2 * rho * x * y + y ** 2) / (nu * r2)))
    return lj - stats.t.logpdf(x, nu) - stats.t.logpdf(y, nu)

def logc_clayton(u, v, th):
    u, v = _clip(u), _clip(v)
    return (np.log1p(th) - (th + 1) * (np.log(u) + np.log(v))
            - (2 + 1 / th) * np.log(u ** -th + v ** -th - 1))

def logc_gumbel(u, v, th):
    u, v = _clip(u), _clip(v)
    x, y = -np.log(u), -np.log(v)
    s = x ** th + y ** th
    A = s ** (1 / th)
    return (-A - np.log(u) - np.log(v) + (th - 1) * (np.log(x) + np.log(y))
            + (1 / th - 2) * np.log(s) + np.log(A + th - 1))

# Gauss-Legendre nodes for C(u,v) = int_0^u h(v|s) ds (Gaussian and t copulas)
_GL_X, _GL_W = np.polynomial.legendre.leggauss(96)

def _cdf_by_h(u, v, hfun):
    u, v = np.atleast_1d(u), np.atleast_1d(v)
    s = 0.5 * (_GL_X[None, :] + 1) * u[:, None]
    w = 0.5 * _GL_W[None, :] * u[:, None]
    return np.sum(w * hfun(np.broadcast_to(v[:, None], s.shape), s), axis=1)

def cdf_gauss(u, v, rho):
    def h(vv, ss):
        return stats.norm.cdf((stats.norm.ppf(_clip(vv)) - rho * stats.norm.ppf(_clip(ss)))
                              / np.sqrt(1 - rho ** 2))
    return _cdf_by_h(u, v, h)

def cdf_t(u, v, rho, nu):
    def h(vv, ss):
        x1, x2 = stats.t.ppf(_clip(ss), nu), stats.t.ppf(_clip(vv), nu)
        return stats.t.cdf((x2 - rho * x1) / np.sqrt((nu + x1 ** 2) * (1 - rho ** 2) / (nu + 1)), nu + 1)
    return _cdf_by_h(u, v, h)

def cdf_clayton(u, v, th):
    u, v = _clip(u), _clip(v)
    return np.maximum(u ** -th + v ** -th - 1, 0) ** (-1 / th)

def cdf_gumbel(u, v, th):
    u, v = _clip(u), _clip(v)
    return np.exp(-(((-np.log(u)) ** th + (-np.log(v)) ** th) ** (1 / th)))

def sample(fam, par, n, rng):
    if fam == "Gaussian":
        z = rng.multivariate_normal([0, 0], [[1, par[0]], [par[0], 1]], n)
        return stats.norm.cdf(z[:, 0]), stats.norm.cdf(z[:, 1])
    if fam == "Student-t":
        rho, nu = par
        z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], n)
        w = np.sqrt(rng.chisquare(nu, n) / nu)
        return stats.t.cdf(z[:, 0] / w, nu), stats.t.cdf(z[:, 1] / w, nu)
    if fam == "Clayton":
        th = par[0]
        u, w = rng.random(n), rng.random(n)
        v = (u ** -th * (w ** (-th / (1 + th)) - 1) + 1) ** (-1 / th)
        return u, v
    if fam == "Gumbel":
        th = par[0]
        if th <= 1 + 1e-8:
            return rng.random(n), rng.random(n)
        a = 1 / th                                  # positive stable (Chambers-Mallows-Stuck)
        U0, E = rng.uniform(0, np.pi, n), rng.exponential(1, n)
        S = (np.sin(a * U0) / np.sin(U0) ** (1 / a)) * (np.sin((1 - a) * U0) / E) ** ((1 - a) / a)
        E1, E2 = rng.exponential(1, n), rng.exponential(1, n)
        return np.exp(-(E1 / S) ** a), np.exp(-(E2 / S) ** a)
    raise ValueError(fam)

FAMILIES = {
    "Gaussian":  dict(p=1, x0=[0.0], bounds=[(-0.98, 0.98)],
                      logc=lambda u, v, q: logc_gauss(u, v, q[0]), cdf=lambda u, v, q: cdf_gauss(u, v, q[0])),
    "Student-t": dict(p=2, x0=[0.0, 8.0], bounds=[(-0.98, 0.98), (1.5, 60.0)],
                      logc=lambda u, v, q: logc_t(u, v, q[0], q[1]), cdf=lambda u, v, q: cdf_t(u, v, q[0], q[1])),
    "Clayton":   dict(p=1, x0=[0.5], bounds=[(1e-4, 30.0)],
                      logc=lambda u, v, q: logc_clayton(u, v, q[0]), cdf=lambda u, v, q: cdf_clayton(u, v, q[0])),
    "Gumbel":    dict(p=1, x0=[1.3], bounds=[(1.0 + 1e-6, 30.0)],
                      logc=lambda u, v, q: logc_gumbel(u, v, q[0]), cdf=lambda u, v, q: cdf_gumbel(u, v, q[0])),
}

def fit_copula(fam, u, v):
    F = FAMILIES[fam]
    nll = lambda q: -np.sum(F["logc"](u, v, q))
    best = None
    starts = [F["x0"]] + ([[0.5, 5.0], [-0.3, 10.0]] if fam == "Student-t" else
                          [[0.5]] if fam == "Gaussian" else [[3.0]] if fam != "Clayton" else [[2.0]])
    for x0 in starts:
        r = optimize.minimize(nll, x0, bounds=F["bounds"], method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun

def empirical_copula_at_points(u, v):
    return np.array([np.mean((u <= ui) & (v <= vi)) for ui, vi in zip(u, v)])

def cvm_stat(fam, par, u, v):
    return np.sum((empirical_copula_at_points(u, v) - FAMILIES[fam]["cdf"](u, v, par)) ** 2)

def pseudo_ranks(x):
    return stats.rankdata(x) / (len(x) + 1)

def gof_bootstrap(fam, u, v, B, rng):
    par, _ = fit_copula(fam, u, v)
    s_obs = cvm_stat(fam, par, u, v)
    n, count = len(u), 0
    for _ in range(B):
        us, vs = sample(fam, par, n, rng)
        us, vs = pseudo_ranks(us), pseudo_ranks(vs)
        pb, _ = fit_copula(fam, us, vs)
        count += cvm_stat(fam, pb, us, vs) >= s_obs
    return par, s_obs, (count + 1) / (B + 1)

def tau_of(fam, par):
    if fam in ("Gaussian", "Student-t"):
        return 2 / np.pi * np.arcsin(par[0])
    if fam == "Clayton":
        return par[0] / (par[0] + 2)
    return 1 - 1 / par[0]

def lambda_u(fam, par):
    if fam == "Gumbel":
        return 2 - 2 ** (1 / par[0])
    if fam == "Student-t":
        rho, nu = par
        return 2 * stats.t.cdf(-np.sqrt((nu + 1) * (1 - rho) / (1 + rho)), nu + 1)
    return 0.0

# ----------------------------------------------------------------------------- operational table
def gumbel_cond_quantile(u_star, q, th):
    """v such that P(V <= v | U > u_star) = q under a Gumbel copula."""
    f = lambda v: (v - cdf_gumbel(np.array([u_star]), np.array([v]), th)[0]) / (1 - u_star) - q
    return optimize.brentq(f, 1e-12, 1 - 1e-12)

def operational_table(mu_T, sig_T, xi_T, mu_D, sig_D, th, name=""):
    gpd = stats.genpareto(xi_T, loc=mu_T, scale=sig_T)
    ln = stats.lognorm(sig_D, loc=0, scale=np.exp(mu_D))
    t50 = gpd.ppf(0.5)
    r80 = ln.ppf(0.8)
    out = dict(name=name, t50=t50, r50=ln.ppf(0.5), r80=r80)
    for p in (0.85, 0.90):
        v = gumbel_cond_quantile(p, 0.8, th)
        C = cdf_gumbel(np.array([p]), np.array([0.8]), th)[0]
        out[p] = dict(t=gpd.ppf(p), t_over_median=gpd.ppf(p) / t50, r80c=ln.ppf(v),
                      mult=ln.ppf(v) / r80, p_exceed=1 - (0.8 - C) / (1 - p))
    C9 = cdf_gumbel(np.array([0.9]), np.array([0.9]), th)[0]
    out["p_joint90"] = (1 - 2 * 0.9 + C9) / 0.1          # P(V > .9 | U > .9); 0.10 under independence
    return out

def print_operational(o):
    print(f"\nOPERATIONAL TABLE {o['name']}")
    print(f"  median interval {o['t50']:.0f} d | median displacement {o['r50']:.1f} km | habitual radius r80 {o['r80']:.1f} km")
    for p in (0.85, 0.90):
        d = o[p]
        print(f"  pause > P{int(p*100)} = {d['t']:.0f} d ({d['t_over_median']:.1f} x median): r80 -> {d['r80c']:.0f} km "
              f"(x{d['mult']:.1f}); P(next D > habitual r80) = {d['p_exceed']:.2f} (0.20 if independent)")
    print(f"  P(D above its P90 | T above its P90) = {o['p_joint90']:.2f} (0.10 if independent)")
    d85, d90 = o[0.85], o[0.90]
    print("  LaTeX row: {} & {:.0f} & {:.0f} & {:.0f} & {:.1f} & {:.0f} & {:.0f} & {:.2f} \\\\".format(
        o["name"], o["t50"], d85["t"], d90["t"], o["r80"], d85["r80c"], d90["r80c"], d90["p_exceed"]))

# ----------------------------------------------------------------------------- main analysis
def analyse(path, B, seed, zero_km, pairs=False):
    rng = np.random.default_rng(seed)
    rows, T, D = load_pairs(path) if pairs else load_series(path)
    n = len(T)
    print(f"File: {path}  offences: {len(rows)}  pairs n = {n}\n")
    print("  i  from -> to                    T_i(d)   D_i(km)")
    for i in range(n):
        print(f"{i+1:3d}  {rows[i][3][:12]:12s} -> {rows[i+1][3][:12]:12s} {T[i]:7.0f}  {D[i]:9.2f}")
    cv = lambda x: x.std(ddof=1) / x.mean()
    print(f"\nT: min {T.min():.0f}  mean {T.mean():.1f}  median {np.median(T):.1f}  max {T.max():.0f}  CV {cv(T):.2f}"
          f"  (sum {T.sum():.0f} d = span of series)")
    print(f"D: min {D.min():.2f}  mean {D.mean():.1f}  median {np.median(D):.1f}  max {D.max():.1f}  CV {cv(D):.2f}")
    if np.any(D <= 0):
        print(f"WARNING: {np.sum(D <= 0)} zero displacements (same coordinates); replaced by {zero_km} km "
              "for the log-normal fit. Report how zeros were handled.")
    Dp = np.where(D <= 0, zero_km, D)
    tau, p_tau = stats.kendalltau(T, D)
    rho, p_rho = stats.spearmanr(T, D)
    z = tau * np.sqrt(9 * n * (n - 1) / (2 * (2 * n + 5)))
    print(f"\nKendall tau = {tau:.3f}  z = {z:.2f}  p(asympt, 2-sided) = {2*stats.norm.sf(abs(z)):.4f}  p(scipy) = {p_tau:.4f}")
    print(f"Spearman rho = {rho:.3f}  p = {p_rho:.4f}")

    ln_n = np.log(n)
    print("\nMARGINALS (p = number of estimated parameters; BIC - AIC = p(ln n - 2))")
    xi, loc, sig = stats.genpareto.fit(T, floc=T.min())
    ll = np.sum(stats.genpareto.logpdf(T, xi, loc, sig)); p = 3
    print(f"  T  GPD      mu={loc:.1f} sigma={sig:.2f} xi={xi:.3f} (tail exponent 1+1/xi = {1+1/xi if xi>0 else np.inf:.2f})"
          f"  AIC {2*p-2*ll:.1f}  BIC {p*ln_n-2*ll:.1f}")
    Tw = np.where(T <= 0, 0.5, T)
    k, _, lam = stats.weibull_min.fit(Tw, floc=0)
    ll = np.sum(stats.weibull_min.logpdf(Tw, k, 0, lam)); p = 2
    print(f"  T  Weibull  k={k:.3f} eta={lam:.1f}  AIC {2*p-2*ll:.1f}  BIC {p*ln_n-2*ll:.1f}")
    _, sc = stats.expon.fit(Dp, floc=0)
    ll = np.sum(stats.expon.logpdf(Dp, 0, sc)); p = 1
    print(f"  D  Exponential lambda={1/sc:.4f}  AIC {2*p-2*ll:.1f}  BIC {p*ln_n-2*ll:.1f}")
    lmu, lsig = np.mean(np.log(Dp)), np.std(np.log(Dp))
    ll = np.sum(stats.lognorm.logpdf(Dp, lsig, 0, np.exp(lmu))); p = 2
    print(f"  D  Log-normal mu={lmu:.3f} sigma={lsig:.3f} (median {np.exp(lmu):.1f} km)  AIC {2*p-2*ll:.1f}  BIC {p*ln_n-2*ll:.1f}")

    U_ifm = _clip(stats.genpareto.cdf(T, xi, loc, sig))
    V_ifm = _clip(stats.lognorm.cdf(Dp, lsig, 0, np.exp(lmu)))
    U_r, V_r = pseudo_ranks(T), pseudo_ranks(D)

    print("\nCOPULAS  (IFM: parametric margins; CML: ranks; GoF: rank-based bootstrap, B = %d)" % B)
    print(f"  {'family':10s} {'IFM par':>16s} {'tau(par)':>8s} {'lamU':>6s} {'AIC':>7s} {'BIC':>7s} {'CML par':>16s} {'S_n':>7s} {'p':>6s}")
    for fam, F in FAMILIES.items():
        par, ll = fit_copula(fam, U_ifm, V_ifm)
        aic, bic = 2 * F["p"] - 2 * ll, F["p"] * ln_n - 2 * ll
        par_c, s_n, pval = gof_bootstrap(fam, U_r, V_r, B, rng)
        fmt = lambda q: ", ".join(f"{x:.3f}" for x in q)
        print(f"  {fam:10s} {fmt(par):>16s} {tau_of(fam, par):8.3f} {lambda_u(fam, par):6.3f} {aic:7.1f} {bic:7.1f}"
              f" {fmt(par_c):>16s} {s_n:7.4f} {pval:6.3f}")
    th_tau = 1 / (1 - tau) if tau > 0 else 1.0
    print(f"\n  Gumbel moment estimate theta_tau = 1/(1 - tau) = {th_tau:.3f}")
    for u0 in (0.8, 0.9):
        a = np.sum(U_r > u0); b = np.sum((U_r > u0) & (V_r > u0))
        print(f"  Empirical: of {a} intervals above their P{int(u0*100)}, {b} were followed by a displacement above its "
              f"P{int(u0*100)} (expected under independence: {a*(1-u0):.1f})")
    print("\n  Pseudo-observations for Fig. 1 (IFM):")
    print("  " + "".join(f"({a:.2f},{b:.2f})" for a, b in zip(U_ifm, V_ifm)))
    print("\n  Empirical copula on the diagonal (rank pseudo-observations):")
    for u0 in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95):
        print(f"   u = {u0:.2f}: C_n = {np.mean((U_r <= u0) & (V_r <= u0)):.3f}")
    return dict(T=T, D=D)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--zero_km", type=float, default=0.1)
    ap.add_argument("--params", type=float, nargs=6, metavar=("muT", "sigT", "xiT", "muD", "sigD", "theta"))
    ap.add_argument("--name", default="")
    ap.add_argument("--pairs", action="store_true", help="CSV has columns T,D instead of label,date,lat,lon")
    a = ap.parse_args()
    if a.params:
        print_operational(operational_table(*a.params, name=a.name))
    elif a.csv:
        analyse(a.csv, a.B, a.seed, a.zero_km, a.pairs)
    else:
        ap.print_help(); sys.exit(1)
