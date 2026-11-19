#!/usr/bin/env python3
"""
revision_checks.py -- analyses added in the second revision of
"Copula-based analysis of space-time dependence in serial homicide: evidence from three
high-volume, long-duration series" (requested by Reviewer 2).

  (1) Table S5: sensitivity of the parametric fits of the displacements and of the IFM copula
      estimates to the value assigned to zero displacements (offences recorded at the same
      locality or recovery cluster). Rank-based results do not use this value.
  (2) Table 7 (bootstrap 95% intervals) and Table S6: precision and stability of the
      operational parameters (median and 85th percentile of the intervals, habitual and
      conditional radii).

USAGE : python3 revision_checks.py      (run in the folder that contains the CSV files and
                                         copula_check.py)
OUTPUT: revision_checks.txt (log) and tables_S5_S6.tex (LaTeX tables used by make_supplement.py)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats, optimize
import copula_check as cc

B, SEED = 2000, 2026
ZERO_VALUES = (0.1, 0.5, 1.0, 2.0, 5.0)          # km; 1.0 km is the value used in the paper
CENSOR = (2.0, 5.0)                               # km; left-censoring thresholds (log-normal only)
SERIES = [("Sutcliffe", "sutcliffe.csv"), ("Chikatilo", "chikatilo.csv"), ("Ridgway", "ridgway_disposal.csv")]
log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
q = lambda x, p: float(np.quantile(x, p))         # empirical quantile (linear interpolation)
def ordinal(k):
    k = int(round(k))
    return f"{k}{'th' if 10 <= k % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(k % 10, 'th')}"

# ----------------------------------------------------------------------------- (1) zero displacements
def parametric_and_ifm(T, D, zero_km):
    Dp = np.where(D <= 0, zero_km, D)
    xi, loc, sig = stats.genpareto.fit(T, floc=T.min())
    _, sc = stats.expon.fit(Dp, floc=0)
    aic_exp = 2 - 2 * np.sum(stats.expon.logpdf(Dp, 0, sc))
    mu, s = np.mean(np.log(Dp)), np.std(np.log(Dp))
    aic_ln = 4 - 2 * np.sum(stats.lognorm.logpdf(Dp, s, 0, np.exp(mu)))
    U = cc._clip(stats.genpareto.cdf(T, xi, loc, sig))
    V = cc._clip(stats.lognorm.cdf(Dp, s, 0, np.exp(mu)))
    ifm = {fam: cc.fit_copula(fam, U, V)[0] for fam in cc.FAMILIES}
    return dict(mu=mu, sigma=s, aic_ln=aic_ln, aic_exp=aic_exp, ifm=ifm)

def censored_lognormal(D, c):
    """Log-normal MLE treating the zero displacements as left-censored below c km."""
    pos, k = D[D > 0], int(np.sum(D <= 0))
    def nll(par):
        mu, s = par[0], np.exp(par[1])
        return -(np.sum(stats.norm.logpdf(np.log(pos), mu, s) - np.log(pos))
                 + k * stats.norm.logcdf((np.log(c) - mu) / s))
    r = optimize.minimize(nll, [np.mean(np.log(pos)), np.log(np.std(np.log(pos)))], method="Nelder-Mead",
                          options=dict(xatol=1e-8, fatol=1e-10, maxiter=4000))
    return r.x[0], np.exp(r.x[1])

S5 = []
log("(1) ZERO DISPLACEMENTS")
for name, f in SERIES[1:]:
    rows, T, D = cc.load_series(f)
    zeros = [f"{rows[i][3]}-{rows[i+1][3]}" for i in range(len(D)) if D[i] <= 0]
    pos = np.sort(D[D > 0])
    tau, rho = stats.kendalltau(T, D)[0], stats.spearmanr(T, D)[0]
    cml = cc.fit_copula("Gumbel", cc.pseudo_ranks(T), cc.pseudo_ranks(D))[0][0]
    log(f"\n{name}: {len(zeros)} zero pairs ({'; '.join(zeros)}); smallest positive displacement {pos[0]:.2f} km")
    log(f"  rank-based (recorded zeros, tied): tau = {tau:.3f}, rho_S = {rho:.3f}, CML Gumbel = {cml:.3f}")
    rows_tex = []
    for z in ZERO_VALUES:
        r = parametric_and_ifm(T, D, z)
        g, t, cl, gu = r["ifm"]["Gaussian"], r["ifm"]["Student-t"], r["ifm"]["Clayton"], r["ifm"]["Gumbel"]
        log(f"  zero -> {z:3.1f} km: LN mu {r['mu']:.2f} sigma {r['sigma']:.2f} AIC {r['aic_ln']:.1f} | Exp AIC {r['aic_exp']:.1f}"
            f" | IFM: Ga {g[0]:.3f}, t {t[0]:.3f} ({t[1]:.1f}), Cl {cl[0]:.3f}, Gu {gu[0]:.3f}")
        bold = (lambda x: r"\textbf{" + x + "}") if z == 1.0 else (lambda x: x)
        rows_tex.append(" & ".join([bold(f"{z:g}~km"), bold(f"{r['mu']:.2f}"), bold(f"{r['sigma']:.2f}"),
                        bold(f"{r['aic_ln']:.1f}"), bold(f"{r['aic_exp']:.1f}"), bold(f"{g[0]:.3f}"),
                        bold(f"{t[0]:.3f} ({t[1]:.1f})"), bold(f"{cl[0]:.3f}"),
                        bold(f"{gu[0]:.3f}")]) + r"\\")
    for c in CENSOR:
        mu, s = censored_lognormal(D, c)
        log(f"  left-censored below {c:.0f} km: LN mu {mu:.2f} sigma {s:.2f}")
        rows_tex.append(f"censored $<{c:g}$~km & {mu:.2f} & {s:.2f} & --- & --- & --- & --- & --- & ---\\\\")
    S5.append((name, len(zeros), pos[0], tau, rho, cml, rows_tex))

# ----------------------------------------------------------------------------- (2) operational parameters
def conditional_radius(T, D, u=0.85, level=0.8):
    th = max(1.0, cc.fit_copula("Gumbel", cc.pseudo_ranks(T), cc.pseudo_ranks(D))[0][0])
    C = lambda a, b: float(cc.cdf_gumbel(np.array([a]), np.array([b]), th)[0])
    v = optimize.brentq(lambda b: (b - C(u, b)) / (1 - u) - level, 0.5, 1 - 1e-9)
    return q(D, v)

log(f"\n(2) OPERATIONAL PARAMETERS: point estimates, bootstrap 95% intervals (B = {B}, seed {SEED}) and stability")
rng = np.random.default_rng(SEED)
S6, T7 = [], []
for name, f in SERIES:
    rows, T, D = cc.load_series(f); n = len(T)
    est = dict(t50=q(T, .5), t85=q(T, .85), r80=q(D, .8), rc=conditional_radius(T, D))
    boot = {k: [] for k in est}
    for _ in range(B):
        i = rng.integers(0, n, n); Tb, Db = T[i], D[i]
        boot["t50"].append(q(Tb, .5)); boot["t85"].append(q(Tb, .85)); boot["r80"].append(q(Db, .8))
        boot["rc"].append(conditional_radius(Tb, Db))
    ci = {k: (np.quantile(v, .025), np.quantile(v, .975)) for k, v in boot.items()}
    log(f"\n{name} (n = {n}):")
    for k, lab in [("t50", "median interval t50 (d)"), ("t85", "85th percentile t85 (d)"),
                   ("r80", "habitual radius r80 (km)"), ("rc", "conditional radius r80(t85) (km)")]:
        log(f"  {lab:34s} {est[k]:8.1f}   95% bootstrap interval [{ci[k][0]:.1f}, {ci[k][1]:.1f}]")
    pct3 = np.mean(T <= 3 * est["t50"])
    log(f"  three times the median = {3*est['t50']:.0f} d, i.e. the {ordinal(100*pct3)} percentile of the observed intervals")
    seq = []
    for k in (15, 20, 25, 30, 40):
        if k - 1 < n:
            seq.append((k, q(T[:k-1], .85), 3 * q(T[:k-1], .5)))
    seq.append((n + 1, est["t85"], 3 * est["t50"]))
    log("  estimates available after k linked offences: " +
        "; ".join(f"k = {k}: t85 = {a:.0f}, 3 x median = {b:.0f}" for k, a, b in seq))
    S6.append((name, n, seq, pct3)); T7.append((name, est, ci))

# ----------------------------------------------------------------------------- LaTeX fragment
fmt = lambda x: f"{x:,.0f}".replace(",", "{,}") if x >= 100 else (f"{x:.0f}" if x >= 10 and float(x).is_integer() else f"{x:.1f}")
tex = [r"\begin{table}[htbp]", r"\centering\footnotesize\setlength{\tabcolsep}{4pt}",
       r"\caption{Sensitivity of the parametric fits of the displacements and of the IFM copula estimates to the value assigned to "
       r"zero displacements (consecutive offences recorded at the same locality or recovery cluster); the value used in the paper (1~km) is in bold. "
       r"`Censored' rows treat these displacements as left-censored below the stated value (log-normal fit only). "
       r"Rank-based results (Kendall's $\tau$, Spearman's $\rho_S$, CML estimates, goodness-of-fit tests and Table~7 of the paper) use the "
       r"recorded zeros as tied smallest values and do not depend on this choice.}",
       r"\begin{tabular}{lcccccccc}", r"\toprule",
       r" & \multicolumn{2}{c}{Log-normal} & \multicolumn{2}{c}{AIC} & \multicolumn{4}{c}{IFM copula estimates}\\",
       r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-9}",
       r"Zero displacement & $\hat\mu$ & $\hat\sigma$ & Log-normal & Exponential & Gaussian $\hat\rho$ & Student-$t$ $\hat\rho$ ($\hat\nu$) & Clayton $\hat\theta$ & Gumbel $\hat\theta$\\",
       r"\midrule"]
for k, (name, nz, p0, tau, rho, cml, rows_tex) in enumerate(S5):
    tex.append(r"\multicolumn{9}{l}{\textit{" + name + f"}} ({nz} zeros; smallest positive displacement {p0:.1f}~km; "
               f"$\\hat\\tau={tau:.3f}$, CML Gumbel $\\hat\\theta={cml:.3f}$)" + r"}\\")
    tex += rows_tex
    if k == 0:
        tex.append(r"\midrule")
tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
tex += [r"\begin{table}[htbp]", r"\centering\small",
        r"\caption{Stability of the trigger. Empirical 85th percentile of the intervals ($t_{85}$) and three times the median interval, "
        r"in days, as they would have been estimated after the first $k$ linked offences ($k-1$ intervals); the last column of each row gives the "
        r"values for the complete series and the percentile of the observed intervals to which three times the median corresponds.}",
        r"\begin{tabular}{llcccccc}", r"\toprule",
        r"Series & & $k=15$ & $k=20$ & $k=25$ & $k=30$ & $k=40$ & Complete series\\", r"\midrule"]
for name, n, seq, pct3 in S6:
    cells85, cells3 = [], []
    d = {k: (a, b) for k, a, b in seq}
    for k in (15, 20, 25, 30, 40):
        if k in d and k != n + 1:
            cells85.append(f"{d[k][0]:.0f}"); cells3.append(f"{d[k][1]:.0f}")
        else:
            cells85.append("---"); cells3.append("---")
    a, b = d[n + 1]
    tex.append(f"{name} ($n={n}$) & $t_{{85}}$ & " + " & ".join(cells85) + f" & {a:.0f}\\\\")
    tex.append(f" & $3\\times t_{{50}}$ & " + " & ".join(cells3) + f" & {b:.0f} ({ordinal(100*pct3)} percentile)\\\\")
tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
open("tables_S5_S6.tex", "w", encoding="utf-8").write("\n".join(tex))
log("\nTable 7 rows (estimate and bootstrap 95% interval):")
for name, est, ci in T7:
    log(f"  {name}: t50 {fmt(est['t50'])} ({fmt(ci['t50'][0])}--{fmt(ci['t50'][1])}); t85 {fmt(est['t85'])} ({fmt(ci['t85'][0])}--{fmt(ci['t85'][1])}); "
        f"r80 {est['r80']:.1f} ({ci['r80'][0]:.1f}--{ci['r80'][1]:.1f}); r80(t85) {est['rc']:.1f} ({ci['rc'][0]:.1f}--{ci['rc'][1]:.1f})")
open("revision_checks.txt", "w", encoding="utf-8").write("\n".join(log_lines) + "\n")
