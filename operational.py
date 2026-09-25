import sys, numpy as np
sys.path.insert(0, '..')
import copula_check as cc
from scipy.optimize import brentq
def series(path):
    rows, T, D = cc.load_series(path); return rows, np.asarray(T, float), np.asarray(D, float)
def gumbel_cml(T, D):
    return cc.fit_copula("Gumbel", cc.pseudo_ranks(T), cc.pseudo_ranks(D))[0][0]
out = {}
for name, f in [("Sutcliffe", "sutcliffe.csv"), ("Chikatilo", "chikatilo.csv"), ("Ridgway", "ridgway_disposal.csv")]:
    rows, T, D = series(f); n = len(T)
    th = max(1.0, gumbel_cml(T, D))
    q = lambda x, p: float(np.quantile(x, p))            # empirical quantiles (linear interpolation)
    t50, t85, t90, r80 = q(T, .5), q(T, .85), q(T, .9), q(D, .8)
    C = lambda u, v: float(cc.cdf_gumbel(np.array([u]), np.array([v]), th)[0])
    def rcond(u):   # distance containing the next displacement with prob. 0.8 given U > u
        g = lambda v: (v - C(u, v)) / (1 - u) - 0.8
        v = brentq(g, 0.5, 1 - 1e-9); return q(D, v), v
    r85, v85 = rcond(.85); r90, v90 = rcond(.90)
    pexc = 1 - (0.8 - C(.9, .8)) / 0.1                      # P(D > r80 | T > t90) under the fitted copula
    pp = lambda u: (1 - 2*u + C(u, u)) / (1 - u)             # P(V>u | U>u)
    sel = T > t85; emp = (np.sum(D[sel] > r80), np.sum(sel))
    print(f"{name:10s} n={n} theta_CML={th:.3f} t50={t50:.0f} t85={t85:.0f} ({t85/t50:.1f}x) t90={t90:.0f} r80={r80:.1f} "
          f"r80|t85={r85:.1f} r80|t90={r90:.1f} P(D>r80|T>t90)={pexc:.2f} P(V>.9|U>.9)={pp(.9):.2f} "
          f"empirical: {emp[0]} of {emp[1]} pauses > t85 followed by D > r80")
    # longest intervals and their displacements
    idx = np.argsort(-T)[:4]
    print("   longest intervals:", "; ".join(f"{rows[i][3]}->{rows[i+1][3]} T={T[i]:.0f} D={D[i]:.1f}" for i in idx))
    idx = np.argsort(-D)[:4]
    print("   largest displacements:", "; ".join(f"{rows[i][3]}->{rows[i+1][3]} T={T[i]:.0f} D={D[i]:.1f}" for i in idx))
