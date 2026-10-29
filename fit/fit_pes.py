"""Fit the 2D PES (pes_model.PES), cross-validate, check gradients, save the model.

Cross-validation (errors on PBE+D3, binned by energy above gas-phase H2):
  heights : leave whole heights out (CV_HEIGHTS present in the data; if none are,
            any interior height whose removal leaves a gap <= 1 A)
  R       : leave every other R point out at all heights
The GP input transform ('none' or 'exp') with the lower CV error is used for the
final model.

Run after collect_data.py:  python3 fit_pes.py   ->  pes_model.joblib, fit_report.txt
"""
import time
import numpy as np
import pandas as pd
import joblib
from pes_model import PES, H_INF, R_INF

CV_HEIGHTS = [1.35, 1.45, 1.55, 1.65, 2.10, 3.50]
WINDOWS = [1.0, 3.0, 5.0, np.inf]       # eV above gas-phase H2
TRANSFORMS = [("none", 1.0), ("exp", 1.0)]

df = pd.read_csv("pes_data.csv")
log = []


def out(s=""):
    print(s, flush=True)
    log.append(s)


def errors(model, test):
    e_dft = model.energy(test.R.values, test.h.values, d3=True)
    e_ref = (test.E_pbe - model.e_ref + test.Edisp - model.d_ref).values
    err = e_dft - e_ref
    above = e_ref - model.gas_ref(True)
    res = []
    for w in WINDOWS:
        m = above < w
        res.append((w, m.sum(), np.sqrt(np.mean(err[m] ** 2)) if m.any() else np.nan,
                    np.max(np.abs(err[m])) if m.any() else np.nan))
    return res


def report(name, res):
    out(f"  {name}")
    for w, n, rmse, mx in res:
        lab = f"E < {w:g} eV" if np.isfinite(w) else "all"
        out(f"    {lab:11s} n = {n:4d}   RMSE {1000 * rmse:7.1f} meV   max {1000 * mx:7.1f} meV")


def keep_baseline(test_mask):
    """Never hold out the h = 8 row or the R = 8 column: they define the baseline."""
    return test_mask & ~np.isclose(df.h, H_INF) & ~np.isclose(df.R, R_INF)


hs = np.array(sorted(df.h.unique()))
cv_h = [h for h in CV_HEIGHTS if np.any(np.isclose(hs, h))]
if not cv_h:
    cv_h = [hs[k] for k in range(1, len(hs) - 1) if hs[k + 1] - hs[k - 1] <= 1.0 + 1e-9]
out(f"data: {len(df)} points, heights " + " ".join(f"{h:.2f}" for h in hs))
out(f"CV heights left out: " + (" ".join(f"{h:.2f}" for h in cv_h) or "none available"))

Rs = np.array(sorted(df.R.unique()))
odd_R = Rs[1:-1:2]
scores = {}
for tr, lam in TRANSFORMS:
    out(f"\n== transform '{tr}'")
    tot = []
    for h in cv_h:
        test = keep_baseline(np.isclose(df.h, h))
        t0 = time.time()
        m = PES(tr, lam, restarts=3).fit(df[~test])
        r = errors(m, df[test])
        report(f"leave out h = {h:.2f}  ({time.time() - t0:.0f} s)", r)
        tot.append(r[1][2])
    test = keep_baseline(np.isin(np.round(df.R, 2), np.round(odd_R, 2)))
    m = PES(tr, lam, restarts=3).fit(df[~test])
    r = errors(m, df[test])
    report("leave out every other R", r)
    tot.append(r[1][2])
    scores[tr] = np.nanmean(tot)
    out(f"  mean RMSE (E < 3 eV): {1000 * scores[tr]:.1f} meV")

best = min(scores, key=scores.get)
lam = dict(TRANSFORMS)[best]
out(f"\n== final model: transform '{best}' (lam = {lam}), all data")
pes = PES(best, lam, restarts=10).fit(df)
out(f"  PBE 3-body GP : {pes.gp.kernel_}")
out(f"  D3 GP         : {pes.gp_d3.kernel_}")
report("training error", errors(pes, df))

# gradients: analytic vs central finite differences
rng = np.random.default_rng(1)
R = rng.uniform(0.3, 8.0, 1000)
h = rng.uniform(hs[0], hs[-1], 1000)
g = pes.gradient(R, h)
eps = 1e-4
fd = np.column_stack([(pes.energy(R + eps, h) - pes.energy(R - eps, h)) / (2 * eps),
                      (pes.energy(R, h + eps) - pes.energy(R, h - eps)) / (2 * eps)])
out(f"\n  gradient check (1000 random points): max |analytic - FD| = "
    f"{1000 * np.max(np.abs(g - fd)):.3f} meV/A")

# gas-phase limit from the model
mu = 0.5 * 1.00794 * 1.66054e-27
for lab, d3 in (("PBE+D3", True), ("PBE   ", False)):
    rr = np.linspace(0.65, 0.85, 2001)
    c = np.polyfit(rr, pes.energy(rr, H_INF, d3=d3), 4)
    r2 = np.linspace(0.65, 0.85, 200001)
    e2 = np.polyval(c, r2)
    k = np.argmin(e2)
    w = np.sqrt(np.polyval(np.polyder(c, 2), r2[k]) * 1.602177e-19 / 1e-20 / mu) / (2 * np.pi * 2.99792458e10)
    out(f"  gas-phase H2 from the model ({lab}): r_e = {r2[k]:.4f} A, De = {-e2[k]:.3f} eV, "
        f"omega_e = {w:.0f} cm-1")

# behaviour off the grid
out("  off-grid cuts (eV, PBE+D3, relative to the atoms; sigma = 3-body GP, 0 where tapered off):")
for lab, RR, HH in [("R = 0.15, h = 1.0", 0.15, 1.0), ("R = 1.0,  h = 0.8", 1.0, 0.8),
                    ("R = 0.15, h = 8  ", 0.15, 8.0), ("R = 0.75, h = 12 ", 0.75, 12.0),
                    ("R = 12,  h = 1.5 ", 12.0, 1.5), ("R = 12,  h = 12  ", 12.0, 12.0)]:
    out(f"    {lab} E = {float(pes.energy(RR, HH)):9.3f}   sigma = {1000 * float(pes.sigma(RR, HH)):6.1f} meV")

joblib.dump(pes, "pes_model.joblib")
open("fit_report.txt", "w").write("\n".join(log) + "\n")
out("\nwrote pes_model.joblib, fit_report.txt")
