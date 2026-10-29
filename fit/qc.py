"""Data checks before fitting. Lists suspicious points; drops nothing.

  1. jumps in the slab's spin moment along R within a height (the cell total minus
     the two H; a free H shows 0.381 muB inside its sphere for 1 muB): leave-one-out
     spline error > MAG_JUMP, so smooth trends are not flagged
  2. outliers in the surface effect E_pbe(h, R) - E_pbe(8, R): leave-one-out cubic
     spline error > LOO_DEV; subtracting the h = 8 row removes the steep wall
  3. baseline assumptions: V_H(z >= 6) ~ 0, gas-phase H2 De / r_e from the h = 8 row

Run after collect_data.py:  python3 qc.py
"""
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

MAG_JUMP = 0.10     # muB, leave-one-out deviation
H_SPHERE = 0.381    # muB inside the H sphere for one free-H spin
LOO_DEV = 0.010     # eV

df = pd.read_csv("pes_data.csv")
nflag = 0

print("== 1. slab-moment jumps (leave-one-out error > %.2f muB)" % MAG_JUMP)
df["m_slab"] = df.mag - (df.mH1 + df.mH2) / H_SPHERE
for h, g in df.groupby("h"):
    g = g.sort_values("R")
    x, ms = g.R.values, g.m_slab.values
    for i in range(1, len(x) - 1):
        dev = ms[i] - CubicSpline(np.delete(x, i), np.delete(ms, i))(x[i])
        if abs(dev) > MAG_JUMP:
            print(f"  h {h:.2f}: R {x[i]:.2f}  m_slab {ms[i]:+.3f}  ({dev:+.3f} off its neighbours)")
            nflag += 1

print("== 2. outliers in E_pbe(h, R) - E_pbe(8, R): leave-one-out spline error > "
      f"{1000 * LOO_DEV:.0f} meV")
gas = df[np.isclose(df.h, 8)].set_index("R").E_pbe
for h, g in df.groupby("h"):
    if np.isclose(h, 8):
        continue
    g = g.sort_values("R")
    x, y = g.R.values, g.E_pbe.values - gas.loc[g.R].values
    loo = np.array([y[i] - CubicSpline(np.delete(x, i), np.delete(y, i))(x[i])
                    for i in range(1, len(x) - 1)])
    bad = np.abs(loo) > LOO_DEV
    print(f"  h {h:.2f}: rms {1000 * np.sqrt(np.mean(loo ** 2)):.1f} meV, {bad.sum()} flagged"
          + "".join(f"\n      R {r:.2f}: {1000 * v:+.1f} meV"
                    for r, v in zip(x[1:-1][bad], loo[bad])))
    nflag += bad.sum()

print("== 3. baseline cuts")
ref = df[np.isclose(df.h, 8) & np.isclose(df.R, 8)].E_pbe.iloc[0]
col = df[np.isclose(df.R, 8)].sort_values("h")
for h, e in zip(col.h, col.E_pbe - ref):
    tag = "  <- should be ~0" if h >= 6 and abs(e) > 0.01 else ""
    print(f"  V_H(z = {h:.2f}) = {e:+.4f} eV{tag}")
row = df[np.isclose(df.h, 8)].sort_values("R")
m = (row.R >= 0.65) & (row.R <= 0.85)
c = np.polyfit(row.R[m], row.E_pbe[m], 4)
rr = np.linspace(0.65, 0.85, 20001)
ee = np.polyval(c, rr)
print(f"  gas-phase H2 (h = 8, PBE): r_e = {rr[np.argmin(ee)]:.3f} A, "
      f"De = {ref - ee.min():.3f} eV")

print(f"\n{nflag} flag(s). Inspect flagged points before trusting the fit there.")
