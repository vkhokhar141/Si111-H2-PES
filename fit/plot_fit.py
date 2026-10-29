"""Plot the fitted PES with the DFT points.

  fit_pes_2d.png : contour map of the fit E(R, h) (PBE+D3, relative to gas-phase H2)
                   with the data points; the fit is not drawn across height gaps
                   > MAXGAP (no data there, the fit is a guess)
  fit_cuts.png   : one panel per computed height, fit (line) vs data (dots)
  fit_pes_3d.png : 3D surface of the fit with the data points, same blank gaps

Run after fit_pes.py:  python3 plot_fit.py
"""
import numpy as np
import pandas as pd
import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

ECAP = 8.0
MAXGAP = 1.01      # A; larger height gaps are left blank on the map
TEXT, TEXT2, GRID, NODATA = "#0b0b0b", "#52514e", "#e4e3df", "#f1f0ec"
BLUE, DOT = "#2a78d6", "#0b0b0b"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("blue_seq", BLUES)

pes = joblib.load("pes_model.joblib")
df = pd.read_csv("pes_data.csv")
df["E"] = df.E_pbe - pes.e_ref + df.Edisp - pes.d_ref - pes.gas_ref(True)
df["fit"] = pes.energy(df.R.values, df.h.values, ref="gas")
H = np.array(sorted(df.h.unique()))
gaps = [(H[k - 1], H[k]) for k in range(1, len(H)) if H[k] - H[k - 1] > MAXGAP]

plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.linewidth": 0.8})

# ---------------- 2D map
RMAX = 6.0
Rg = np.linspace(0.25, RMAX, 300)
hg = np.linspace(H[0], H[-1], 400)
RR, HH = np.meshgrid(Rg, hg)
EE = pes.energy(RR, HH, ref="gas")
levels = np.arange(0, ECAP + 1e-9, 0.5)
norm = BoundaryNorm(levels, CMAP.N, extend="both")

for lo, hi in gaps:
    EE[(HH > lo + 1e-6) & (HH < hi - 1e-6)] = np.nan
fig, ax = plt.subplots(figsize=(7.6, 6.6))
ax.set_facecolor(NODATA)
cf = ax.contourf(RR, HH, np.clip(EE, -0.5, None), levels=levels, cmap=CMAP, norm=norm, extend="both")
ax.contour(RR, HH, EE, levels=np.arange(1, ECAP + 1e-9, 1), colors="white", linewidths=0.8)
for lo, hi in gaps:
    ax.text(RMAX / 2, (lo + hi) / 2, f"no data between {lo:.2f} and {hi:.2f} Å\n(fit not drawn)",
            color=TEXT2, fontsize=9, ha="center", va="center", zorder=3)
d = df[df.R <= RMAX]
ax.scatter(d.R, d.h, s=9, c="white", edgecolors=DOT, linewidths=0.6, zorder=4, label="DFT points")
ax.text(0.78, 7.55, "H$_2$(g)", color=TEXT, fontsize=9, zorder=5)
ax.text(4.75, 1.72, "H + H–Si  ↓", color=TEXT, fontsize=9, zorder=5, ha="center",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
ax.set_xlim(0.25, RMAX)
ax.set_ylim(H[0], H[-1])
ax.set_xlabel("H–H distance R (Å)")
ax.set_ylabel("height of lower H above Si (Å)")
ax.set_title("Fitted 2D PES (GPR, PBE+D3) with DFT points", color=TEXT, fontsize=11, loc="left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=True, facecolor="white", edgecolor="none", loc="upper right", fontsize=8.5,
          labelcolor=TEXT).set_zorder(6)
cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=np.arange(0, ECAP + 1e-9, 1))
cb.set_label("E − E(H$_2$ gas)  (eV)", color=TEXT)
cb.outline.set_visible(False)
fig.savefig("fit_pes_2d.png", dpi=200, bbox_inches="tight", facecolor="white")

# ---------------- cuts per height
n = len(H)
ncol = 4
nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(12, 2.7 * nrow), sharex=True, sharey=True)
Rf = np.linspace(0.25, 8.0, 800)
YLO, YHI = -0.6, 8.0
for ax, h in zip(axs.flat, H):
    g = df[np.isclose(df.h, h)].sort_values("R")
    ax.plot(Rf, pes.energy(Rf, h, ref="gas"), color=BLUE, lw=1.8, label="fit")
    on = g.E <= YHI
    ax.scatter(g.R[on], g.E[on], s=14, c="white", edgecolors=DOT, linewidths=0.7, zorder=3, label="DFT")
    m = g.E < 3.0 + g.E.min()
    rms = 1000 * np.sqrt(np.mean((g.fit - g.E)[m] ** 2))
    ax.text(0.97, 0.95, f"h = {h:.2f} Å\nfit error {rms:.1f} meV", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color=TEXT)
    ax.axhline(0, color=TEXT2, lw=0.6)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
for ax in axs.flat[n:]:
    ax.set_visible(False)
axs.flat[0].set_ylim(YLO, YHI)
axs.flat[0].legend(frameon=False, loc="center right", fontsize=8.5, labelcolor=TEXT)
for ax in axs[-1]:
    ax.set_xlabel("R (Å)")
for ax in axs[:, 0]:
    ax.set_ylabel("E − E(H$_2$ gas) (eV)")
fig.suptitle("Fit vs DFT at each computed height  (fit error = RMS over points within 3 eV of that "
             "height's minimum; points above 8 eV off scale)", x=0.01, ha="left", fontsize=10, color=TEXT)
fig.tight_layout(rect=(0, 0, 1, 1 - 0.35 / fig.get_figheight()))
fig.savefig("fit_cuts.png", dpi=170, bbox_inches="tight", facecolor="white")
# ---------------- 3D surface
R3MIN, Z3MAX = 0.55, 9.0      # every data point in this R range is below Z3MAX (max 8.2 eV)
blocks, cur = [], [H[0]]
for k in range(1, len(H)):
    if H[k] - H[k - 1] > MAXGAP:
        blocks.append(cur)
        cur = []
    cur.append(H[k])
blocks.append(cur)
fig = plt.figure(figsize=(9, 7.2))
ax = fig.add_subplot(projection="3d", computed_zorder=False)
R3 = np.linspace(R3MIN, RMAX, 160)
for b in blocks:
    if len(b) < 2:
        continue
    hb = np.linspace(b[0], b[-1], max(20, int((b[-1] - b[0]) / 0.02)))
    RR3, HH3 = np.meshgrid(R3, hb)
    EE3 = np.clip(pes.energy(RR3, HH3, ref="gas"), -0.5, Z3MAX)
    ax.plot_surface(RR3, HH3, EE3, cmap=CMAP, norm=norm, rstride=2, cstride=2, linewidth=0,
                    antialiased=True, alpha=0.85, zorder=1)
d3 = df[(df.R >= R3MIN - 1e-9) & (df.R <= RMAX + 1e-9)]
ax.scatter(d3.R, d3.h, d3.E, s=7, c=DOT, depthshade=False, zorder=3, label="DFT points")
ax.set_xlim(R3MIN, RMAX)
ax.set_ylim(H[0], H[-1])
ax.set_zlim(-0.5, Z3MAX)
ax.set_xlabel("R (Å)")
ax.set_ylabel("height h (Å)")
ax.set_zlabel("E − E(H$_2$ gas) (eV)")
ax.view_init(elev=24, azim=-128)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor("white")
    a.pane.set_edgecolor(GRID)
    a._axinfo["grid"]["color"] = GRID
ax.set_title("Fitted 3D PES (GPR, PBE+D3) with DFT points", color=TEXT, fontsize=11, loc="left")
ax.legend(frameon=False, loc="upper right", fontsize=8.5, labelcolor=TEXT)
gtxt = "; ".join(f"{lo:.2f}–{hi:.2f} Å" for lo, hi in gaps)
fig.text(0.02, 0.03, f"surface drawn only where heights were computed (no data {gtxt}); "
         f"R ≥ {R3MIN} Å shown, all points in this range are below {Z3MAX:g} eV",
         color=TEXT2, fontsize=7.5)
fig.savefig("fit_pes_3d.png", dpi=200, bbox_inches="tight", facecolor="white")
print("wrote fit_pes_2d.png, fit_cuts.png, fit_pes_3d.png")
