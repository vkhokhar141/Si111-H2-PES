"""2D (contour) and 3D (surface) PES E0(R, h) from the data used for the fit.

Reads fit/pes_data.csv (scan values, replaced by the reruns where those were accepted,
excluded points removed; run fit/collect_data.py first). E0 = PBE+D3.
R = H-H distance, h = height of the lower H above Si atom 13.
Energies are relative to gas-phase H2: the R = 0.75 A point of the highest height.
Heights are interpolated only across gaps <= MAXGAP; larger gaps stay blank. A point
missing from a height (excluded) is filled for drawing only, by a cubic spline along R.

Run from the project directory:  python3 plot_pes2d3d.py   ->  pes_2d.png, pes_3d.png
"""
import os
import re
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

TOP = os.path.dirname(os.path.abspath(__file__))
MAXGAP = 0.51      # A; interpolate between heights at most this far apart
STRIP = 0.10       # A; half-width of the band drawn for an isolated height
ECAP = 8.0         # eV; colour scale / 3D axis top, higher values are capped
RMAX = 5.0         # A; the curves are flat (< 3 meV) beyond this
RMIN3D = 0.45      # A; 3D view starts here (R <= 0.35 is a > 9 eV wall)
ZMAX3D = 11.0      # eV; 3D axis top

TEXT, TEXT2, GRID, NODATA = "#0b0b0b", "#52514e", "#e4e3df", "#f1f0ec"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("blue_seq", BLUES)


def load_heights():
    """{h: array (R, E0)} on the common R grid; also the list of filled points."""
    df = pd.read_csv(os.path.join(TOP, "fit", "pes_data.csv"))
    grid = np.array(sorted(df.R.unique()))
    data, filled = {}, []
    for h, g in df.groupby("h"):
        g = g.sort_values("R")
        e = np.interp(grid, g.R, g.E0)
        miss = ~np.isin(np.round(grid, 3), np.round(g.R.values, 3))
        if miss.any():
            e[miss] = CubicSpline(g.R.values, g.E0.values)(grid[miss])
            filled += [(h, r) for r in grid[miss]]
        data[float(h)] = np.column_stack([grid, e])
    src = df.src.value_counts().to_dict()
    return data, filled, src


data, filled, src = load_heights()
H = np.array(sorted(data))
R = data[H[0]][:, 0]
E = np.array([data[h][:, 1] for h in H])                       # (nh, nR)
Eref = E[-1][np.argmin(abs(R - 0.75))]
E = E - Eref
print("heights:", " ".join(f"{h:.2f}" for h in H), f"| reference E0 = {Eref:.4f} eV")
src_txt = " + ".join(f"{n} {k}" for k, n in src.items())
fill_txt = ("; filled along R for drawing: " + ", ".join(f"{h:.2f}/R_{r:.2f}" for h, r in filled)
            if filled else "")

# fine height grid; NaN where no data
hf = np.round(np.arange(H[0] - STRIP, H[-1] + STRIP + 1e-9, 0.01), 3)
Ef = np.full((hf.size, R.size), np.nan)
for i, h in enumerate(hf):
    j = np.searchsorted(H, h)
    if 0 < j < H.size and H[j] - H[j - 1] <= MAXGAP:
        t = (h - H[j - 1]) / (H[j] - H[j - 1])
        Ef[i] = (1 - t) * E[j - 1] + t * E[j]
    else:
        k = np.argmin(abs(H - h))
        if abs(H[k] - h) <= STRIP + 1e-9:
            Ef[i] = E[k]

# contiguous height blocks (for the 3D surfaces)
blocks, cur = [], [0]
for k in range(1, H.size):
    if H[k] - H[k - 1] <= MAXGAP:
        cur.append(k)
    else:
        blocks.append(cur)
        cur = [k]
blocks.append(cur)

plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.linewidth": 0.8,
                     "hatch.color": "#d6d5cf", "hatch.linewidth": 0.6})

# ---------------- 2D contour map
levels = np.arange(0, ECAP + 1e-9, 0.5)
norm = BoundaryNorm(levels, CMAP.N, extend="max")
fig, ax = plt.subplots(figsize=(7.5, 6.2))
ax.patch.set_facecolor(NODATA)
ax.patch.set_hatch("///")
m = R <= RMAX + 1e-9
cf = ax.contourf(R[m], hf, np.clip(Ef[:, m], 0, None), levels=levels, cmap=CMAP,
                 norm=norm, extend="max")
ax.contour(R[m], hf, Ef[:, m], levels=np.arange(1, ECAP + 1e-9, 1), colors="white",
                linewidths=0.8, alpha=0.9)
for h in H:
    ax.plot([R[m][0], R[m][-1]], [h, h], color="white", lw=0.6, alpha=0.7)
ticks = [h for h in H if h > 2.21 or round(h * 100) % 20 == 0]   # every other below 2.2 A
ax.set_yticks(ticks)
ax.set_yticklabels([f"{h:.2f}" for h in ticks], fontsize=8)
ax.set_xlim(R[0], RMAX)
ax.set_ylim(hf[0], hf[-1])
ax.set_xlabel("H–H distance R (Å)")
ax.set_ylabel("height of lower H above Si (Å)")
ax.set_title("H$_2$ on Si(111): 2D PES  E(R, h) relative to gas-phase H$_2$",
             color=TEXT, fontsize=11, loc="left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
gaps = [(H[k - 1], H[k]) for k in range(1, H.size) if H[k] - H[k - 1] > MAXGAP]
for lo, hi in gaps:
    ax.text(RMAX / 2, (lo + hi) / 2, "no data", color=TEXT2, fontsize=9,
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.9))
cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=np.arange(0, ECAP + 1e-9, 1))
cb.set_label("E − E(H$_2$ gas)  (eV)", color=TEXT)
cb.outline.set_visible(False)
fig.text(0.01, -0.03, f"data used for the fit (fit/pes_data.csv, PBE+D3): {src_txt} points"
         f"{fill_txt}\nwhite lines: the {H.size} computed heights ({H[0]:.2f}-{H[-1]:.2f} Å)  ·  "
         f"interpolated between heights ≤ 0.5 Å apart  ·  colours capped at {ECAP:g} eV",
         color=TEXT2, fontsize=7.5)
fig.savefig(os.path.join(TOP, "pes_2d.png"), dpi=200, bbox_inches="tight", facecolor="white")

# ---------------- 3D surface
fig = plt.figure(figsize=(8.5, 6.8))
ax = fig.add_subplot(projection="3d", computed_zorder=False)
m = (R >= RMIN3D - 1e-9) & (R <= RMAX + 1e-9)
Rm = R[m]
for b in blocks:
    if len(b) < 2:
        continue
    hb = hf[(hf >= H[b[0]] - 1e-9) & (hf <= H[b[-1]] + 1e-9)]
    Eb = Ef[(hf >= H[b[0]] - 1e-9) & (hf <= H[b[-1]] + 1e-9)][:, m]
    RR, HH = np.meshgrid(Rm, hb)
    ax.plot_surface(RR, HH, np.clip(Eb, 0, ZMAX3D), cmap=CMAP, norm=norm, rstride=5, cstride=1,
                    linewidth=0, antialiased=True, alpha=0.95)
for k, h in enumerate(H):
    ax.plot(Rm, np.full(Rm.size, h), np.clip(E[k, m], 0, ZMAX3D), color="#0d366b", lw=1.4)
ax.set_xlim(RMIN3D, RMAX)
ax.set_ylim(H[0], H[-1])
ax.set_zlim(0, ZMAX3D)
ax.set_xlabel("R (Å)")
ax.set_ylabel("height h (Å)")
ax.set_zlabel("E − E(H$_2$ gas) (eV)")
ax.view_init(elev=24, azim=-128)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor("white")
    a.pane.set_edgecolor(GRID)
    a._axinfo["grid"]["color"] = GRID
ax.set_title("H$_2$ on Si(111): 3D PES", color=TEXT, fontsize=11, loc="left")
fig.text(0.02, 0.02, f"data used for the fit (fit/pes_data.csv): {src_txt} points\n{fill_txt[2:]}; "
         f"lines: computed heights; surface only where heights are ≤ 0.5 Å apart; colour scale "
         f"as in pes_2d.png", color=TEXT2, fontsize=7.5)
fig.savefig(os.path.join(TOP, "pes_3d.png"), dpi=200, bbox_inches="tight", facecolor="white")
print("wrote pes_2d.png, pes_3d.png")
