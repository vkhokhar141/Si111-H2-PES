"""Plot the H-H scan at h = 6.00 A: E0(R) relative to R = 8 A, and the H spin moments.

Run from 6.00:  python3 plot_pes.py   ->  pes.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#8d8c86"
YMAX = 2.0   # eV; the wall at R <= 0.35 A is drawn as off-scale markers

rows = [l.split() for l in open("6.00.dat") if not l.startswith("#")]
d = np.array([[float(r[0]), float(r[3]), float(r[5]), float(r[6])] for r in rows])
d = d[np.argsort(d[:, 0])]
R, E0, mH1, mH2 = d.T
Einf = E0[R == 8.0][0]
dE = E0 - Einf

# old spin-polarized run (6.00s, ENCUT 330, default MAGMOM) for comparison
old = []
for l in open("../../old/6.00/6.00s.dat"):
    p = l.split()
    if "E0=" in p:
        old.append((float(p[1]), float(p[p.index("E0=") + 1])))
old = np.array(sorted(old))
old_dE = old[:, 1] - old[old[:, 0] == 8.0, 1][0]

# equilibrium from a quartic fit around the minimum
m = (R >= 0.65) & (R <= 0.85)
c = np.polyfit(R[m], E0[m], 4)
rr = np.linspace(0.65, 0.85, 20001)
ee = np.polyval(c, rr)
re_, De = rr[np.argmin(ee)], Einf - ee.min()

plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.linewidth": 0.8})
fig, (ax, bx) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True,
                             gridspec_kw={"height_ratios": [2.3, 1], "hspace": 0.08})
for a in (ax, bx):
    a.grid(True, color=GRID, lw=0.6)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)

# --- energy panel
ax.axhline(0, color=TEXT2, lw=0.8, ls=":")
ax.plot(old[:, 0], old_dE, color=GRAY, lw=1.5, ls="--", marker="o", ms=4, mfc="white",
        label="old ISPIN=2 run (6.00s): both H spins parallel")
ax.plot(R, dE, color=BLUE, lw=2, marker="o", ms=5, mec="white", mew=0.8,
        label="new scan: opposite H spins")
off = dE > YMAX
note = "off scale: " + ", ".join(f"R = {r:.2f} Å → {e:+.1f} eV" for r, e in zip(R[off], dE[off]))
ax.text(0.62, YMAX - 0.12, note, color=TEXT2, fontsize=8, va="top")
ax.annotate(f"$r_e$ = {re_:.3f} Å\n$D_e$ = {De:.2f} eV", xy=(re_, -De), xytext=(1.25, -4.2),
            color=TEXT, fontsize=9, arrowprops=dict(arrowstyle="->", color=TEXT2, lw=0.8))
ax.set_ylim(-5.0, YMAX)
ax.set_ylabel("E0(R) − E0(8 Å)  (eV)")
ax.set_title("H$_2$ 6 Å above Si(111): H–H stretch (PBE-D3BJ, ENCUT 400)",
             color=TEXT, fontsize=11, loc="left")
ax.legend(frameon=False, loc="lower right", fontsize=8.5, labelcolor=TEXT)

# --- spin panel
bx.axhline(0, color=TEXT2, lw=0.8)
bx.plot(R, mH1, color=BLUE, lw=2, marker="o", ms=4, mec="white", mew=0.6)
bx.plot(R, mH2, color=ORANGE, lw=2, marker="o", ms=4, mec="white", mew=0.6)
bx.text(8.05, mH1[-1], "lower H", color=TEXT, va="center", fontsize=8.5)
bx.text(8.05, mH2[-1], "upper H", color=TEXT, va="center", fontsize=8.5)
bx.set_ylim(-0.5, 0.5)
bx.set_ylabel("H moment (μB)")
bx.set_xlabel("H–H distance R (Å)")
bx.set_xlim(0, 8.7)

fig.savefig("pes.png", dpi=200, bbox_inches="tight", facecolor="white")
print(f"r_e = {re_:.3f} A, De = {De:.3f} eV; wrote pes.png")
