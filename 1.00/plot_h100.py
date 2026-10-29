"""Analysis figure for the h = 1.00 A scan (lower H 1.00 A above Si 13).

(a) E0(R) vs gas-phase H2, compared with h = 8.00
(b) h = 1.00 relative to its R = 8 A limit, split into DFT and D3 dispersion
(c) spin moments of the two H atoms and the cell total

Run from 1.00:  python3 plot_h100.py   ->  pes_1.00.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#8d8c86"


def load(dat):
    rows = [l.split() for l in open(dat) if not l.startswith("#")]
    a = np.array([[float(r[0]), float(r[3]), float(r[4]), float(r[5]), float(r[6])] for r in rows])
    return a[np.argsort(a[:, 0])]


def edisp(R):
    lines = [l for l in open(f"R_{R:.2f}/OUTCAR") if "Edisp" in l]
    return float(lines[-1].split()[-1])


d1 = load("1.00.dat")
d8 = load("../8.00/8.00.dat")
R, E0, mag, mH1, mH2 = d1.T
Eref = d8[np.argmin(abs(d8[:, 0] - 0.75)), 1]          # gas-phase H2 at equilibrium
D = np.array([edisp(r) for r in R])
i8 = np.argmin(abs(R - 8.0))
dE, dD = E0 - E0[i8], D - D[i8]
dDFT = dE - dD

plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.linewidth": 0.8})
fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(7.2, 9.6), sharex=True,
                                 gridspec_kw={"height_ratios": [1.5, 1.2, 1], "hspace": 0.12})
for a in (ax, bx, cx):
    a.grid(True, color=GRID, lw=0.6)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
mk = dict(marker="o", ms=4.5, mec="white", mew=0.7)

# (a) vs gas-phase H2
YMAX = 12.0
ax.axhline(0, color=TEXT2, lw=0.8)
ax.plot(d8[:, 0], d8[:, 1] - Eref, color=GRAY, lw=1.6, **mk, label="h = 8.00 Å (gas-phase H$_2$)")
ax.plot(R, E0 - Eref, color=BLUE, lw=2, **mk, label="h = 1.00 Å")
ax.axhline(E0[i8] - Eref, color=BLUE, lw=0.8, alpha=0.5)
ax.annotate("", xy=(6.6, E0[i8] - Eref), xytext=(6.6, 0),
            arrowprops=dict(arrowstyle="<->", color=TEXT2, lw=0.9))
ax.text(6.7, (E0[i8] - Eref) / 2, f"{E0[i8] - Eref:.2f} eV\nH on Si (1.00 Å)\n+ free H",
        color=TEXT, fontsize=8.5, va="center")
i75 = np.argmin(abs(R - 0.75))
ax.annotate(f"+{E0[i75] - Eref:.2f} eV at R = 0.75 Å\n(no H$_2$ well)", xy=(0.75, E0[i75] - Eref),
            xytext=(1.5, 9.3), color=TEXT, fontsize=8.5,
            arrowprops=dict(arrowstyle="->", color=TEXT2, lw=0.8))
ax.annotate("gas H$_2$ well, −4.53 eV", xy=(0.75, d8[np.argmin(abs(d8[:, 0] - 0.75)), 1] - Eref),
            xytext=(1.5, -3.6), color=TEXT, fontsize=8.5,
            arrowprops=dict(arrowstyle="->", color=TEXT2, lw=0.8))
off = (E0 - Eref) > YMAX
ax.text(0.03, 0.97, "off scale: " + ", ".join(f"R = {r:.2f} Å → +{e:.1f} eV"
        for r, e in zip(R[off], E0[off] - Eref)), transform=ax.transAxes,
        color=TEXT2, fontsize=7.5, va="top")
ax.set_ylim(-5.2, YMAX)
ax.set_ylabel("E0 − E(H$_2$ gas)  (eV)")
ax.set_title("(a)  Lower H held 1.00 Å above Si: H–H stretch", color=TEXT, fontsize=10.5, loc="left")
ax.legend(frameon=False, loc="lower right", fontsize=8.5, labelcolor=TEXT)

# (b) DFT vs D3, zoomed
z = R >= 1.45
bx.axhline(0, color=TEXT2, lw=0.8)
bx.axvspan(1.72, 2.45, color="#f1f0ec", zorder=0)
bx.text(2.08, 285, "~10 meV scatter:\nslab spin shifting", color=TEXT2, fontsize=7.5, ha="center", va="top")
bx.plot(R[z], 1000 * dE[z], color=BLUE, lw=2, **mk, label="total E0")
bx.plot(R[z], 1000 * dDFT[z], color=ORANGE, lw=1.6, **mk, label="DFT part")
bx.plot(R[z], 1000 * dD[z], color=AQUA, lw=1.6, **mk, label="D3 dispersion")
j = np.argmin(dE)
bx.annotate(f"{1000 * dE[j]:.0f} meV dip at {R[j]:.2f} Å\n= D3 attraction", xy=(R[j], 1000 * dE[j]),
            xytext=(3.4, 120), color=TEXT, fontsize=8.5,
            arrowprops=dict(arrowstyle="->", color=TEXT2, lw=0.8))
bx.set_ylim(-80, 300)
bx.set_ylabel("E − E(R = 8 Å)  (meV)")
bx.set_title("(b)  h = 1.00 Å relative to its R = 8 Å limit (R ≥ 1.5 Å)", color=TEXT,
             fontsize=10.5, loc="left")
bx.legend(frameon=False, loc="upper right", fontsize=8.5, labelcolor=TEXT)

# (c) spin
cx.axhline(0, color=TEXT2, lw=0.8)
cx.plot(R, mag, color=GRAY, lw=1.6, **mk, label="cell total")
cx.plot(R, mH1, color=BLUE, lw=2, **mk, label="lower H (on Si)")
cx.plot(R, mH2, color=ORANGE, lw=2, **mk, label="upper H")
cx.set_ylim(-1.15, 0.15)
cx.set_ylabel("moment (μB)")
cx.set_xlabel("H–H distance R (Å)")
cx.set_title("(c)  spin moments", color=TEXT, fontsize=10.5, loc="left")
cx.legend(frameon=False, loc="center right", fontsize=8.5, labelcolor=TEXT)
cx.set_xlim(0.2, 8.2)

fig.savefig("pes_1.00.png", dpi=200, bbox_inches="tight", facecolor="white")
print("wrote pes_1.00.png")
