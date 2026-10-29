"""Side view (x-z) of both relaxations: start (scan slab, open circles) vs relaxed (filled),
bonds of the relaxed structure, dangling-bond atoms labelled, new/broken bonds marked.

Run in relax_Si111_fixbot/ (needs its sibling relax_SiH[_fixbot] too):
    python3 structure_change.py -> structure_change.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#b4b2aa"
DB = (1, 2, 13, 14, 9, 10, 23, 24)            # dangling-bond Si of the scan slab (1-based)


def read(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    n = sum(int(x) for x in L[6].split())
    k = 8 if L[7].strip()[0] in "sS" else 7
    P = np.array([[float(x) for x in l.split()[:3]] for l in L[k + 1:k + 1 + n]])
    return cell, (P @ cell if L[k].strip()[0] in "dD" else P * s)


def bonds(cell, P):
    out = {}
    L = np.diag(cell)
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            d = P[j] - P[i]
            d -= np.round(d / L) * L
            if np.linalg.norm(d) < (1.8 if max(i, j) >= 24 else 2.65):
                out[(i, j)] = d
    return out


fig, axs = plt.subplots(1, 2, figsize=(12, 7.2), sharey=True)
SIB = os.path.basename(HERE).replace("relax_Si111", "relax_SiH")
FIXED = "fixbot" in os.path.basename(HERE)
for ax, d, title in zip(axs, (os.path.basename(HERE), SIB), ("bare slab", "slab + H on Si 13")):
    cell, P0 = read(os.path.join(TOP, d, "POSCAR.orig"))
    _, P1 = read(os.path.join(TOP, d, "CONTCAR"))
    b0, b1 = bonds(cell, P0), bonds(cell, P1)
    for (i, j), v in b1.items():                       # relaxed bonds (minimum image)
        new = (i, j) not in b0
        ax.plot([P1[i, 0], P1[i, 0] + v[0]], [P1[i, 2], P1[i, 2] + v[2]],
                color=ORANGE if new else "#9ec5f4", lw=3 if new else 1.6, zorder=1,
                solid_capstyle="round")
    for (i, j), v in b0.items():                       # bonds that broke
        if (i, j) not in b1:
            ax.plot([P1[i, 0], P1[j, 0]], [P1[i, 2], P1[j, 2]], color=ORANGE, lw=1.6,
                    ls=(0, (2, 2)), zorder=1)
    for i in range(len(P0)):
        h = i >= 24
        ax.scatter(P0[i, 0], P0[i, 2], s=70 if h else 150, facecolors="none", edgecolors=GRAY,
                   linewidths=1.2, zorder=2)
        ax.annotate("", (P1[i, 0], P1[i, 2]), (P0[i, 0], P0[i, 2]),
                    arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.8), zorder=2)
        ax.scatter(P1[i, 0], P1[i, 2], s=70 if h else 150, color=GREEN if h else BLUE,
                   edgecolors="white", linewidths=0.8, zorder=3)
        if (i + 1) in DB or h:
            dz = P1[i, 2] - P0[i, 2]
            ax.text(P1[i, 0] + 0.22, P1[i, 2], f"{'H' if h else i + 1}  {dz:+.2f}", fontsize=8,
                    color=TEXT, va="center", zorder=4,
                    bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))
    ax.set_title(f"{title}", color=TEXT, fontsize=11, loc="left")
    ax.set_xlabel("x (Å)")
    ax.set_xlim(-0.6, cell[0, 0] + 1.3)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_aspect("equal")
axs[0].set_ylabel("z (Å)")
axs[0].set_ylim(-0.2, 10.2)
fig.text(0.01, 0.01, "open grey: start (scan slab); filled: relaxed; grey line: displacement; "
         "orange: new bond (solid) / broken bond (dashed); labels: dangling-bond Si and H with dz (Å)",
         color=TEXT2, fontsize=8)
fig.suptitle("What moved in the relaxations" + (" (bottom bilayer fixed)" if FIXED else " (all atoms free)")
             + "  -  side view x-z, bonds < 2.65 Å", x=0.01, ha="left",
             color=TEXT, fontsize=12)
fig.savefig(os.path.join(HERE, "structure_change.png"), dpi=170, bbox_inches="tight",
            facecolor="white")
print("wrote structure_change.png")
