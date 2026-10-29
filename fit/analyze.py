"""Minimum-energy path and energetics on the fitted 2D PES (pes_model.joblib).

Reaction (collinear, H2 upright above Si 13):
    H2(g)  ->  H + H-Si          (dissociation with one H bonded to Si 13, the other free)
and its reverse, Eley-Rideal abstraction  H(g) + H-Si -> H2(g) + Si.

MEP: simplified string method (E, Ren, Vanden-Eijnden 2007) in the Cartesian heights of the
two H atoms, z1 = h (lower H), z2 = h + R (upper H). Both atoms are H, so this is the
mass-weighted frame and steepest descent there is the physical MEP (in (R, h) the metric is
skewed). End points are fixed in the two asymptotic valleys: H2(g) at (R = r_e, h = 8) and
H + H-Si at (R = 8, h = h_SiH). The highest image is refined to a stationary point (Newton on
the analytic gradient) and classified by its Hessian. Everything is done for PBE+D3 and PBE.

Output: analyze_report.txt, mep_2d.png, mep_profile.png, sigma_map.png, residuals_map.png
Run after fit_pes.py:  python3 analyze.py
"""
import numpy as np
import pandas as pd
import joblib
from scipy.optimize import minimize, root, minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, TwoSlopeNorm, LogNorm
from pes_model import H_INF, R_INF

NIMG = 241            # images on the string
DT = 0.004            # A^2/eV, string time step
STEP_MAX = 0.01       # A, largest move of an image per iteration
TOL = 1e-7            # A, convergence: largest change of the string per iteration
MAXIT = 40000
M_H = 1.00794         # amu
NOISE = 0.010         # eV; stationary points less prominent than this are data/fit ripples
TEXT, TEXT2, GRID, DOT = "#0b0b0b", "#52514e", "#e4e3df", "#0b0b0b"
BLUE, ORANGE, RED, MID = "#2a78d6", "#eb6834", "#e34948", "#f0efec"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("blue_seq", BLUES)
DIVERGE = LinearSegmentedColormap.from_list("blue_gray_red", ["#184f95", BLUE, MID, RED, "#a32a2a"])

pes = joblib.load("pes_model.joblib")
log = []


def out(s=""):
    print(s, flush=True)
    log.append(s)


# ---------------- PES in the Cartesian frame z = (z1, z2) = (h, h + R)
def V(z, d3):
    z = np.atleast_2d(z)
    return pes.energy(z[:, 1] - z[:, 0], z[:, 0], d3=d3, ref="gas")


def G(z, d3):
    z = np.atleast_2d(z)
    g = pes.gradient(z[:, 1] - z[:, 0], z[:, 0], d3=d3)          # (dV/dR, dV/dh)
    return np.column_stack([g[:, 1] - g[:, 0], g[:, 0]])          # (dV/dz1, dV/dz2)


def hessian(z, d3, eps=1e-4):
    H = np.empty((2, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = eps
        H[:, j] = (G(z + e, d3)[0] - G(z - e, d3)[0]) / (2 * eps)
    return 0.5 * (H + H.T)


def wavenumber(lam):
    """eV/A^2 eigenvalue of the Cartesian Hessian (H mass) -> cm^-1 (negative = imaginary)."""
    w = np.sqrt(abs(lam) * 1.602177e-19 / 1e-20 / (M_H * 1.66054e-27)) / (2 * np.pi * 2.99792458e10)
    return np.sign(lam) * w


def to_z(R, h):
    return np.column_stack([np.ravel(h), np.ravel(h) + np.ravel(R)])


def reparam(z):
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(z, axis=0), axis=1))])
    t = np.linspace(0, s[-1], len(z))
    return np.column_stack([np.interp(t, s, z[:, 0]), np.interp(t, s, z[:, 1])])


def string_method(z0, d3):
    z = reparam(z0)
    for it in range(MAXIT):
        step = -DT * G(z[1:-1], d3)
        n = np.linalg.norm(step, axis=1, keepdims=True)
        step *= np.minimum(1, STEP_MAX / np.maximum(n, 1e-30))
        znew = z.copy()
        znew[1:-1] += step
        znew = reparam(znew)
        dz = np.max(np.abs(znew - z))
        z = znew
        if dz < TOL:
            break
    return z, it + 1, dz


# ---------------- reference points
def gas_re(d3):
    r = minimize_scalar(lambda x: float(pes.energy(x, H_INF, d3=d3)), bounds=(0.6, 0.9),
                        method="bounded", options={"xatol": 1e-7})
    return r.x, -r.fun


def sih_min(d3):
    """Minimum of one H against the surface (the other H far away, R = 8)."""
    r = minimize_scalar(lambda x: float(pes.energy(R_INF, x, d3=d3)), bounds=(1.2, 2.0),
                        method="bounded", options={"xatol": 1e-7})
    e = lambda x: float(pes.energy(R_INF, x, d3=d3))
    k = (e(r.x + 1e-3) - 2 * e(r.x) + e(r.x - 1e-3)) / 1e-6
    return r.x, r.fun, wavenumber(k)                              # Si frozen: mass m_H


res = {}
for d3, lab in ((True, "PBE+D3"), (False, "PBE")):
    re_, De = gas_re(d3)
    hS, eS_atoms, wS = sih_min(d3)
    zA, zB = to_z(re_, H_INF)[0], to_z(R_INF, hS)[0]
    # initial string: down the H2 valley at R = r_e, then out along the Si-H valley
    corner = to_z(re_, hS)[0]
    n1 = NIMG // 2
    z0 = np.vstack([np.linspace(zA, corner, n1, endpoint=False), np.linspace(corner, zB, NIMG - n1)])
    z, nit, dz = string_method(z0, d3)
    E = V(z, d3)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(z, axis=0), axis=1))])
    R, h = z[:, 1] - z[:, 0], z[:, 0]
    sig = pes.sigma(R, h)

    # stationary points along the path: interior local maxima and minima of E(s), keeping
    # only those more prominent than NOISE (adjacent min/max pairs closer than NOISE, or an
    # extremum within NOISE of an end point, are removed repeatedly, smallest first)
    ext = [0] + [i for i in range(1, NIMG - 1) if (E[i] > E[i - 1] and E[i] >= E[i + 1])
                 or (E[i] < E[i - 1] and E[i] <= E[i + 1])] + [NIMG - 1]
    n_all = len(ext) - 2
    while True:
        gaps = [(abs(E[ext[k + 1]] - E[ext[k]]), k) for k in range(len(ext) - 1)]
        gaps = [(g, k) for g, k in gaps if g < NOISE]
        if not gaps or len(ext) == 2:
            break
        g, k = min(gaps)
        if k == 0:
            del ext[1]                     # extremum next to the start point
        elif k + 1 == len(ext) - 1:
            del ext[k]                     # extremum next to the end point
        else:
            del ext[k:k + 2]               # an adjacent min/max pair
    stat = []
    for i in ext[1:-1]:
        sol = root(lambda x: G(x, d3)[0], z[i], jac=lambda x: hessian(x, d3), method="hybr")
        zs = sol.x
        lam, vec = np.linalg.eigh(hessian(zs, d3))
        kind = {0: "minimum", 1: "saddle"}.get(int(np.sum(lam < 0)), "maximum")
        stat.append(dict(i=i, z=zs, R=zs[1] - zs[0], h=zs[0], E=float(V(zs, d3)[0]),
                         gnorm=float(np.linalg.norm(G(zs, d3))), lam=lam, kind=kind,
                         w=wavenumber(lam), sigma=float(pes.sigma(zs[1] - zs[0], zs[0])),
                         ok=sol.success))
    res[lab] = dict(re=re_, De=De, hS=hS, eS=eS_atoms, wS=wS, z=z, E=E, s=s, R=R, h=h, sig=sig,
                    stat=stat, nit=nit, dz=dz, Eprod=float(V(zB, d3)[0]), n_ripple=n_all - len(stat))
    # physisorption well of intact H2: local minimum started from the lowest image of the path
    ph = minimize(lambda x: float(V(x, d3)[0]), z[int(np.argmin(E))], jac=lambda x: G(x, d3)[0],
                  method="BFGS", options={"gtol": 1e-8})
    res[lab]["phys"] = dict(R=ph.x[1] - ph.x[0], h=ph.x[0], E=ph.fun)

# ---------------- report
out("MEP and energetics on the fitted PES (pes_model.joblib); energies in eV relative to")
out("gas-phase H2 at r_e unless stated. Coordinates: R = H-H, h = lower H above Si 13 (A).")
out("Caveat: rows h = 3.0-4.0 keep scan values with 10-20 meV slab-state offsets, so features")
out("< ~20 meV there (e.g. the PBE physisorption well) are within the data uncertainty.")
for lab, r in res.items():
    out(f"\n== {lab}")
    out(f"  string: {NIMG} images, {r['nit']} iterations, last change {r['dz']:.1e} A,"
        f" path length {r['s'][-1]:.2f} A")
    out(f"  H2(g)            r_e = {r['re']:.4f} A, De = {r['De']:.3f} eV")
    out(f"  Si-H (H far)     h = {r['hS']:.3f} A, binding to Si 13 = {-r['eS']:.3f} eV"
        f" (vs free H), omega = {r['wS']:.0f} cm-1 (Si frozen)")
    out(f"  H + H-Si         E = {r['Eprod']:+.3f} eV  (dissociation endothermic by this much;"
        f" abstraction H + H-Si -> H2 exothermic by it)")
    p = r["phys"]
    out(f"  H2 physisorbed   R = {p['R']:.3f} A, h = {p['h']:.3f} A, E = {1000 * p['E']:+.1f} meV")
    for st in r["stat"]:
        w = ", ".join(f"{x:.0f}" if x >= 0 else f"{-x:.0f}i" for x in st["w"])
        out(f"  {st['kind']:8s} on MEP  R = {st['R']:.3f} A, h = {st['h']:.3f} A, E = {st['E']:+.3f} eV,"
            f" |grad| = {st['gnorm']:.1e} eV/A, frequencies {w} cm-1, GP sigma {1000 * st['sigma']:.1f} meV")
    out(f"  ({r['n_ripple']} more stationary points along the path are ripples < {1000 * NOISE:.0f} meV"
        f" deep, in the flat asymptotic valleys: data/fit noise, not reported)")
    sad = [st for st in r["stat"] if st["kind"] == "saddle"]
    emax, imax = r["E"].max(), int(np.argmax(r["E"]))
    if sad:
        b = max(sad, key=lambda x: x["E"])
        out(f"  barrier: dissociation {b['E']:+.3f} eV (above H2 gas), abstraction"
            f" {b['E'] - r['Eprod']:+.3f} eV (above H + H-Si)")
    else:
        out(f"  no barrier beyond {1000 * NOISE:.0f} meV: from the physisorption well E rises monotonically"
            f" to H + H-Si (highest image {emax:+.3f} eV at R = {r['R'][imax]:.2f} A, h ="
            f" {r['h'][imax]:.2f} A; product {r['Eprod']:+.3f} eV).")
        out(f"  -> dissociation needs only the endothermicity {r['Eprod']:+.3f} eV; Eley-Rideal abstraction"
            f" H + H-Si -> H2 is barrierless here (within {1000 * NOISE:.0f} meV).")
    out(f"  GP sigma along the path: max {1000 * r['sig'].max():.1f} meV")
open("analyze_report.txt", "w").write("\n".join(log) + "\n")

# ---------------- plots
plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.linewidth": 0.8})
df = pd.read_csv("pes_data.csv")
df["E"] = df.E_pbe - pes.e_ref + df.Edisp - pes.d_ref - pes.gas_ref(True)
df["fit"] = pes.energy(df.R.values, df.h.values, ref="gas")
main, alt = res["PBE+D3"], res["PBE"]


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def mark_stationary(ax, r, xkey="R", ykey="h"):
    for st in r["stat"]:
        if st["kind"] == "saddle":
            ax.plot(st[xkey], st[ykey], marker="X", ms=10, color=ORANGE, mec="white", mew=1.2, zorder=6)


# 1. MEP on the surface
RG, HG = np.meshgrid(np.linspace(0.4, 8, 400), np.linspace(1.0, 8, 400))
EG = pes.energy(RG, HG, ref="gas")
levels = np.arange(-0.25, 3.0 + 1e-9, 0.25)
norm = BoundaryNorm(levels, CMAP.N, extend="max")
fig, ax = plt.subplots(figsize=(7.8, 6.6))
cf = ax.contourf(RG, HG, EG, levels=levels, cmap=CMAP, norm=norm, extend="max")
ax.contour(RG, HG, EG, levels=np.arange(0.5, 3.0 + 1e-9, 0.5), colors="white", linewidths=0.7)
ax.plot(main["R"], main["h"], color=ORANGE, lw=2.4, label="MEP, PBE+D3", zorder=5)
ax.plot(alt["R"], alt["h"], color=TEXT, lw=1.3, ls=(0, (4, 3)), label="MEP, PBE", zorder=5)
mark_stationary(ax, main)
for st in main["stat"]:
    if st["kind"] == "saddle":
        ax.annotate(f"saddle {st['E']:+.2f} eV", (st["R"], st["h"]), xytext=(14, 10),
                    textcoords="offset points", color=TEXT, fontsize=9,
                    bbox=dict(fc="white", ec="none", pad=1.5))
ax.text(main["re"] + 0.12, 7.6, "H$_2$(g)", color=TEXT, fontsize=9,
        bbox=dict(fc="white", ec="none", pad=1.5))
ax.text(7.9, main["hS"] + 0.15, f"H + H–Si  {main['Eprod']:+.2f} eV", color=TEXT, fontsize=9,
        ha="right", bbox=dict(fc="white", ec="none", pad=1.5))
ax.set_xlim(0.4, 8)
ax.set_ylim(1.0, 8)
ax.set_xlabel("H–H distance R (Å)")
ax.set_ylabel("height of lower H above Si h (Å)")
ax.set_title("Minimum-energy path on the fitted PES", color=TEXT, fontsize=11, loc="left")
style(ax)
ax.legend(frameon=True, facecolor="white", edgecolor="none", loc="upper right", fontsize=8.5,
          labelcolor=TEXT)
cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=np.arange(0, 3.0 + 1e-9, 0.5))
cb.set_label("E − E(H$_2$ gas), PBE+D3  (eV)", color=TEXT)
cb.outline.set_visible(False)
fig.text(0.01, 0.005, "string method in the H Cartesian heights (mass-weighted frame); colours"
         " capped at 3 eV", color=TEXT2, fontsize=7.5)
fig.savefig("mep_2d.png", dpi=200, bbox_inches="tight", facecolor="white")

# 2. energy profile along the path (left: whole path, right: the reaction zone), with R and h
#    below (same unit, one axis)
z0, z1 = main["s"][np.argmax(main["R"] > 0.8)] - 1.0, main["s"][np.argmax(main["h"] < main["hS"] + 0.02)] + 1.0
fig, axs = plt.subplots(2, 2, figsize=(11, 6.6), sharex="col",
                        gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.08, wspace=0.18))
for col, (lo, hi) in enumerate(((None, None), (z0, z1))):
    a1, a2 = axs[0, col], axs[1, col]
    for r, lab, c, ls in ((main, "PBE+D3", BLUE, "-"), (alt, "PBE", ORANGE, (0, (4, 3)))):
        a1.plot(r["s"], r["E"], color=c, lw=2, ls=ls, label=lab)
    for st in main["stat"]:
        x = main["s"][st["i"]]
        if lo is not None and not lo <= x <= hi:
            continue
        a1.plot(x, st["E"], marker="X" if st["kind"] == "saddle" else "o", ms=8, color=TEXT,
                mec="white", mew=1.2, zorder=5)
        a1.annotate(f"physisorbed H$_2$ (PBE+D3)  {1000 * st['E']:+.0f} meV\nR {st['R']:.2f}, h {st['h']:.2f} Å"
                    if st["kind"] == "minimum" else f"saddle {st['E']:+.3f} eV",
                    (x, st["E"]), xytext=(-150, 22), textcoords="offset points", fontsize=8,
                    color=TEXT, arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.6))
    a1.axhline(0, color=TEXT2, lw=0.6)
    a1.grid(True, color=GRID, lw=0.5)
    a1.set_axisbelow(True)
    a2.plot(main["s"], main["R"], color=TEXT, lw=1.6)
    a2.plot(main["s"], main["h"], color=TEXT, lw=1.6, ls=(0, (4, 3)))
    a2.grid(True, color=GRID, lw=0.5)
    a2.set_axisbelow(True)
    if lo is not None:
        a1.set_xlim(lo, hi)
        a2.set_ylim(0, 3.5)
        a1.set_title("reaction zone (zoom)", color=TEXT2, fontsize=9.5, loc="left")
    else:
        a1.set_title("whole path", color=TEXT2, fontsize=9.5, loc="left")
    for ax in (a1, a2):
        style(ax)
    a2.set_xlabel("distance along the path s (Å)")
a = axs[0, 0]
a.set_ylabel("E − E(H$_2$ gas) (eV)")
a.legend(frameon=False, loc="upper left", fontsize=8.5, labelcolor=TEXT)
a.axvspan(z0, z1, color=GRID, alpha=0.6, lw=0)
fig.suptitle("Energy along the minimum-energy path  (H$_2$(g) → H + H–Si)", color=TEXT,
             fontsize=11, x=0.125, ha="left", y=0.97)
a.text(main["s"][-1], 1.13, f"H + H–Si\nPBE+D3 {main['Eprod']:+.3f} eV\nPBE {alt['Eprod']:+.3f} eV",
       color=TEXT, fontsize=8, ha="right", va="top")
b = axs[1, 0]
b.set_ylabel("Å")
b.text(main["s"][6], main["R"][6] + 0.5, "R (H–H)", fontsize=8.5, color=TEXT)
b.text(main["s"][6], main["h"][6] - 1.3, "h (lower H above Si)", fontsize=8.5, color=TEXT)
fig.text(0.01, -0.02, f"s: arc length in the Cartesian heights of the two H atoms (mass-weighted frame)  ·  "
         f"GP 1σ along the path ≤ {1000 * main['sig'].max():.0f} meV (too narrow to draw)  ·  no barrier "
         f"beyond {1000 * NOISE:.0f} meV in either functional", color=TEXT2, fontsize=7.5)
fig.savefig("mep_profile.png", dpi=200, bbox_inches="tight", facecolor="white")

# 3. GP uncertainty map (3-body part, PBE)
RS, HS = np.meshgrid(np.linspace(0.25, 8, 300), np.linspace(1.0, 8, 300))
SG = 1000 * pes.sigma(RS, HS)
fig, ax = plt.subplots(figsize=(7.8, 6.6))
lv = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50])
pm = ax.contourf(RS, HS, np.clip(SG, lv[0], None), levels=lv, cmap=CMAP, norm=LogNorm(lv[0], lv[-1]),
                 extend="max")
ax.scatter(df.R, df.h, s=4, c="white", edgecolors=DOT, linewidths=0.3, zorder=4, label="DFT points")
ax.plot(main["R"], main["h"], color=ORANGE, lw=2, label="MEP, PBE+D3", zorder=5)
ax.set_xlim(0.25, 8)
ax.set_ylim(1.0, 8)
ax.set_xlabel("H–H distance R (Å)")
ax.set_ylabel("height of lower H above Si h (Å)")
ax.set_title("GP uncertainty of the fit (3-body part, 1σ)", color=TEXT, fontsize=11, loc="left")
style(ax)
ax.legend(frameon=True, facecolor="white", edgecolor="none", loc="upper right", fontsize=8.5,
          labelcolor=TEXT)
cb = fig.colorbar(pm, ax=ax, pad=0.02, ticks=lv)
cb.ax.set_yticklabels([f"{x:g}" for x in lv])
cb.set_label("σ (meV)", color=TEXT)
cb.outline.set_visible(False)
fig.text(0.01, 0.005, "σ is 0 where the 3-body term is tapered off (R or h > 5.5–8 Å): there the"
         " surface is the 1D cuts themselves", color=TEXT2, fontsize=7.5)
fig.savefig("sigma_map.png", dpi=200, bbox_inches="tight", facecolor="white")

# 4. residuals fit - DFT at every data point
d = df[df.R <= 8]
resid = 1000 * (d.fit - d.E)
fig, ax = plt.subplots(figsize=(7.8, 6.6))
lim = 5.0
sc = ax.scatter(d.R, d.h, c=np.clip(resid, -lim, lim), cmap=DIVERGE, norm=TwoSlopeNorm(0, -lim, lim),
                s=np.where(d.E < 3, 26, 10), edgecolors="none", zorder=3)
ax.set_xlim(0.15, 8.1)
ax.set_ylim(0.9, 8.1)
ax.set_xlabel("H–H distance R (Å)")
ax.set_ylabel("height of lower H above Si h (Å)")
m3 = d.E < 3
ax.set_title(f"Fit − DFT at every data point  (E < 3 eV: RMS {np.sqrt(np.mean(resid[m3] ** 2)):.1f} meV,"
             f" max {np.abs(resid[m3]).max():.1f} meV)", color=TEXT, fontsize=11, loc="left")
style(ax)
cb = fig.colorbar(sc, ax=ax, pad=0.02, extend="both")
cb.set_label("fit − DFT (meV), clipped at ±5", color=TEXT)
cb.outline.set_visible(False)
fig.text(0.01, 0.005, "large dots: E < 3 eV above H$_2$ gas; small dots: higher (wall) points",
         color=TEXT2, fontsize=7.5)
fig.savefig("residuals_map.png", dpi=200, bbox_inches="tight", facecolor="white")
print("wrote analyze_report.txt, mep_2d.png, mep_profile.png, sigma_map.png, residuals_map.png")
