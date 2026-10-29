"""Collect all converged scan points into one table for fitting.

Reads SCF_data/<h>/<h>.dat (E0, moments, steps) and Edisp from SCF_data/<h>/R_*/OUTCAR. Only points
marked ok whose OUTCAR reached EDIFF are used; heights with missing points are
listed, and their converged points are kept. Converged reruns in ../reruns, ../reruns2, ...
(in that order, so a later round overrides an earlier one) replace scan values according to
their 'use' rule in list.txt there (column src = scan / reruns / reruns2 ...), and
only if the slab's dangling bonds stayed non-magnetic (sum |m| over Si 1 2 14 9 10 23 24
< DB_MAX): the polarized-dangling-bond state is ~0.2 eV lower and must not be mixed in.
Points in EXCLUDE are dropped (listed with the reason when run). Rerun points merged into
<h>/R_<R>/ by ../merge_reruns.py are read from there (src from their source.txt), so the
rerun directories can be deleted. Also writes pes_r1r2_E0.csv (r1 = R, r2 = h, E0 only).

Run from fit/:  python3 collect_data.py   ->  pes_data.csv
"""
import os
import re
import numpy as np
import pandas as pd

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCF = os.path.join(TOP, "SCF_data")   # one directory per height
NPTS = 40
DB_ATOMS = (1, 2, 14, 9, 10, 23, 24)   # dangling-bond Si except 13, which bonds the H
DB_MAX = 0.15                          # muB, sum of |m| over DB_ATOMS
# points left out of the fit: no trustworthy non-magnetic energy after all rerun rounds
EXCLUDE = {
    (3.00, 2.60): "scan +26 meV off its neighbours; reruns2/3 polarized the dangling bonds"
                  " (rejected); the row is smooth (< 8 meV) without it",
}


def edisp(outcar):
    val = None
    with open(outcar) as f:
        for line in f:
            if "Edisp" in line:
                val = float(line.split()[-1])
    return val


def db_moment(outcar):
    """Sum of |m| on the dangling-bond Si atoms in the last magnetization (x) table."""
    m, block = {}, False
    with open(outcar) as f:
        for line in f:
            if "magnetization (x)" in line:
                m, block = {}, True
            elif block:
                p = line.split()
                if p and p[0].isdigit():
                    m[int(p[0])] = float(p[-1])
                elif p and p[0] == "tot":
                    block = False
    return sum(abs(m[i]) for i in DB_ATOMS)


def converged(outcar):
    if not os.path.isfile(outcar):
        return False
    with open(outcar) as f:
        return any("aborting loop because EDIFF is reached" in line for line in f)


rows = []
for d in sorted(os.listdir(SCF)):
    if not re.fullmatch(r"\d+\.\d\d", d):
        continue
    dat = os.path.join(SCF, d, d + ".dat")
    if not os.path.isfile(dat):
        continue
    lines = [l.split() for l in open(dat) if not l.startswith("#")]
    good = [l for l in lines if l[-1] == "ok"
            and converged(os.path.join(SCF, d, f"R_{l[0]}", "OUTCAR"))]
    if not good:
        continue
    if len(good) < NPTS:
        print(f"{d}: {len(good)} of {NPTS} points converged, using those")
    h = float(d)
    for l in good:
        R = float(l[0])
        ed = edisp(os.path.join(SCF, d, f"R_{l[0]}", "OUTCAR"))
        if ed is None:
            raise SystemExit(f"{d}/R_{l[0]}: no Edisp in OUTCAR")
        rows.append(dict(h=h, R=R, z_low=h, z_up=h + R, F=float(l[2]), E0=float(l[3]),
                         Edisp=ed, E_pbe=float(l[3]) - ed, mag=float(l[4]),
                         mH1=float(l[5]), mH2=float(l[6]), steps=int(l[7]), status=l[8]))

df = pd.DataFrame(rows)
# rerun points already merged into <h>/R_<R>/ (merge_reruns.py) name their round in source.txt
df["src"] = [open(p).readline().strip() if os.path.isfile(p) else "scan"
             for p in (os.path.join(SCF, f"{h:.2f}", f"R_{R:.2f}", "source.txt")
                       for h, R in zip(df.h, df.R))]

# reruns (../reruns, ../reruns2, ...; see list.txt there): 'lower' replaces the current
# value if the rerun is lower with the slab in the same state, 'always' replaces it, 'never'
# is a control; values are read from each rerun's own OSZICAR/OUTCAR (reruns.dat is only
# for reading)
RERUN_DIRS = sorted((x for x in os.listdir(TOP) if re.fullmatch(r"reruns\d*", x)),
                    key=lambda x: int(x[6:] or 1))


def h_moments(outcar):
    """m(H1), m(H2) (ions 25, 26) from the last magnetization (x) table."""
    m, block = {}, False
    with open(outcar) as f:
        for line in f:
            if "magnetization (x)" in line:
                m, block = {}, True
            elif block:
                p = line.split()
                if p and p[0].isdigit():
                    m[int(p[0])] = float(p[-1])
                elif p and p[0] == "tot":
                    block = False
    return m[25], m[26]


for rdir in RERUN_DIRS:
    LIST = os.path.join(TOP, rdir, "list.txt")
    if not os.path.isfile(LIST):
        continue
    for l in open(LIST):
        p = l.split()
        if l.startswith("#") or len(p) < 4:
            continue
        h, R, use = float(p[0]), float(p[1]), p[3]
        d = os.path.join(TOP, rdir, f"{p[0]}_R_{p[1]}")
        k = df.index[np.isclose(df.h, h) & np.isclose(df.R, R)]
        if len(k) != 1 or not converged(os.path.join(d, "OUTCAR")):
            continue
        k = k[0]
        with open(os.path.join(d, "OSZICAR")) as f:
            last = [x for x in f if " F= " in x][-1].split()
        with open(os.path.join(d, "OSZICAR")) as f:
            steps = sum(1 for x in f if x[:4] in ("DAV:", "RMM:", "CG :"))
        F, E0, mag = float(last[2]), float(last[4]), float(last[-1])
        mH1, mH2 = h_moments(os.path.join(d, "OUTCAR"))
        dE = 1000 * (E0 - df.at[k, "E0"])
        dms = (mag - (mH1 + mH2) / 0.381) - (df.at[k, "mag"] - (df.at[k, "mH1"] + df.at[k, "mH2"]) / 0.381)
        db = db_moment(os.path.join(d, "OUTCAR"))
        take = (use == "always" or (use == "lower" and dE < 0 and abs(dms) < 0.15)) and db < DB_MAX
        why = "used" if take else ("scan value kept (dangling bonds polarized)" if db >= DB_MAX
                                   else "scan value kept")
        print(f"{rdir} {p[0]}/R_{p[1]} ({use}): {dE:+.1f} meV, slab moment {dms:+.2f} muB vs scan,"
              f" DB moment {db:.2f} -> {why}")
        if take:
            ed = edisp(os.path.join(d, "OUTCAR"))
            df.loc[k, ["F", "E0", "Edisp", "E_pbe", "mag", "mH1", "mH2", "steps", "src"]] = \
                [F, E0, ed, E0 - ed, mag, mH1, mH2, steps, rdir]

for (h, R), why in EXCLUDE.items():
    k = df.index[np.isclose(df.h, h) & np.isclose(df.R, R)]
    if len(k):
        print(f"excluded {h:.2f}/R_{R:.2f}: {why}")
        df = df.drop(k)

df = df.sort_values(["h", "R"]).reset_index(drop=True)
df.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pes_data.csv"),
          index=False, float_format="%.8f")
# plain table: r1 = R (H1-H2), r2 = h (H1 above Si 13), both A; E0 in eV (PBE+D3, absolute)
df.rename(columns={"R": "r1", "h": "r2"})[["r1", "r2", "E0"]].to_csv(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "pes_r1r2_E0.csv"),
    index=False, float_format="%.8f")
hs = sorted(df.h.unique())
print(f"{len(df)} points, {len(hs)} heights: " + " ".join(f"{h:.2f}" for h in hs))
