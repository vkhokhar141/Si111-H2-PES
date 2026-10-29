"""Summary of a slab relaxation (relax_SiH_fixbot/ or relax_Si111_fixbot/; standard library only).

Final E0/F and steps, largest remaining force, relaxation energy (final E0 minus the first
ionic step, which is the frozen scan slab geometry, as in the PES), and how far the atoms
moved from POSCAR.orig. With an H on Si 13 (25 atoms): Si 13-H bond, H height above Si 13
and tilt. The same file is used in both directories.

Run in the relaxation directory:  python3 summary.py
"""
import glob
import math
import os

SI13, H = 13, 25                     # 1-based ion numbers
TOP_DB, BOT_DB = (1, 2, 13, 14), (9, 10, 23, 24)
PES_SIH = 1.501                      # A, Si-H height from the frozen-slab fit (fit/analyze.py)


def read_poscar(path):
    lines = open(path).read().splitlines()
    scale = float(lines[1].split()[0])
    cell = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    counts = [int(x) for x in lines[6].split()]
    k = 7
    if lines[k].strip()[0] in "sS":          # selective dynamics
        k += 1
    direct = lines[k].strip()[0] in "dD"
    n = sum(counts)
    pos = []
    for line in lines[k + 1:k + 1 + n]:
        v = [float(x) for x in line.split()[:3]]
        if direct:
            v = [sum(v[j] * cell[j][i] for j in range(3)) for i in range(3)]
        else:
            v = [scale * x for x in v]
        pos.append(v)
    return cell, pos


def mic(d, cell):
    """Minimum-image vector for an orthorhombic cell."""
    return [d[i] - cell[i][i] * round(d[i] / cell[i][i]) for i in range(3)]


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def energies(oszicar):
    out = []
    for line in open(oszicar):
        if " F=" in line:
            p = line.split()
            out.append((float(p[2]), float(p[4])))     # F, E0
    return out


def last_forces(outcar):
    f, block, cur = None, False, []
    for line in open(outcar):
        if "TOTAL-FORCE" in line:
            block, cur = True, []
            continue
        if block:
            p = line.split()
            if len(p) == 6:
                cur.append([float(x) for x in p[3:]])
            elif cur and line.strip().startswith("---"):
                block, f = False, cur
    return f


runs = sorted(glob.glob("run_*/OSZICAR"), key=lambda p: int(p.split("_")[1].split("/")[0]))
steps = [e for r in runs for e in energies(r)] + (energies("OSZICAR") if os.path.isfile("OSZICAR") else [])
if not steps:
    raise SystemExit("no ionic steps yet")
relaxed = os.path.isfile("OUTCAR") and "reached required accuracy" in open("OUTCAR").read()
cell, p0 = read_poscar("POSCAR.orig")
_, p1 = read_poscar("CONTCAR" if os.path.getsize("CONTCAR") > 0 else "POSCAR")

print(f"status: {'relaxed (EDIFFG reached)' if relaxed else 'NOT converged yet'};"
      f" {len(steps)} ionic steps over {len(runs) + 1} run(s)")
print(f"E0 first step (frozen scan slab geometry)     = {steps[0][1]:.6f} eV")
print(f"E0 final                                      = {steps[-1][1]:.6f} eV   (F = {steps[-1][0]:.6f})")
print(f"relaxation energy E0(final) - E0(first)       = {steps[-1][1] - steps[0][1]:+.4f} eV")
f = last_forces("OUTCAR") if os.path.isfile("OUTCAR") else None
if f:
    # selective dynamics: fixed atoms keep their forces, so only the free ones count
    pl = open("POSCAR").read().splitlines()
    sel = pl[7].strip()[:1] in "sS"
    free = [not (sel and "F" in l.split()[3:6]) for l in pl[9:9 + len(f)]]
    fm = [norm(v) for v in f]
    i = max((k for k in range(len(fm)) if free[k]), key=fm.__getitem__)
    print(f"largest remaining force on free atoms {fm[i]:.4f} eV/A on ion {i + 1}"
          + (f" ({free.count(False)} atoms fixed)" if not all(free) else ""))

if len(p0) >= H:
    for lab, p in (("start", p0), ("final", p1)):
        b = mic([p[H - 1][k] - p[SI13 - 1][k] for k in range(3)], cell)
        tilt = math.degrees(math.acos(max(-1, min(1, b[2] / norm(b)))))
        print(f"{lab}: Si13-H = {norm(b):.4f} A, H height above Si13 = {b[2]:.4f} A, tilt {tilt:.2f} deg")
    print(f"(frozen-slab PES fit: Si-H height {PES_SIH:.3f} A)")

disp = [norm(mic([p1[i][k] - p0[i][k] for k in range(3)], cell)) for i in range(len(p0))]
dz = [p1[i][2] - p0[i][2] for i in range(len(p0))]
print("dangling-bond Si, displacement |dr| (dz) in A:")
for lab, ids in (("top", TOP_DB), ("bottom", BOT_DB)):
    print(f"  {lab:6s} " + "  ".join(f"{i}: {disp[i - 1]:.3f} ({dz[i - 1]:+.3f})" for i in ids))
order = sorted(range(len(disp)), key=lambda i: -disp[i])[:6]
print("largest displacements: " + "  ".join(f"ion {i + 1}: {disp[i]:.3f}" for i in order))
