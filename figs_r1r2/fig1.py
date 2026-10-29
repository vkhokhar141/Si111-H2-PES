"""Fig. 1 analogue: coordinate system r1, r2 for H2 upright above Si* (Si 13) of Si(111),
rendered with VMD (TachyonInternal, no display needed) and converted to PNG with Pillow.

r1 = H1-H2 distance (R in the PES), r2 = height of H1 above Si* (h in the PES).
Writes fig1_structure.xyz (slab with its periodic images around Si*, plus H1, H2) and
fig1.tcl, runs VMD, and saves fig1_coords.png (cropped).

Run in figs_r1r2/:  python3 fig1.py
"""
import os
import subprocess
import numpy as np
from PIL import Image, ImageChops

HERE = os.path.dirname(os.path.abspath(__file__))
POSCAR = os.path.join(HERE, "..", "1.50", "R_2.00", "POSCAR")   # any scan point: slab geometry
H_SHOW, R_SHOW = 1.50, 2.00      # A, geometry drawn: H1 at the Si-H minimum, H2 2 A above it
SI13 = 12                        # 0-based index of Si*
VMD = os.environ.get("VMD", "vmd")

# ---- slab + images within a window around Si*
lines = open(POSCAR).read().splitlines()
cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
si = np.array([[float(x) for x in l.split()[:3]] for l in lines[8:32]])
c = si[SI13]
atoms = []
for i in (-1, 0, 1):
    for j in (-1, 0, 1):
        for k, p in enumerate(si):
            q = p + i * cell[0] + j * cell[1]
            if abs(q[0] - c[0]) < 4.7 and abs(q[1] - c[1]) < 2.8:
                atoms.append(("Si", q, i == 0 and j == 0 and k == SI13))
h1 = np.array([c[0], c[1], c[2] + H_SHOW])
h2 = h1 + [0, 0, R_SHOW]
atoms += [("H", h1, False), ("H", h2, False)]
star = [n for n, a in enumerate(atoms) if a[2]][0]
with open(os.path.join(HERE, "fig1_structure.xyz"), "w") as f:
    f.write(f"{len(atoms)}\nSi(111) slab around Si* + H1 (h = {H_SHOW} A) + H2 (R = {R_SHOW} A)\n")
    for n, p, _ in atoms:
        f.write(f"{n} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
nH1, nH2 = len(atoms) - 2, len(atoms) - 1


def v(p):
    return "{" + " ".join(f"{x:.4f}" for x in p) + "}"


def dashed(p, q, r=0.06, n=9):
    """Dashed cylinder p -> q as Tcl draw commands."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    out = []
    for k in range(n):
        a = p + (q - p) * (k / n)
        b = p + (q - p) * ((k + 0.55) / n)
        out.append(f"draw cylinder {v(a)} {v(b)} radius {r} resolution 20")
    return out


def double_arrow(p, q, r=0.045, head=0.28, hr=0.13):
    """Dashed double-headed arrow p <-> q."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    u = (q - p) / np.linalg.norm(q - p)
    out = dashed(p + u * head, q - u * head, r=r, n=7)
    out.append(f"draw cone {v(p + u * head)} {v(p)} radius {hr} resolution 24")
    out.append(f"draw cone {v(q - u * head)} {v(q)} radius {hr} resolution 24")
    return out


ax = c[0] + 1.15                 # x of the arrows, right of the H2 axis
tcl = f"""color Display Background white
display projection Orthographic
display depthcue off
display shadows off
display ambientocclusion on
axes location off
color change rgb 0 0.27 0.45 0.77
color change rgb 7 0.13 0.69 0.30
mol new {{{os.path.join(HERE, 'fig1_structure.xyz')}}} type xyz waitfor all
set all [atomselect top all]
[atomselect top "name Si"] set radius 2.1
[atomselect top "name H"] set radius 1.0
mol bondsrecalc top
mol delrep 0 top
mol representation CPK 0.9 0.32 40 40
mol selection "name Si"
mol color ColorID 0
mol material AOChalky
mol addrep top
mol representation VDW 0.42 40
mol selection "index {nH1} {nH2}"
mol color ColorID 7
mol material AOChalky
mol addrep top
mol representation Bonds 0.16 40
mol selection "index {star} {nH1}"
mol color ColorID 0
mol material AOChalky
mol addrep top
# view: look along -y (z up), slightly turned for depth
rotate x by -90
rotate y by 10
scale by 1.45
translate by -0.05 -0.42 0
draw material AOChalky
draw color 7
{chr(10).join(dashed(h1 + [0, 0, 0.42], h2 - [0, 0, 0.42]))}
draw color black
{chr(10).join(double_arrow([ax, c[1], h1[2]], [ax, c[1], h2[2]]))}
{chr(10).join(double_arrow([ax, c[1], c[2]], [ax, c[1], h1[2]]))}
draw text {v([ax + 0.22, c[1] - 3, (h1[2] + h2[2]) / 2 - 0.15])} "r1" size 1.9 thickness 4
draw text {v([ax + 0.22, c[1] - 3, (c[2] + h1[2]) / 2 - 0.15])} "r2" size 1.9 thickness 4
draw text {v([h2[0] - 1.9, c[1] - 3, h2[2] - 0.15])} "H2" size 1.6 thickness 4
draw text {v([h1[0] - 1.9, c[1] - 3, h1[2] - 0.15])} "H1" size 1.6 thickness 4
draw text {v([c[0] - 2.3, c[1] - 3, c[2] + 0.35])} "Si*" size 1.6 thickness 4
render TachyonInternal {{{os.path.join(HERE, 'fig1_coords.tga')}}}
quit
"""
open(os.path.join(HERE, "fig1.tcl"), "w").write(tcl)
# empty startup file: skip ~/.vmdrc (its Tk GUI plugins fail in text mode). No 'display resize':
# it kills VMD 2.0 in text mode, so the render has the default size (672 x 848).
open(os.path.join(HERE, ".vmd_startup_empty.tcl"), "w").close()
subprocess.run([VMD, "-dispdev", "text", "-startup", os.path.join(HERE, ".vmd_startup_empty.tcl"),
                "-e", os.path.join(HERE, "fig1.tcl")],
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
im = Image.open(os.path.join(HERE, "fig1_coords.tga")).convert("RGB")
bbox = ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox()
pad = 20
im = im.crop((max(bbox[0] - pad, 0), max(bbox[1] - pad, 0),
              min(bbox[2] + pad, im.size[0]), min(bbox[3] + pad, im.size[1])))
im.save(os.path.join(HERE, "fig1_coords.png"))
os.remove(os.path.join(HERE, "fig1_coords.tga"))
print(f"wrote fig1_coords.png {im.size}, fig1_structure.xyz ({len(atoms)} atoms), fig1.tcl")
