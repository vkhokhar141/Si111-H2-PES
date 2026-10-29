"""Data files for the r1/r2 figures (gnuplot), from the final fitted PES and the fit data.

r1 = H1-H2 distance (= R), r2 = height of H1 above Si* (Si 13) (= h). Energies in eV.
Units of r1, r2 set by UNIT ('bohr' like the reference figure, or 'A').

  pes_grid.dat      r1 r2 V      fitted PES (PBE+D3) relative to gas-phase H2, gnuplot grid
                                 blocks (one per r2), down to small r1 and r2 so the repulsive
                                 walls show; values as fitted (no cut), plot top VTOP above max
  cut_chem.dat      r1 V         H1 chemisorbed: r2 at the Si-H minimum, V relative to r1 -> inf
  cut_away.dat      r1 V         H1 away from the surface (r2 = 8 A): gas-phase H2 curve, same zero
  pts_chem.dat, pts_away.dat     the DFT points on those two rows (h = 1.50 A and 8 A)
  vh_sudden.dat     r2 V         H1 above Si* with H2 far (R = 8 A), frozen slab, relative to H1 far
  vh_pts.dat        r2 V         DFT points of that column
  vh_relaxed.dat    (not written here) a relaxed-slab H1 curve, if one is computed later,
                    is plotted in panel (c) when this file exists: columns r2 V, same zero
  labels.gp         gnuplot variables (unit label, ranges, Si-H minimum)

Run in figs_r1r2/:  python3 make_data.py
"""
import os
import sys
import numpy as np
import pandas as pd
import joblib

FIT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fit")
sys.path.insert(0, FIT)
UNIT = "bohr"                        # 'bohr' or 'A'
BOHR = 0.529177210903                # A
f = 1 / BOHR if UNIT == "bohr" else 1.0
lab = "bohr" if UNIT == "bohr" else "Å"

pes = joblib.load(os.path.join(FIT, "pes_model.joblib"))
df = pd.read_csv(os.path.join(FIT, "pes_data.csv"))
df["V"] = df.E_pbe - pes.e_ref + df.Edisp - pes.d_ref          # relative to the atoms (PBE+D3)


def save(name, cols, header):
    np.savetxt(name, np.column_stack(cols), fmt="%.6f", header=header)


# (a) 3D grid down to small r1, r2 (inside the DFT data: smallest computed R = 0.25 A,
#     h = 1.00 A) so the repulsive walls show; values written as they are (no cut): the plot's
#     top (VTOP) is set above the largest value on the grid
R1MIN = 1.0                                  # bohr; 0.5 reaches the smallest computed R but the
                                             # wall then goes to 27 eV and flattens the rest
r1 = np.linspace(R1MIN, 12.0, 111) / f
r2 = np.linspace(1.9, 12.0, 102) / f
RR, HH = np.meshgrid(r1, r2)                 # rows: r2 (gnuplot blocks), columns: r1
V = pes.energy(RR, HH, ref="gas")
VTOP = float(np.ceil(V.max()))
with open("pes_grid.dat", "w") as fh:
    fh.write(f"# r1({lab}) r2({lab}) V(eV, relative to gas-phase H2, PBE+D3), max {V.max():.2f}\n")
    for i, h in enumerate(r2):
        for a, b in zip(r1, V[i]):
            fh.write(f"{a * f:.5f} {h * f:.5f} {b:.6f}\n")
        fh.write("\n")

# (b) V(r1) with H1 chemisorbed (r2 = Si-H minimum) and H1 away (r2 = 8 A), each relative to r1 -> inf
hs = np.linspace(1.3, 1.8, 50001)
h_sih = hs[np.argmin(pes.energy(8.0, hs))]
rr = np.linspace(0.25, 8.0, 800)          # from the smallest computed R, so the walls show
for name, h in (("chem", h_sih), ("away", 8.0)):
    v = pes.energy(rr, h) - pes.energy(8.0, h)
    save(f"cut_{name}.dat", [rr * f, v], f"r1({lab}) V(eV) at r2 = {h:.4f} A, relative to r1 -> inf")
for name, h in (("chem", 1.50), ("away", 8.0)):
    g = df[np.isclose(df.h, h)].sort_values("R")
    ref = float(g[np.isclose(g.R, 8.0)].V.iloc[0])
    save(f"pts_{name}.dat", [g.R * f, g.V - ref], f"r1({lab}) V(eV) DFT at h = {h} A, relative to R = 8")

# (c) H1 above Si* with H2 far away (R = 8 A): frozen slab ("sudden"), relative to H1 far (h = 8)
hh = np.linspace(0.75, 8.0, 800)
v = pes.energy(8.0, hh) - pes.energy(8.0, 8.0)
save("vh_sudden.dat", [hh * f, v], f"r2({lab}) V(eV) frozen slab, R = 8 A, relative to h = 8 A")
g = df[np.isclose(df.R, 8.0)].sort_values("h")
save("vh_pts.dat", [g.h * f, g.V - float(g[np.isclose(g.h, 8.0)].V.iloc[0])], f"r2({lab}) V(eV) DFT")

with open("labels.gp", "w") as fh:
    fh.write(f'unit = "{lab}"\nr2_sih = {h_sih * f:.4f}\nvtop = {VTOP}\nr1min = {R1MIN}\n')
    fh.write(f'fit_note = "fitted PES, PBE+D3 ({len(df)} DFT points)"\n')
print(f"wrote data in {lab}: Si-H minimum r2 = {h_sih * f:.3f} {lab} ({h_sih:.4f} A)")
