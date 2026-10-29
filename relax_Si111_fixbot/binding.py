"""Si-H binding energy on Si 13, frozen and relaxed (standard library only).

Uses the two relaxations and one scan point:
  this directory   bare slab:            first ionic step = frozen scan slab, last = relaxed
  its sibling      slab + H on Si 13:    first ionic step = frozen slab, H 1.50 A, last = relaxed
                   (relax_Si111_fixbot <-> relax_SiH_fixbot)
  ../8.00/8.00.dat R = 8 point: frozen slab + two free H atoms (non-magnetic slab)
E(H)  = [E0(h = 8, R = 8) - E0(bare slab, frozen)] / 2
E_b   = E0(slab) + E(H) - E0(slab + H)      (> 0: bound), frozen and relaxed
Check: the frozen E_b must reproduce the PES value -V_H(1.50 A) = 3.254 eV (PBE+D3), which
also confirms that the ISPIN = 1 relaxations match the ISPIN = 2 PES (non-magnetic slab).

Run in relax_Si111_fixbot/ (after both relaxations):  python3 binding.py
"""
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
SIH = os.path.join(TOP, os.path.basename(HERE).replace("relax_Si111", "relax_SiH"))
PES_EB = 3.254           # eV, -V_H at the frozen-slab minimum, fit/analyze_report.txt (PBE+D3)


def steps(d):
    runs = sorted(glob.glob(os.path.join(d, "run_*", "OSZICAR")),
                  key=lambda p: int(os.path.basename(os.path.dirname(p)).split("_")[1]))
    files = runs + ([os.path.join(d, "OSZICAR")] if os.path.isfile(os.path.join(d, "OSZICAR")) else [])
    e0 = [float(l.split()[4]) for f in files for l in open(f) if " F=" in l]
    outcar = os.path.join(d, "OUTCAR")
    if not e0 and os.path.isfile(outcar):          # no OSZICAR copied: E0 per ionic step from OUTCAR
        e0 = [float(l.split()[-1]) for l in open(outcar) if "energy(sigma->0)" in l
              and "energy  without entropy" in l]
    done = os.path.isfile(outcar) and "reached required accuracy" in open(outcar).read()
    return e0, done


slab, slab_ok = steps(HERE)
sih, sih_ok = steps(SIH)
if not slab or not sih:
    raise SystemExit(f"need at least one ionic step in both {HERE} and {SIH}")
e88 = [float(l.split()[3]) for l in open(os.path.join(TOP, "8.00", "8.00.dat"))
       if not l.startswith("#") and l.split()[0] == "8.00"][0]
eH = (e88 - slab[0]) / 2

print(f"{os.path.basename(HERE)} + {os.path.basename(SIH)}")
print(f"E0 bare slab   frozen {slab[0]:.6f}   final {slab[-1]:.6f} eV"
      f"   ({'relaxed' if slab_ok else 'NOT converged'}, relaxation {slab[-1] - slab[0]:+.4f} eV)")
print(f"E0 slab + H    frozen {sih[0]:.6f}   final {sih[-1]:.6f} eV"
      f"   ({'relaxed' if sih_ok else 'NOT converged'}, relaxation {sih[-1] - sih[0]:+.4f} eV)")
print(f"E(H) from the h = 8, R = 8 scan point = {eH:.6f} eV")
eb0 = slab[0] + eH - sih[0]
ebr = slab[-1] + eH - sih[-1]
print(f"Si-H binding, frozen slab : {eb0:.3f} eV   (PES fit {PES_EB:.3f} eV; difference"
      f" {1000 * (eb0 - PES_EB):+.0f} meV)")
print(f"Si-H binding, relaxed     : {ebr:.3f} eV   (relaxation changes it by {ebr - eb0:+.3f} eV)")
if not (slab_ok and sih_ok):
    print("note: a relaxation has not converged yet; the 'relaxed' numbers are provisional")
