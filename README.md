# 2D PES for H2 on Si(111) (VASP PBE-D3)

Potential energy surface V(r1, r2) for an H2 molecule standing upright above one surface Si
atom (Si*) of a frozen Si(111) slab, from 1159 VASP single points (PBE + D3), with a fitted
surface (physical baseline + Gaussian-process regression, analytic gradients) for dynamics
and reaction-path work.

- r1 = H1-H2 distance, r2 = height of the lower atom H1 above Si* (both in Angstrom; called
  R and h in the scripts).
- **Main data: [`fit/pes_r1r2_E0.csv`](fit/pes_r1r2_E0.csv)**: columns r1, r2 (Angstrom), E0
  (eV, absolute VASP energy, PBE+D3). Zero references: free atoms E(8, 8) = -129.37379 eV,
  gas-phase H2 E(0.75, 8) = -133.90415 eV.
- Full table with the D3 part, PBE-only energy and spin moments: `fit/pes_data.csv`.
- **Read first: [`report/PES_report.pdf`](report/PES_report.pdf)** (5 pages: method, data, fit,
  files, usage; `.docx` alongside). Detailed notes, settings and pitfalls: [`CLAUDE.md`](CLAUDE.md).

## Using the fitted PES

```python
import joblib, pes_model                     # run in fit/ (pes_model.py must be importable)
pes = joblib.load("pes_model.joblib")
E = pes.energy(r1, r2, d3=True, ref="gas")   # eV vs gas-phase H2; ref="atoms": vs free atoms
g = pes.gradient(r1, r2)                     # (dV/dr1, dV/dr2), eV/A
s = pes.sigma(r1, r2)                        # GP uncertainty, eV
```

Made with Python 3.10, numpy 2.2, scipy 1.15, scikit-learn 1.7, pandas 2.2, matplotlib 3.10,
joblib 1.5. The pickled model needs a compatible scikit-learn; otherwise refit with
`python3 fit_pes.py` (about 20 min) from `fit/pes_data.csv`.

## Layout

| Path | Content |
|---|---|
| `SCF_data/<r2>/` (29 heights, e.g. `SCF_data/1.50/`) | one scan: `INCAR KPOINTS scan.sh`, table `<r2>.dat`, one `R_<r1>/` per point |
| `SCF_data/<r2>/R_<r1>/` | `INCAR KPOINTS POSCAR OSZICAR` of each point; `source.txt` if the point is a rerun |
| `fit/` | data collection, QC, fit, MEP analysis, the CSV files and the fitted model |
| `figs_r1r2/` | paper-style figures (bohr) |
| `relax_SiH_fixbot/`, `relax_Si111_fixbot/` | slab relaxations (bottom bilayer fixed), Si-H binding |
| `report/` | the report and its generator `make_report.py` |
| `submit.example` | sample SLURM script to run a height scan (`scan.sh` needs `VASP_CMD`) |
| `old/` | the first, discarded spin-polarized scan (kept for comparison) |

## Not in this repository

- Cluster job scripts (`submit`, `submit_all`); use `submit.example` as a template.
- **POTCAR files** (VASP license). To rerun a point, add your own PAW_PBE `Si 05Jan2001` +
  `H 15Jun2001` POTCAR (in that order).
- Large VASP outputs (OUTCAR, CONTCAR, XDATCAR, vasp.out, ...) of the 1159 single points.
  Their energies, D3 terms and moments are in `fit/pes_data.csv`, so the fit and all figures
  can be rebuilt from the repository, but `fit/collect_data.py` (which reads the OUTCARs) can
  only be rerun after recomputing the points.
