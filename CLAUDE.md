# Si(111) + H2 potential energy surface (VASP)

2D PES E(R, h) for an upright H2 above one Si atom of a frozen Si(111) slab, built from
single-point VASP SCF calculations (PBE-D3), fitted with a GPR model for dynamics and MEP work.

## Geometry and coordinates
- Slab: 24 Si, frozen; cell 7.680 x 6.651 x 40.406 A. Top surface reconstructed (atoms 1, 2, 13,
  14 are 3-coordinated); bottom layer (atoms 9, 10, 23, 24) is **not** H-passivated -> metallic,
  open-shell slab. This is the main source of SCF difficulty (see "Slab spin states").
- H2 is collinear above **Si atom 13** (z = 7.944237044 A), lower H = ion 25, upper H = ion 26.
- **h** = height of the lower H above Si 13 = the height directory name (`1.35/` -> 1.35 A).
- **R** = H-H distance = the point directory (`R_0.75`). Upper H at z = 7.944237044 + h + R.
- Paper-style names (figs_r1r2): r1 = R (H1-H2), r2 = h (H1 above Si* = Si 13).
- Energy references: two far atoms = E(h = 8, R = 8); gas-phase H2 = E(h = 8, R = 0.75).
- Units: A and eV everywhere, except `figs_r1r2/` (bohr, like the reference paper).

## Status (2026-10-28)
- All 29 heights x 40 R computed; 3 rerun rounds done and merged. Fit data: 1159 points =
  984 scan + 170 reruns + 4 reruns2 + 1 reruns3; 3.00/R_2.60 excluded. QC: no energy outliers,
  per-height leave-one-out rms <= 3.1 meV. Fit, MEP analysis and all figures are current.
- Relaxations done: `relax_SiH_fixbot/`, `relax_Si111_fixbot/` (bottom bilayer fixed):
  relaxed Si-H binding 3.159 eV vs 3.254 frozen. All-atoms-free runs were unusable (the
  unpassivated bottom reconstructs) and were deleted (2026-10-29).
- Not done: a relaxed-slab H1-Si* curve (the "relaxed" curve of figs_r1r2 panel c): relax the
  slab at a series of fixed H heights (bottom bilayer fixed).
- Key results (PBE+D3, vs gas-phase H2): H + H-Si +1.277 eV (PBE +1.245) at h = 1.50 A;
  Si-H binding 3.25 eV frozen / 3.16 eV relaxed, 2022 cm-1 (Si frozen); physisorbed H2
  -52 meV at h = 2.85 A (mostly D3); MEP has no barrier beyond 10 meV (dissociation costs only
  the endothermicity; Eley-Rideal abstraction H + H-Si -> H2 is barrierless in this collinear
  geometry).

## Layout
- `<h>/` (29 heights: 1.00-1.30 step 0.1, 1.30-1.70 step 0.05, 1.70-2.20 step 0.1,
  2.50-8.00 step 0.5): `INCAR KPOINTS POTCAR scan.sh submit`, results in `<h>.dat`, `scan.log`,
  one `R_<R>/` per point (40 R values, 8.00 -> 0.25).
- `submit_all`: submits every unfinished height as a SLURM chain (`--dependency=afterany`).
- `reruns/`, `reruns2/`, `reruns3/`: recomputed single points, one round each (see "Reruns").
  The 175 used ones are merged into `<h>/R_<R>/` (`merge_reruns.py`); the rerun directories
  were then deleted locally (2026-10-29).
- `fit/`: fitting pipeline and MEP analysis (see "Fitting pipeline").
- `plot_pes2d3d.py`: data maps `pes_2d.png`, `pes_3d.png` from `fit/pes_data.csv` (the fit
  data; the excluded point is filled along R for drawing only).
- `figs_r1r2/`: paper-style figures in r1/r2 (see "Figures").
- `report/`: `PES_report.docx` (5-page user report: method, data, fit, files, usage) built by
  `python3 make_report.py` (python-docx; uses the figures in place); `PES_report.pdf` = render.
- `relax_SiH_fixbot/`, `relax_Si111_fixbot/`: relaxations of slab + H on Si 13 and of the bare
  slab with the bottom bilayer fixed (see "Relaxations").
- `6.00/plot_pes.py`, `1.00/plot_h100.py`: single-height figures. They exist **only locally**:
  copying height directories back from the cluster replaced them once (restored from the
  session transcript), so re-check after every copy. The cluster's `6.00/` may still have the
  old `scan.sh`/`submit` and `6.00new.dat` (locally fixed: `6.00.dat`, current scripts).
- `old/`: the original, broken spin-polarized scans (kept for comparison only).

## Running scans (on the cluster)
- `cd <h> && bash submit [JOBID]` or `bash submit_all` from the top. 32 MPI ranks,
  `OMP_NUM_THREADS=1`, 4 h limit; restartable (converged points are skipped).
- `bash <h>/scan.sh --check` -> exit 0 if all 40 points converged.
- Regenerate a stale `<h>.dat` without running VASP: `VASP_CMD=false bash <h>/scan.sh`.
- Each scan runs **from R = 8 inward**, restarting every point from the previous point's
  WAVECAR/CHGCAR. Do not reverse this: starting from the closed-shell 0.75 A point and going
  outward keeps the restricted (non-magnetic) solution past the H2 spin-symmetry break.
- After a copy from the cluster: wait until no file was written for a minute, then check.

## Electronic settings that matter (INCAR) and why
- `ISPIN = 2`, `MAGMOM = 24*0 1 -1` (slab non-magnetic, H atoms opposite = broken-symmetry
  singlet). The old default MAGMOM (1 per atom) put H2 in the repulsive triplet state and caused
  250-750-step SCF runs with magnetic sloshing.
- `AMIX_MAG = 0.8`, `BMIX_MAG = 0.0001`; default charge mixing is fine.
- `ENCUT = 400`, `PREC = Normal`, `IVDW = 12` (D3-BJ), `ISMEAR = 0`, `SIGMA = 0.2`; use **E0**.
- `IBRION = -1`, `NSW = 0` (single points; never relax: the upper H is the scan coordinate).
- `KPAR = 16`, `NCORE = 2` (16 irreducible k-points). `LORBIT = 11` for the OUTCAR per-atom
  magnetization table (it also writes PROCAR, which nothing reads).
- `NUPDOWN` (reruns only) fixes the total moment. **`NUPDOWN = -1` is VASP's "no constraint"**;
  for a total of -1 use NUPDOWN = 1 with flipped H spins (`MAGMOM = 24*0 -1 1`), same energy.
  NUPDOWN = 0 does **not** prevent the antiferromagnetic dangling-bond state (total 0).
- D3 is a non-self-consistent additive term: `Edisp` in OUTCAR, PBE energy = E0 - Edisp.
- INCAR comment lines must not contain `=` (VASP parses them as tags); use `!` comments.

## Slab spin states (why the reruns exist)
- The slab's real ground state has spin-polarized dangling bonds (top Si 1, 2, 13, 14 ~ -0.15,
  bottom 9, 10, 23, 24 ~ +0.16 muB per sphere, total 0, smearing entropy -0.10 vs -0.26 eV),
  ~0.2 eV below the non-magnetic slab that almost all scan points are in. The non-magnetic state
  is metastable; SCF runs drift partly toward polarization at random (sum |m| on the DB atoms
  0.3-0.6 -> 10-30 meV lower): this caused 10-30 meV row offsets, and ~20 single points sat on
  a 6-42 meV higher SCF solution.
- **Decision (2026-10-28, kept after finding the polarized ground state): non-magnetic slab
  everywhere**, as one consistent (metastable) reference. Test per point: sum |m| over Si
  1 2 14 9 10 23 24 < 0.15 muB (Si 13 excluded: it bonds the H). Total moment: 0 when H2 is
  intact or both H are free (h >= 3), -1 when the lower H is bonded and the upper free (h <= 2.5,
  large R). Rows h <= 2.2 were already non-magnetic; rows h = 3.0-4.0 keep scan values (their
  reruns polarized) and carry 10-20 meV offsets.
- No CHGCARs survive (local or cluster), so states can't be inherited from old points.

## Reruns (`reruns/`, `reruns2/`, `reruns3/`)
- Each round: `list.txt` (h, R, NUPDOWN, use rule, why), one `<h>_R_<R>/` per point, `run.sh [H]`
  (chains grid-neighbour points of a height with the same NUPDOWN like scan.sh, otherwise
  fresh; `--check [H]`, `--heights`), `submit H [JOBID]` (job `H2_rr_<h>`, `H2_rr2_<h>`, ...
  from the directory name; log `run_<h>.log`), `submit_all [--parallel]` (skips complete
  heights and heights whose job is in `squeue`). `reruns.dat` is a readable summary only.
- Use rules: `always` replaces the scan value, `lower` only if lower with the slab in the same
  state, `never` = control. `fit/collect_data.py` reads the rounds in order (later overrides)
  from each point's own OSZICAR/OUTCAR, and only if the dangling bonds stayed non-magnetic.
- New round: make `reruns4/` the same way; never reuse a finished round (run.sh skips converged
  points). To redo a point, delete its outputs here **and on the cluster**.
- **Merged (2026-10-29):** `python3 merge_reruns.py [--do]` copies every rerun point the fit
  uses (src != scan in `fit/pes_data.csv`) into `<h>/R_<R>/` (scan outputs kept in
  `R_<R>/scan_old/`, round named in `R_<R>/source.txt`) and rebuilds `<h>.dat`; collect_data.py
  then gives the same 1159 points without any rerun directory. Done locally only: copying
  height directories back from the cluster would undo it (rerun merge_reruns.py there, or copy
  the local height directories to the cluster). For a new round, run collect_data.py, then
  merge_reruns.py --do.
- History: round 1 (297 points): 6 controls within 0.7 meV of the scan; 19 outliers fixed
  (6-42 meV lower); 171 merged, 119 rejected (polarized, h = 3.0/3.5/4.0); QC flags 61 -> 15.
  Round 2 (5): 1.00/R_2.10, 1.55/R_2.20, 3.00/R_2.40 used (-16 to -21 meV); 2.50/R_8.00 equal
  to the scan; 2.50/R_7.00 kept. Round 3 (3): 2.50/R_7.00 -11 meV used; 3.00/R_2.60, R_2.80
  polarized -> rejected; 3.00/R_2.60 then excluded (`EXCLUDE` in collect_data.py: scan +26 meV
  off, no non-magnetic rerun; the row is smooth < 8 meV without it). Locally, reruns/2.50_R_7.00,
  2.50_R_8.00, 3.00_R_2.40 have no outputs (still on the cluster; not needed).

## Fitting pipeline (`fit/`)
Run in order after new data arrive (fit ~20 min, the rest ~1 min):
`python3 collect_data.py && python3 qc.py > qc_report.txt && python3 fit_pes.py && python3 plot_fit.py && python3 analyze.py && (cd .. && python3 plot_pes2d3d.py && bash figs_r1r2/make_figs.sh)`
- `collect_data.py` -> `pes_data.csv` (h, R, E0, Edisp, E_pbe = E0 - Edisp, moments, steps,
  src = scan / reruns / reruns2 ...): converged points, rerun merge, `EXCLUDE` list. Also
  `pes_r1r2_E0.csv`: just r1 = R, r2 = h (A) and E0 (eV, absolute) for the same 1159 points.
- `qc.py`: leave-one-out checks on the slab moment and on E(h, R) - E(8, R); flags, never drops.
- `pes_model.py`: `PES` class. V = V_HH(R) + V_H(h) + V_H(h+R) + taper * GPR(3-body)
  [+ GPR(D3)]. V_HH = h = 8 row, V_H = R = 8 column (cubic splines, 3-parameter exponential
  wall below the grid). 3-body GPR (Matern 5/2 + White, `normalize_y=False` so the far field
  is exactly 0) is tapered to 0 over 5.5-8 A in R and h, so both asymptotic channels equal
  the 1D cuts exactly. `energy`, `gradient` (analytic), `sigma`, `energy_z`, `energy_com`.
- `fit_pes.py`: leave-height-out (1.35, 1.45, 1.55, 1.65, 2.10, 3.50) and leave-every-other-R
  CV, picks the GP input transform, gradient FD check, gas-phase check, writes
  `pes_model.joblib`, `fit_report.txt`. Scales ~n^3.
  Last fit: transform 'exp'; CV RMSE (E < 3 eV) 1.0-3.4 meV leave-height-out (14 meV at
  h = 3.50, the slab-state offsets) and 3.4 meV (max 14) every-other-R; training RMSE 0.9 meV,
  max 6.2 meV below 3 eV; gas H2 r_e 0.7509 A, De 4.531 eV, omega_e 4276 cm-1.
- `plot_fit.py`: `fit_pes_2d.png`, `fit_cuts.png`, `fit_pes_3d.png`.
- `analyze.py`: MEP H2(g) -> H + H-Si by the simplified string method (241 images) in the H
  Cartesian heights (z1 = h, z2 = h + R: equal masses, so the mass-weighted frame; never in
  (R, h)), end points fixed at (R = r_e, h = 8) and (R = 8, h = h_SiH); stationary points refined
  by Newton on the analytic gradient + Hessian frequencies; extrema less prominent than
  NOISE = 10 meV are only counted (the flat valleys ripple by 1-5 meV). PBE+D3 and PBE.
  Writes `analyze_report.txt`, `mep_2d.png`, `mep_profile.png`, `sigma_map.png`,
  `residuals_map.png`. The PBE "physisorption well" (-17 meV at h = 4.0) sits on a row with a
  slab-state offset: not meaningful.

## Figures (`figs_r1r2/`, bohr)
- `bash make_figs.sh`: `make_data.py` (gnuplot data from the fit; `UNIT` bohr/A), `fig2.gp`
  (gnuplot), `fig1.py` (VMD).
- `fig2a_pes3d.png`: mesh + base contours from r1 = R1MIN = 1.0, r2 = 1.9 bohr, view 42/36.
  Values as fitted, **never cut or edited**: the z top VTOP is set above the grid maximum
  (9 eV; R1MIN = 0.5 reaches the smallest computed R but needs 27 eV and flattens the rest).
- `fig2b_cuts.png`: V(r1) with H1 chemisorbed (r2 at the Si-H minimum) and H1 away, each
  relative to r1 -> inf, with the DFT points. `fig2c_vh.png`: H1-Si* frozen-slab ("sudden")
  curve; a "relaxed" curve is added if `vh_relaxed.dat` (r2 V, zero at H1 far) exists.
  `fig2.png`: all three, laid out like the paper.
- `fig1_coords.png`: coordinate system rendered with VMD. VMD 2.0 text mode: run with an
  empty `-startup` file (the user's ~/.vmdrc loads a Tk plugin that fails), and never use
  `display resize` (kills VMD in text mode; default render 672x848). TGA -> PNG with Pillow.
- gnuplot: old-style `if (c) cmd` swallows the rest of the line (fatal inside one-line macros);
  use `if (c) { cmd }`.

## Relaxations (`relax_SiH_fixbot/`, `relax_Si111_fixbot/`)
- Two directories, bottom bilayer fixed: slab + one H on Si 13 (`relax_SiH_fixbot`, 25 atoms)
  and bare slab (`relax_Si111_fixbot`, 24 atoms, Si-only POTCAR = first dataset of the Si+H
  POTCAR). Never relax with all atoms free (bottom artifact, see below).
- Settings: PES settings, ISPIN = 1 (non-magnetic slab, as for the PES; no free H, so nothing
  needs spin; slab + H has 97 electrons, fine with smearing, as in the PES where its moment was
  ~0), IBRION = 2, EDIFFG = -0.01, EDIFF = 1E-6, cell fixed, no dipole correction (as in the
  PES). Start = scan slab (+ H 1.50 A above Si 13), kept as `POSCAR.orig`. POSCARs:
  `Selective dynamics`, Si 3 4 5 6 9 10 23 24 `F F F` (z 0.90 and 1.59 A; next layer at
  4.02 A), rest `T T T`.
- Scripts: `bash submit` (jobs SiH_relax_fb, Si111_relax_fb; 4 h; log relax.log); `run.sh`
  continues from CONTCAR/WAVECAR after a time limit (old outputs -> `run_<n>/`; its CONTCAR
  check counts a Selective dynamics line); `run.sh --check`. `summary.py` (same file in both,
  stdlib only) -> summary: relaxation energy E0(final) - E0(first ionic step = frozen
  scan geometry, as in the PES), largest force on free atoms (fixed atoms keep theirs), atom
  moves, Si-H bond/height/tilt. In `relax_Si111_fixbot`: `binding.py` (pairs with
  `relax_SiH_fixbot` by name): E(H) = [E0(h = 8, R = 8 scan) - E0(bare frozen)]/2, Si-H binding
  E_b = E(slab) + E(H) - E(slab + H), frozen and relaxed; `structure_change.py` ->
  `structure_change.png` (side view, start vs relaxed). Copies from the cluster may lack
  OSZICAR: summary.py and binding.py then read the energies from OUTCAR.
- Results, bottom bilayer fixed (2026-10-28; SiH 27 steps, bare 64; free-atom forces
  < 0.010 eV/A): no bonds formed or broken, bonds 2.32-2.42 A. Bare slab -0.132 eV: Si 13 (Si*)
  sinks 0.60 A and Si 14 0.23 A (empty-dangling-bond flattening, like the paper's C* puckering),
  Si 1, 2 up 0.06 A. Slab + H -0.037 eV: top bilayer ~0.1 A inward, H follows, Si-H 1.501 A
  unchanged; H holds Si* 0.50 A further out than on the relaxed bare surface (z 7.85 vs 7.35).
  Si-H binding: frozen 3.254 eV (= PES exactly, so ISPIN = 1 is consistent with the PES),
  relaxed 3.159 eV (-0.095 eV = the bare slab's extra relaxation).
- Results, all atoms free (old `relax_SiH/`, `relax_Si111/`, deleted 2026-10-29; do not
  repeat): the unpassivated bottom reconstructs. Bare slab
  -0.947 eV (Si 9, 10 down 0.33, Si 23, 24 up 0.93 A; bonds 23-11, 24-12 form, 3-11, 4-12
  break); slab + H -0.276 eV (bottom buckles, no bond switching). "Relaxed binding" 2.583 eV
  is dominated by this artifact.

## Disk space
Each completed height writes ~4 GB (40 CHGCAR + WAVECAR/CHG). PROCAR (2 GB), vasprun.xml
(3 GB) and DOSCAR (0.8 GB) over all heights are never read either. Once a height is complete
(`scan.sh --check`), CHGCAR, WAVECAR and CHG can be deleted; energies/moments live in
OUTCAR/OSZICAR/`<h>.dat`. Before deleting in a freshly copied tree, confirm the copy has
finished (no files written in the last minute). For an incomplete height, keep the last
converged point's CHGCAR/WAVECAR so the next point can restart.

## Conventions
- Plots: PNG only (no PDF). Sequential single-hue blue colour scale for energies.
- Never interpolate or draw the surface across height gaps > ~1 A; leave them blank.
- Never cut, clip or edit PES values to make a figure look better: choose ranges and axis
  limits so the data fit.
- Compare energies before any refit when new reruns arrive; refit only when asked.
