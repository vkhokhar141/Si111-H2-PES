"""Build PES_report.docx: short report on the H2/Si(111) PES data, method, fit and files.

Run in report/:  python3 make_report.py   (needs python-docx; figures from ../fit, ../figs_r1r2)
"""
import os
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
ACCENT = RGBColor(0x1F, 0x4E, 0x8C)

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Inches(8.5), Inches(11)
for side in ("left_margin", "right_margin"):
    setattr(sec, side, Inches(0.9))
sec.top_margin = sec.bottom_margin = Inches(0.8)

st = doc.styles["Normal"]
st.font.name = "Calibri"
st.font.size = Pt(10.5)
st.paragraph_format.space_after = Pt(4)
st.paragraph_format.line_spacing = 1.08
for name, size in (("Heading 1", 13.5), ("Heading 2", 11.5)):
    h = doc.styles[name]
    h.font.name, h.font.size, h.font.bold = "Calibri", Pt(size), True
    h.font.color.rgb = ACCENT
    h.paragraph_format.space_before = Pt(10 if name == "Heading 1" else 6)
    h.paragraph_format.space_after = Pt(3)
for name in ("List Bullet",):
    doc.styles[name].paragraph_format.space_after = Pt(2)


def para(text="", bold_lead=None, size=None, italic=False, align=None, after=None):
    p = doc.add_paragraph()
    if bold_lead:
        p.add_run(bold_lead).bold = True
    add_rich(p, text)
    for r in p.runs:
        if size:
            r.font.size = Pt(size)
        if italic:
            r.italic = True
    if align:
        p.alignment = align
    if after is not None:
        p.paragraph_format.space_after = Pt(after)
    return p


def add_rich(p, text):
    """`code` spans in Consolas, the rest plain."""
    parts = text.split("`")
    for i, s in enumerate(parts):
        if not s:
            continue
        r = p.add_run(s)
        if i % 2:
            r.font.name = "Consolas"
            r.font.size = Pt(9.5)


def bullet(text, lead=None):
    p = doc.add_paragraph(style="List Bullet")
    if lead:
        p.add_run(lead).bold = True
    add_rich(p, text)
    return p


def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:color"), "auto")
    sh.set(qn("w:fill"), hexcolor)
    tcPr.append(sh)


def set_widths(t, widths):
    """Column widths that Word and LibreOffice both honour: fixed layout, tblGrid and every cell."""
    t.autofit = False
    tblPr = t._tbl.tblPr
    lay = OxmlElement("w:tblLayout")
    lay.set(qn("w:type"), "fixed")
    tblPr.append(lay)
    tw = OxmlElement("w:tblW")
    tw.set(qn("w:w"), str(int(sum(widths) * 1440)))
    tw.set(qn("w:type"), "dxa")
    for old in tblPr.findall(qn("w:tblW")):
        tblPr.remove(old)
    tblPr.append(tw)
    for gc, w in zip(t._tbl.tblGrid.findall(qn("w:gridCol")), widths):
        gc.set(qn("w:w"), str(int(w * 1440)))
    for row in t.rows:
        row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
        for c, w in zip(row.cells, widths):
            c.width = Inches(w)


def table(rows, widths, header=True, size=9.5):
    t = doc.add_table(rows=len(rows), cols=len(rows[0]))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        for j, txt in enumerate(row):
            c = t.cell(i, j)
            c.width = Inches(widths[j])
            c.text = ""
            p = c.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            add_rich(p, txt)
            for r in p.runs:
                r.font.size = Pt(size)
                if header and i == 0:
                    r.bold = True
            if header and i == 0:
                shade(c, "DCE6F2")
        if header and i == 0:                         # repeat the header row after a page break
            t.rows[0]._tr.get_or_add_trPr().append(OxmlElement("w:tblHeader"))
    set_widths(t, widths)
    doc.add_paragraph().paragraph_format.space_after = Pt(0)
    return t


def figure(path, width, caption):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(1)
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(os.path.join(TOP, path), width=Inches(width))
    caption_par(caption)


def split_caption(cap):
    """'Fig. 3. Text' -> ('Fig. 3', ' Text'); the label is set in bold."""
    k = cap.index(".", cap.index(".") + 1)
    return cap[:k], cap[k + 1:]


def caption_par(caption):
    c = doc.add_paragraph()
    c.alignment = WD_ALIGN_PARAGRAPH.CENTER
    lead, rest = split_caption(caption)
    r = c.add_run(lead + ".")
    r.bold, r.font.size = True, Pt(9)
    r = c.add_run(rest)
    r.font.size = Pt(9)
    c.paragraph_format.space_after = Pt(8)


def figures_side(items, width):
    """Two figures next to each other in a borderless table, captions below each."""
    t = doc.add_table(rows=1, cols=len(items))       # one row: figure and caption stay together
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, (path, cap) in enumerate(items):
        p = t.cell(0, j).paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(os.path.join(TOP, path), width=Inches(width))
        c = t.cell(0, j).add_paragraph()
        c.alignment = WD_ALIGN_PARAGRAPH.CENTER
        lead, rest = split_caption(cap)
        r = c.add_run(lead + ".")
        r.bold, r.font.size = True, Pt(9)
        r = c.add_run(rest)
        r.font.size = Pt(9)
    set_widths(t, [6.7 / len(items)] * len(items))
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def code(lines):
    for k, l in enumerate(lines):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(0 if k < len(lines) - 1 else 6)
        p.paragraph_format.left_indent = Inches(0.25)
        r = p.add_run(l)
        r.font.name, r.font.size = "Consolas", Pt(9)


# ===================================================================== title
t = doc.add_paragraph()
r = t.add_run("2D potential energy surface for H₂ on Si(111)")
r.bold, r.font.size, r.font.color.rgb = True, Pt(18), ACCENT
t.paragraph_format.space_after = Pt(0)
para("PBE-D3 single points, upright H₂ above a surface Si atom — data, method, fit and files",
     size=11, italic=True, after=2)
para("V. Khokhar · October 2026 · project directory `Si111_H2/`", size=9.5, after=8)

para("This note describes a two-dimensional potential energy surface (PES) for an H₂ molecule "
     "standing upright above one surface Si atom of Si(111), computed with VASP (PBE + D3). "
     "The main data file is `fit/pes_r1r2_E0.csv` (1159 points). A smooth fitted surface with "
     "analytic gradients (`fit/pes_model.joblib`) is also available for dynamics and reaction-path "
     "work. Sections 1–3 explain how the data were made, section 4 the fit, section 5 the main "
     "results, and section 6 where everything is and how to use it.")

# ===================================================================== 1 system
doc.add_heading("1. System and coordinates", level=1)
lay = doc.add_table(rows=1, cols=2)
lay.alignment = WD_TABLE_ALIGNMENT.CENTER
left, right = lay.cell(0, 0), lay.cell(0, 1)
set_widths(lay, [4.6, 2.1])


def cell_par(cell, text, lead=None, first=False):
    p = cell.paragraphs[0] if first else cell.add_paragraph()
    p.paragraph_format.space_after = Pt(4)
    if lead:
        p.add_run(lead).bold = True
    add_rich(p, text)
    return p


cell_par(left, "Si(111) slab, 24 Si atoms, periodic cell 7.680 × 6.651 × 40.406 Å "
         "(about 32 Å of vacuum). The top surface is reconstructed (Si 1, 2, 13, 14 are "
         "three-coordinated); the bottom layer is not H-passivated. All Si atoms are frozen.",
         "Surface. ", first=True)
cell_par(left, "H₂ is collinear with the surface normal, directly above Si 13 (Si*), "
         "with both H atoms on that axis (a Si*–H1–H2 line).", "Molecule. ")
cell_par(left, "r1 = H1–H2 distance; r2 = height of the lower atom H1 above Si*. "
         "H2 then sits at r1 + r2 above Si*. In the scripts the same coordinates are called "
         "R (= r1) and h (= r2). All distances are in Å and energies in eV (only the figures in "
         "`figs_r1r2/` use bohr).", "Coordinates (Fig. 1). ")
cell_par(left, "r1 large, r2 large: two free H atoms. r1 = 0.75 Å, r2 large: gas-phase H₂. "
         "r1 large, r2 ≈ 1.5 Å: one H bonded to Si* (Si–H) and one free H. "
         "The PES connects H₂ approaching the surface with dissociative adsorption and, read "
         "backwards, with Eley–Rideal abstraction H + H–Si → H₂.", "Limits. ")
pr = right.paragraphs[0]
pr.alignment = WD_ALIGN_PARAGRAPH.CENTER
pr.add_run().add_picture(os.path.join(TOP, "figs_r1r2/fig1_coords.png"), width=Inches(1.9))
cp = right.add_paragraph()
cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = cp.add_run("Fig. 1.")
r.bold, r.font.size = True, Pt(9)
r = cp.add_run(" Coordinates r1 and r2 (Si* = Si 13).")
r.font.size = Pt(9)

# ===================================================================== 2 method
doc.add_heading("2. How the data were computed", level=1)
para("Every point is a separate single-point (no relaxation) self-consistent DFT calculation. "
     "Key VASP settings:")
table([
    ["Setting", "Value", "Why / note"],
    ["Functional", "PBE + D3(BJ) (`IVDW = 12`)", "D3 is an additive, non-self-consistent term; "
     "PBE alone = E0 − Edisp (Edisp in OUTCAR)"],
    ["Basis, precision", "PAW, `ENCUT = 400` eV, `PREC = Normal`", ""],
    ["k-points / smearing", "16 irreducible k-points; Gaussian, `SIGMA = 0.2` eV", "energy used: "
     "E0 (sigma → 0)"],
    ["Spin", "`ISPIN = 2`, `MAGMOM = 24*0 1 -1`", "slab non-magnetic, the two H opposite "
     "(broken-symmetry singlet): correct H₂ dissociation to two atoms"],
    ["Geometry", "`IBRION = -1`, `NSW = 0`", "frozen slab, H positions fixed by (r1, r2)"],
], widths=[1.35, 2.35, 3.0])

doc.add_heading("Grid and scan procedure", level=2)
para("29 heights r2 = 1.00–1.30 (step 0.1), 1.30–1.70 (0.05), 1.70–2.20 (0.1), "
     "2.50–8.00 Å (0.5), each with 40 H–H distances r1 from 8.00 down to 0.25 Å "
     "(dense around 0.6–2.4 Å). Each height is one scan (directory `<r2>/`, e.g. "
     "`1.50/`): it starts at r1 = 8 Å and walks inward, each point restarting from the "
     "previous point's wavefunction and charge density. Going inward matters: starting from "
     "closed-shell H₂ and walking outward keeps the restricted solution past the point where "
     "the H₂ spin symmetry breaks and gives energies that are far too high.")

figure("report/fig_pes_2d.png", 3.4, "Fig. 2. The DFT data (PBE+D3) relative to gas-phase H₂. White lines: "
       "the 29 computed heights. The light valley on the left is intact H₂; the one along the "
       "bottom (r2 ≈ 1.5 Å) is H bonded to Si* with the other H free. Colours capped at "
       "8 eV; the r1 = 0.25 Å wall reaches about 30 eV.")
doc.add_heading("Slab spin state and reruns", level=2)
para("The unpassivated slab has two nearly degenerate electronic states: a non-magnetic one and "
     "one with spin-polarized dangling bonds (about 0.2 eV lower). Almost all scan points sit in "
     "the non-magnetic state, but some SCF runs drifted partly toward polarization, which gave "
     "10–30 meV offsets for whole rows and about 20 single points 6–42 meV too high. "
     "Decision: use the non-magnetic slab everywhere, as one consistent reference. Suspicious "
     "points were recomputed in three rerun rounds with the total spin fixed (`NUPDOWN`), and "
     "a rerun was accepted only if the dangling-bond Si atoms stayed non-magnetic "
     "(Σ|m| < 0.15 μ_B). 175 rerun points replaced scan values; they now live in the "
     "height directories (old scan outputs kept in `R_<r1>/scan_old/`). One point "
     "(r2 = 3.00, r1 = 2.60) has no trustworthy non-magnetic energy and is left out, so the data "
     "set has 1159 of 1160 points. A leave-one-out check finds no outliers (rms per row "
     "≤ 3.1 meV).")

# ===================================================================== 3 data
doc.add_heading("3. The data: fit/pes_r1r2_E0.csv", level=1)
table([
    ["Column", "Meaning", "Unit"],
    ["r1", "H1–H2 distance (= R)", "Å"],
    ["r2", "height of H1 above Si* (= h)", "Å"],
    ["E0", "total energy of the slab + 2 H cell, PBE+D3, absolute (VASP E0)", "eV"],
], widths=[0.8, 4.6, 0.7])
para("Energies are absolute; subtract a reference before use:")
bullet("E(r1 = 8, r2 = 8) = −129.37379 eV: two H atoms far apart and far from the surface "
       "(the zero of the fitted model).", "Free atoms: ")
bullet("E(r1 = 0.75, r2 = 8) = −133.90415 eV, i.e. 4.530 eV below the atoms "
       "(≈ D_e of H₂). Fig. 2 and section 5 use this zero.", "Gas-phase H₂: ")
para("The full table `fit/pes_data.csv` has the same points with more columns: E0, Edisp "
     "(D3 part), E_pbe = E0 − Edisp, the total and H magnetic moments, SCF steps and the "
     "source (scan or rerun round). Use it to fit PBE without dispersion. Caveat: rows "
     "r2 = 3.0–4.0 Å keep 10–20 meV slab-state offsets, so features smaller than "
     "≈ 20 meV there are within the data uncertainty.")

# ===================================================================== 4 fit
doc.add_heading("4. Fitting the surface", level=1)
para("The energy range (−4.5 to +30 eV) and the flat asymptotes are hard for a single "
     "machine-learning model, so the fit learns only what is left after a physical baseline "
     "(Δ-learning). With energies relative to the free atoms:")
p = para("V(r1, r2) = V_HH(r1) + V_H(r2) + V_H(r1 + r2) + T(r1) T(r2) · I_GP(r1, r2) "
         "[+ D_GP(r1, r2)]", align=WD_ALIGN_PARAGRAPH.CENTER)
for rr in p.runs:
    rr.italic = True
bullet("the gas-phase H₂ curve, i.e. the r2 = 8 row of the data;", "V_HH(r1): ")
bullet("one H atom against the surface, the r1 = 8 column (Si–H curve). Both 1D curves are "
       "cubic splines with an exponential repulsive wall below the grid and constant above it;",
       "V_H(z): ")
bullet("Gaussian-process regression (scikit-learn, Matérn 5/2 + white-noise kernel, "
       "inputs exp(−r1), exp(−r2)) on the remaining three-body part. It is exactly zero on "
       "both 1D cuts and is switched off smoothly by T between 5.5 and 8 Å, so both asymptotic "
       "channels equal the 1D curves exactly and the surface stays physical outside the grid;",
       "I_GP: ")
bullet("a second, small GP for the D3 term, so dispersion can be switched off (PBE only).",
       "D_GP: ")
para("Gradients are analytic (spline derivatives plus the closed-form kernel derivative; checked "
     "against finite differences to 0.01 meV/Å). The GP also returns an uncertainty σ. "
     "Quality, for energies below 3 eV above gas-phase H₂:")
table([
    ["Test", "RMSE", "Max"],
    ["Training error (all 1159 points)", "0.9 meV", "6.2 meV"],
    ["Leave a whole height out (r2 = 1.35, 1.45, 1.55, 1.65, 2.10)", "1.0–3.4 meV", "≤ 13.4 meV"],
    ["Leave out r2 = 3.50 (row with a slab-state offset)", "14 meV", "23 meV"],
    ["Leave out every other r1", "3.4 meV", "14 meV"],
    ["Gas-phase H₂ from the fit", "r_e 0.7509 Å, D_e 4.531 eV, ω_e 4276 cm⁻¹", ""],
], widths=[3.9, 1.7, 0.9])
figures_side([
    ("fit/fit_pes_3d.png", "Fig. 3. Fitted surface (PBE+D3) with the DFT points."),
    ("fit/residuals_map.png", "Fig. 4. Fit − DFT at every point (meV)."),
], width=3.25)

# ===================================================================== 5 results
doc.add_heading("5. Main results (PBE+D3, relative to gas-phase H₂)", level=1)
table([
    ["Quantity", "Value"],
    ["H + H–Si (r1 → ∞, r2 = 1.50 Å)", "+1.277 eV (PBE: +1.245 eV)"],
    ["Si–H bond (frozen slab)", "1.501 Å, binding 3.254 eV vs free H, 2022 cm⁻¹"],
    ["Si–H binding with the slab relaxed (bottom bilayer fixed)", "3.159 eV (−0.095 eV)"],
    ["Physisorbed H₂", "−52 meV at r2 = 2.85 Å (mostly D3)"],
    ["Barrier on the minimum-energy path", "none beyond 10 meV"],
], widths=[3.6, 2.9])
para("The minimum-energy path (string method on the fitted surface, Fig. 5) has no barrier: "
     "dissociating H₂ onto Si* costs only the endothermicity of +1.28 eV, and Eley–Rideal "
     "abstraction H + H–Si → H₂ is barrierless and exothermic by the same amount in this "
     "collinear geometry. The relaxed-slab number comes from two separate relaxations "
     "(`relax_*_fixbot/`); the PES itself is frozen-slab.")
figures_side([
    ("fit/mep_2d.png", "Fig. 5. Minimum-energy path on the fitted PES."),
    ("fit/mep_profile.png", "Fig. 6. Energy and coordinates along the path."),
], width=3.05)

# ===================================================================== 6 files
doc.add_heading("6. Where things are and how to use them", level=1)
table([
    ["Location", "Content / purpose"],
    ["`fit/pes_r1r2_E0.csv`", "main data: r1, r2, E0 (1159 points)"],
    ["`fit/pes_data.csv`", "same points, all columns (Edisp, PBE energy, moments, source)"],
    ["`fit/pes_model.joblib`, `pes_model.py`", "fitted PES (load with joblib; class `PES`)"],
    ["`<r2>/R_<r1>/`", "VASP inputs/outputs of every point (INCAR, POSCAR, OUTCAR, OSZICAR); "
     "`source.txt` if from a rerun"],
    ["`<r2>/scan.sh`, `submit`, `submit_all`", "run a height scan on the cluster (SLURM); "
     "`scan.sh --check` = complete?"],
    ["`<r2>/<r2>.dat`", "per-height table: r1, F, E0, moments, SCF steps, status"],
    ["`fit/collect_data.py`", "gathers converged points → both CSV files"],
    ["`fit/qc.py`", "leave-one-out outlier and spin-state checks → `qc_report.txt`"],
    ["`fit/fit_pes.py`", "fit + cross-validation → `pes_model.joblib`, `fit_report.txt`"],
    ["`fit/plot_fit.py`, `analyze.py`", "fit figures; MEP and energetics → `analyze_report.txt`"],
    ["`plot_pes2d3d.py`, `figs_r1r2/`", "data maps; paper-style figures in bohr (`make_figs.sh`)"],
    ["`relax_SiH_fixbot/`, `relax_Si111_fixbot/`", "slab relaxations, Si–H binding "
     "(`summary.py`, `binding.py`)"],
    ["`CLAUDE.md`", "detailed project notes (all settings, decisions and pitfalls)"],
], widths=[2.55, 4.15], size=9)
para("Rebuild everything after new data: in `fit/`, run `collect_data.py`, `qc.py`, "
     "`fit_pes.py` (≈ 20 min), `plot_fit.py`, `analyze.py`. Using the fitted PES in Python "
     "(run from `fit/`, or put it on the path):")
code([
    "import joblib, pes_model                     # pes_model.py must be importable",
    "pes = joblib.load('pes_model.joblib')",
    "E  = pes.energy(r1, r2, d3=True, ref='gas')  # eV; ref='atoms' -> zero at free atoms",
    "g  = pes.gradient(r1, r2)                    # (dV/dr1, dV/dr2) in eV/A",
    "s  = pes.sigma(r1, r2)                       # GP uncertainty (eV)",
    "E2 = pes.energy_com(r1, Z)                   # from the H2 centre-of-mass height Z",
])
para("Limits of validity: frozen slab, upright H₂ on a single site (2D), PBE+D3 with a "
     "non-magnetic slab; the data cover r1 = 0.25–8 Å and r2 = 1.0–8 Å. Outside "
     "the grid the model falls back to the 1D curves (well behaved, but not checked against "
     "DFT); the gap between r2 = 2.2 and 2.5 Å and the sparser rows above are the least "
     "certain regions (see `fit/sigma_map.png`).", size=9.5)

out = os.path.join(HERE, "PES_report.docx")
doc.save(out)
print("wrote", out)
