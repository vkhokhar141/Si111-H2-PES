# Fig. 2 analogue for H2 on Si(111): r1 = H1-H2, r2 = H1 above Si* (Si 13).
#   (a) fitted PES V(r1, r2) as a mesh with contours on the base
#   (b) V(r1) with H1 chemisorbed (r2 at the Si-H minimum) and H1 away from the surface
#   (c) H1-Si* potential with H2 far away: frozen slab ("sudden"); a relaxed-slab curve is
#       added if vh_relaxed.dat exists
# Data from make_data.py.   Run in figs_r1r2/:  gnuplot fig2.gp
#   -> fig2a_pes3d.png, fig2b_cuts.png, fig2c_vh.png, fig2.png (all three)

load "labels.gp"
PURPLE = "#9400d3"
RED    = "#e8000b"
GREEN  = "#009e73"
have_relaxed = system("test -s vh_relaxed.dat && echo 1 || echo 0") + 0
FONT = "DejaVu Sans,13"

# ---------------------------------------------------------------- panel macros
# every panel starts with reset; its placement and (a)/(b)/(c) tag come from these
# variables, which reset does not clear (full canvas and no tag for the single panels)
place = 'set origin ox, oy; set size sx, sy; if (strlen(tag) > 0) { set label 99 tag at screen tx, screen ty center font ",15" }'
ox = 0; oy = 0; sx = 1; sy = 1; tag = ""; tx = 0; ty = 0
panel_a = ' \
  reset; load "labels.gp"; @place; \
  set hidden3d front; set view 42, 36, 1, 1.15; \
  set xrange [r1min:12]; set yrange [1.9:12]; set zrange [0:vtop]; \
  set xyplane at -0.6; set ticslevel 0.15; \
  set xtics 2, 2, 12 offset 0, -0.4; set ytics 2, 2, 12 offset 0.6, -0.4; set ztics 0, (vtop <= 12 ? 1.5 : 3), vtop; \
  set xlabel "r_1/".unit offset 0, -1.2 rotate parallel; set ylabel "r_2/".unit offset 1, -0.6 rotate parallel; \
  set zlabel "V (eV)" rotate by 90 offset 1.5, 0; \
  set border 31 lw 1; unset key; \
  set contour base; set cntrparam levels discrete 0.1, 0.3, 0.6, 1.0, 1.4, 1.8, 2.4, 3.0, 3.6, 4.2; \
  splot "pes_grid.dat" u 1:2:3 w l lc rgb PURPLE lw 0.6 nocontours, \
        "pes_grid.dat" u 1:2:3 w l lw 1.2 nosurface'

panel_b = ' \
  reset; load "labels.gp"; @place; \
  set xrange [0.5:10]; set yrange [-4.8:1.5]; set xtics 2; set ytics 1.5; \
  set xlabel "r_1/".unit; set ylabel "V (eV)"; set grid lc rgb "#e4e3df"; \
  set key at graph 0.97, graph 0.42 right box opaque samplen 2.5 spacing 1.2; \
  plot "cut_chem.dat" u 1:2 w l lc rgb PURPLE lw 2 t "H1 chemisorbed", \
       "pts_chem.dat" u 1:2 w p pt 6 ps 0.8 lc rgb PURPLE notitle, \
       "cut_away.dat" u 1:2 w l lc rgb RED lw 2 t "H1 away", \
       "pts_away.dat" u 1:2 w p pt 1 ps 0.9 lc rgb RED notitle'

panel_c = ' \
  reset; load "labels.gp"; @place; \
  set xrange [1.5:15.5]; set yrange [-3.6:1.0]; set xtics 1.5, 2, 15.5; set ytics 0.5; \
  set xlabel "Distance of H1 above the surface Si* (".unit.")"; set ylabel "V (eV)"; \
  set grid lc rgb "#e4e3df"; set key right center box opaque samplen 2.5 spacing 1.2; \
  plot "vh_sudden.dat" u 1:2 w l lc rgb GREEN lw 2 t "sudden (frozen slab)", \
       "vh_pts.dat" u 1:2 w p pt 6 ps 0.8 lc rgb GREEN notitle' . \
  (have_relaxed ? ', "vh_relaxed.dat" u 1:2 w l lc rgb PURPLE lw 2 t "relaxed"' : '')

# ---------------------------------------------------------------- single panels
set terminal pngcairo size 1000, 800 enhanced font FONT
set output "fig2a_pes3d.png"
@panel_a
set terminal pngcairo size 760, 560 enhanced font FONT
set output "fig2b_cuts.png"
@panel_b
set output "fig2c_vh.png"
@panel_c

# ---------------------------------------------------------------- all three, laid out as Fig. 2
set terminal pngcairo size 1500, 1500 enhanced font FONT
set output "fig2.png"
set multiplot
ox = 0.0; oy = 0.40; sx = 1.0; sy = 0.60; tag = "(a)"; tx = 0.40; ty = 0.415
@panel_a
ox = 0.0; oy = 0.035; sx = 0.5; sy = 0.37; tag = "(b)"; tx = 0.27; ty = 0.012
@panel_b
ox = 0.5; oy = 0.035; sx = 0.5; sy = 0.37; tag = "(c)"; tx = 0.77; ty = 0.012
@panel_c
unset multiplot
unset output
