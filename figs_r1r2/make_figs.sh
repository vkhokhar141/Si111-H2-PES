#!/bin/bash
# Figures in r1/r2 coordinates (r1 = H1-H2 = R, r2 = H1 above Si* = h), in bohr.
#   fig1_coords.png  coordinate system (VMD, TachyonInternal)
#   fig2a_pes3d.png, fig2b_cuts.png, fig2c_vh.png, fig2.png  (gnuplot)
# Needs fit/pes_model.joblib and fit/pes_data.csv. Units: UNIT in make_data.py.
set -e
cd "$(dirname "$0")"
python3 make_data.py
gnuplot fig2.gp
python3 fig1.py
