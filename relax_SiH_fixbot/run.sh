#!/bin/bash
# Relaxation of the Si(111) slab + one H on Si 13 (bottom bilayer Si 3 4 5 6 9 10 23 24 fixed, rest free, cell fixed).
# Restartable: if a previous run stopped before EDIFFG was reached (time limit), its
# outputs are kept in run_<n>/ and the relaxation continues from CONTCAR (+ WAVECAR).
#
#   bash run.sh            run / continue (needs VASP_CMD from the job script)
#   bash run.sh --check    exit 0 if the relaxation has converged
#   POSCAR.orig is the starting structure (selective dynamics; scan slab geometry, H 1.50 A above Si 13).

set -u
cd "$(dirname "$0")"

relaxed() { [ -f OUTCAR ] && grep -q "reached required accuracy" OUTCAR; }

if [ "${1:-}" = "--check" ]; then relaxed; exit; fi
: "${VASP_CMD:?set VASP_CMD to the full mpi + vasp_std command}"

if relaxed; then
    echo "already relaxed"
    exit 0
fi

if [ -f OUTCAR ]; then
    # keep the previous run, continue from its last geometry
    n=1; while [ -e "run_$n" ]; do n=$((n + 1)); done
    mkdir "run_$n"
    cp POSCAR INCAR OUTCAR OSZICAR vasp.out "run_$n/" 2>/dev/null
    cp XDATCAR CONTCAR "run_$n/" 2>/dev/null
    natoms=$(awk 'NR == 7 {s = 0; for (i = 1; i <= NF; i++) s += $i; print s}' POSCAR.orig)
    head=$(awk 'NR == 8 {print ($1 ~ /^[sS]/) ? 9 : 8}' POSCAR.orig)   # + "Selective dynamics"
    if [ "$(grep -c . CONTCAR 2>/dev/null || echo 0)" -ge $((head + natoms)) ]; then
        cp CONTCAR POSCAR
        echo "continuing from CONTCAR of run_$n"
    fi
    if [ -s WAVECAR ]; then
        sed -i -e 's/^ISTART *=.*/ISTART =  1/' -e 's/^ICHARG *=.*/ICHARG =  0/' INCAR
    elif [ -s CHGCAR ]; then
        sed -i -e 's/^ISTART *=.*/ISTART =  0/' -e 's/^ICHARG *=.*/ICHARG =  1/' INCAR
    fi
fi

echo "start $(date)"
$VASP_CMD > vasp.out 2>&1 < /dev/null
if relaxed; then
    echo "relaxed $(date)"
    rm -f CHG
    python3 summary.py > summary.txt 2>&1 && cat summary.txt
else
    echo "not converged yet $(date): resubmit (bash submit) to continue"
fi
