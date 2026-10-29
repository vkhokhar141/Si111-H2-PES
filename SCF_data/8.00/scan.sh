#!/bin/bash
# H-H stretch scan for H2 standing upright above Si atom 13 (rest atom) of the
# frozen Si(111) slab. The height of the lower H above that Si is taken from the
# name of this directory (e.g. 6.00/ -> 6.00 A).
# Runs from large R inward; every point restarts from the WAVECAR/CHGCAR of
# the previous converged point, so the slab keeps one magnetic state and the
# H2 follows the broken-symmetry solution down to the closed-shell region.
# Re-running skips converged points and continues from the last one.
#
#   bash scan.sh           run the scan (needs VASP_CMD from the job script,
#                          e.g. "mpiexec.hydra -np 32 .../vasp_std")
#   bash scan.sh --check   exit 0 if every point is already converged

set -u

SCANDIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCANDIR"

H=$(basename "$SCANDIR")
[[ $H =~ ^[0-9]+\.[0-9][0-9]$ ]] || { echo "directory name '$H' is not a height like 6.00" >&2; exit 2; }
ZSI=7.944237044           # z of Si atom 13, directly below the H2
ZH1=$(bc -l <<< "$ZSI+$H")   # z of the fixed (lower) H; the upper H sits at ZH1+R
RLIST="8.00 7.00 6.00 5.00 4.40 3.90 3.50 3.10 2.80 2.60 2.40 2.30 2.20 2.10
       2.00 1.95 1.90 1.85 1.80 1.75 1.70 1.65 1.50 1.40 1.30 1.20 1.10 1.00
       0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.45 0.35 0.25"
DAT="$SCANDIR/$H.dat"

converged() { [ -f "$1/OUTCAR" ] && grep -q "aborting loop because EDIFF is reached" "$1/OUTCAR"; }

if [ "${1:-}" = "--check" ]; then
    for R in $RLIST; do converged "R_$R" || exit 1; done
    exit 0
fi

: "${VASP_CMD:?set VASP_CMD to the full mpi + vasp_std command}"

write_poscar() {   # $1 = z of the upper H
cat > POSCAR <<!
silicon24 hydrogen2
1.0
        7.6797451973         0.0000000000         0.0000000000
        0.0000000000         6.6508541107         0.0000000000
        0.0000000000         0.0000000000        40.4057273865
   Si    H
   24    2
Cartesian
     3.459346970         1.687071863         8.065362953
     7.278459988         1.687071863         8.065362953
     5.368903136         0.531619865         1.590662867
     1.529030537         0.531619865         1.590662867
     3.448967294         3.857046920         1.590662867
     7.288839893         3.857046920         1.590662867
     3.448967294         3.854631126         4.022291171
     7.288839893         3.854631126         4.022291171
     3.448967294         1.639490941         0.896186703
     7.288839893         1.639490941         0.896186703
     5.368903136         0.529204219         4.022291171
     1.529030537         0.529204219         4.022291171
     5.368903136         4.954664348         7.944237044
     1.529030537         4.985048493         7.570820630
     5.368903136         2.722188439         7.201458706
     1.529030537         2.687076370         7.213625184
     3.487373537         6.167982631         7.212074194
     7.250433422         6.167982631         7.212074194
     5.368903136         2.746287912         4.795579329
     1.529030537         2.746287912         4.795579329
     3.448967294         6.071715165         4.795579329
     7.288839893         6.071715165         4.795579329
     5.368903136         4.964917798         0.896186703
     1.529030537         4.964917798         0.896186703
     5.368903136         4.974062852        $ZH1
     5.368903136         4.974062852        $1
!
}

write_table() {
    printf "%-6s %-14s %-16s %-16s %-9s %-8s %-8s %-6s %s\n" \
        "#R(A)" "z_H2(A)" "F(eV)" "E0(eV)" "mag" "m(H1)" "m(H2)" "steps" "status" > "$DAT"
    for R in $RLIST; do
        d="R_$R"
        [ -f "$d/OSZICAR" ] || continue
        line=$(grep " F= " "$d/OSZICAR" | tail -1)
        [ -n "$line" ] || continue
        F=$(awk '{print $3}' <<< "$line")
        E0=$(awk '{print $5}' <<< "$line")
        mag=$(awk '{print $NF}' <<< "$line")
        mH=$(awk '/magnetization \(x\)/{b=1} b && $1=="25"{a=$NF} b && $1=="26"{c=$NF; b=0} END{print (a==""?"-":a), (c==""?"-":c)}' "$d/OUTCAR")
        n=$(grep -cE "^(DAV|RMM|CG ):" "$d/OSZICAR")
        st=$(converged "$d" && echo ok || echo NOT_CONVERGED)
        printf "%-6s %-14s %-16s %-16s %-9s %-8s %-8s %-6s %s\n" \
            "$R" "$(bc -l <<< "$ZH1+$R")" "$F" "$E0" "$mag" $mH "$n" "$st" >> "$DAT"
    done
}

prev=""
for R in $RLIST; do
    d="R_$R"
    if converged "$d"; then
        echo "R = $R already converged, skipping"
        prev="$d"
        continue
    fi

    mkdir -p "$d"
    cp INCAR KPOINTS POTCAR "$d/"
    cd "$d"
    write_poscar "$(bc -l <<< "$ZH1+$R")"

    rm -f WAVECAR CHGCAR CHG
    if [ -n "$prev" ] && [ -s "../$prev/CHGCAR" ]; then
        # CHGCAR carries the spin density; WAVECAR is only kept for the newest point
        cp "../$prev/CHGCAR" .
        sed -i 's/^ICHARG *=.*/ICHARG =  1/' INCAR
        if [ -s "../$prev/WAVECAR" ]; then
            cp "../$prev/WAVECAR" .
            sed -i 's/^ISTART *=.*/ISTART =  1/' INCAR
        fi
        echo "R = $R  restart from $prev   $(date)"
    else
        echo "R = $R  fresh start (MAGMOM guess)   $(date)"
    fi

    $VASP_CMD > vasp.out 2>&1
    cd "$SCANDIR"

    if converged "$d"; then
        # only the newest converged point needs its WAVECAR (~0.5 GB); CHGCARs are kept
        [ -n "$prev" ] && rm -f "$prev/WAVECAR" "$prev/CHG"
        prev="$d"
    else
        echo "WARNING: R = $R did not converge in NELM steps; next point restarts from $prev"
    fi
    write_table
done

write_table
echo "scan finished $(date)"
