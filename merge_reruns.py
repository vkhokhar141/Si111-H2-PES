"""Move the rerun points that the fit uses into their height directories.

For every point of fit/pes_data.csv whose src is a rerun round (reruns, reruns2, ...), the
scan outputs in <h>/R_<R>/ go to <h>/R_<R>/scan_old/, the rerun's files replace them, and
<h>/R_<R>/source.txt records where they came from. Then <h>/<h>.dat is rebuilt with
VASP_CMD=false bash <h>/scan.sh. Afterwards the rerun directories are no longer needed:
collect_data.py takes src from source.txt.

Run in the project directory after fit/collect_data.py:
  python3 merge_reruns.py          dry run (lists the moves)
  python3 merge_reruns.py --do     move
"""
import os
import shutil
import subprocess
import sys
import pandas as pd

TOP = os.path.dirname(os.path.abspath(__file__))
DO = "--do" in sys.argv

df = pd.read_csv(os.path.join(TOP, "fit", "pes_data.csv"))
heights = set()
for _, r in df[df.src != "scan"].iterrows():
    hs, Rs = f"{r.h:.2f}", f"{r.R:.2f}"
    src = os.path.join(TOP, r.src, f"{hs}_R_{Rs}")
    dst = os.path.join(TOP, hs, f"R_{Rs}")
    if os.path.isfile(os.path.join(dst, "source.txt")):
        print(f"{hs}/R_{Rs}: already merged")
        continue
    if not os.path.isfile(os.path.join(src, "OUTCAR")):
        raise SystemExit(f"{src}: no OUTCAR")
    print(f"{r.src}/{hs}_R_{Rs} -> {hs}/R_{Rs}   (scan outputs -> {hs}/R_{Rs}/scan_old/)")
    if not DO:
        continue
    old = os.path.join(dst, "scan_old")
    os.makedirs(old)
    for f in os.listdir(dst):
        if f != "scan_old":
            shutil.move(os.path.join(dst, f), old)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy2(p, dst)
    nup = [l.split("=")[1].split()[0] for l in open(os.path.join(src, "INCAR"))
           if l.split("=")[0].strip() == "NUPDOWN"]
    with open(os.path.join(dst, "source.txt"), "w") as f:
        f.write(f"{r.src}\n# rerun {r.src}/{hs}_R_{Rs}"
                f" (NUPDOWN {nup[0] if nup else 'unset'}), merged by merge_reruns.py;"
                f" scan outputs in scan_old/\n")
    heights.add(hs)

for hs in sorted(heights):
    subprocess.run(["bash", os.path.join(TOP, hs, "scan.sh")], env={**os.environ, "VASP_CMD": "false"},
                   stdout=subprocess.DEVNULL, check=True)
    print(f"rebuilt {hs}/{hs}.dat")
if not DO:
    print("dry run; use --do to move")
