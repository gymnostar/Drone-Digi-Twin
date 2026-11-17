#!/bin/sh
# Overnight: run coarse -> medium -> fine one after the other.
# Start from the repo root with:
#   nohup caffeinate -i /Applications/OpenFOAM-v2306.app/Contents/MacOS/openfoam -c "sh DLR/meshStudy/run_all.sh" > DLR/meshStudy/run_all.out 2>&1 &
# (caffeinate -i stops the Mac from sleeping; keep it plugged in.)
cd "$(dirname "$0")" || exit 1
for c in coarse medium fine; do
    echo "=== $c start $(date)"
    ./$c/Allrun || echo "!!! $c FAILED, see $c/log.*"
    echo "=== $c end   $(date)"
done
echo "ALL DONE $(date). Next: .venv/bin/python DLR/meshStudy/make_table.py"
