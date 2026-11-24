#!/bin/sh
# Overnight wind atlas: one mesh, 12 directions + one 10 m/s run, solved one after the other.
# Start from the repo root:
#   nohup caffeinate -i sh DLR/atlas/run_atlas.sh > DLR/atlas/run_atlas.out 2>&1 &
# Each run: split inlet/outlet for its direction -> simpleFoam (450 it, 8 cores) -> keep the last
# time -> resample onto the shared 5 m grid (grid.noindex/<run>.npz). A run with a DONE file is skipped.
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
A="$ROOT/DLR/atlas"; CASES="$A/cases.noindex"
OF=/Applications/OpenFOAM-v2306.app/Contents/Resources/etc/openfoam
PV=/Applications/ParaView-6.2.0.app/Contents/bin/pvpython
PY="$ROOT/.venv/bin/python3"
mkdir -p "$A/grid.noindex"

"$PY" "$A/setup_atlas.py" base || exit 1
B="$CASES/base"
if [ ! -f "$B/BASE_OK" ]; then
    "$OF" -c "cd '$B' && createPatch -overwrite -dict system/createPatchDict.merge > log.createPatch.merge 2>&1"
    grep -q "^End" "$B/log.createPatch.merge" && touch "$B/BASE_OK" || { echo "base merge FAILED"; exit 1; }
fi
"$PY" "$A/setup_atlas.py" || exit 1

for run in $("$PY" -c "import sys; sys.path.insert(0,'$A'); import setup_atlas as s; print(' '.join(s.name(t,u) for t,u in s.RUNS))"); do
    C="$CASES/$run"
    [ -f "$C/DONE" ] && { echo "== $run already done"; continue; }
    echo "== $run start $(date '+%H:%M:%S')  free disk $(df -h "$ROOT" | awk 'NR==2{print $4}')"
    "$OF" -c "cd '$C' && . \$WM_PROJECT_DIR/bin/tools/RunFunctions && \
        runApplication -s inletOutlet topoSet -dict system/topoSetDict.inletOutlet && \
        runApplication createPatch -overwrite && \
        rm -rf 0 && cp -r 0.orig 0 && \
        runApplication decomposePar && \
        runParallel simpleFoam && \
        runApplication reconstructPar -latestTime && rm -rf processor*" \
      || { echo "!!! $run solver FAILED (see $C/log.*)"; continue; }
    grep -q "^End" "$C/log.simpleFoam" || { echo "!!! $run simpleFoam did not finish"; continue; }
    "$PV" "$A/sample_grid.py" "$C" "$A/grid.noindex/$run.npz" > "$C/log.sample" 2>&1 \
      || echo "!!! $run sampling FAILED (see $C/log.sample)"
    touch "$C/DONE"
    echo "== $run end   $(date '+%H:%M:%S')"
done
echo "ALL DONE $(date)"
