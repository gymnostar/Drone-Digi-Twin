#!/bin/sh
# Iteration-error check: continue a 350-iteration atlas run to 600 iterations and compare the fields.
#   sh DLR/atlas/iteration_check.sh [run=d060_u05] [end=600]
# Waits for the atlas queue (run_atlas.out: ALL DONE) first, so it never competes with it.
# Output: DLR/atlas/results/iteration_check.json (|dU|/Uref between iteration 350 and `end`).
cd "$(dirname "$0")/../.." || exit 1
RUN="${1:-d060_u05}"; END="${2:-600}"
C=DLR/atlas/cases.noindex/$RUN; W=DLR/atlas/cases.noindex/iter_$RUN
OF=/Applications/OpenFOAM-v2306.app/Contents/Resources/etc/openfoam
PV=/Applications/ParaView-6.2.0.app/Contents/bin/pvpython
until grep -q "ALL DONE" DLR/atlas/run_atlas.out 2>/dev/null; do sleep 60; done
echo "== iteration check $RUN -> $END start $(date '+%H:%M')"
rm -rf "$W"; cp -r "$C" "$W" && rm -f "$W/DONE"
LAST=$(ls "$W" | grep -E '^[0-9]+$' | grep -v '^0$' | sort -n | tail -1)
.venv/bin/python - "$W/system/controlDict" "$END" <<'PY'
import re, sys
p, end = sys.argv[1], sys.argv[2]; t = open(p).read()
t = re.sub(r"startFrom\s+\w+;", "startFrom       latestTime;", t)
t = re.sub(r"endTime\s+\d+;", f"endTime         {end};", t)
t = re.sub(r"writeControl\s+\w+;", "writeControl    runTime;", t, count=1)
t = re.sub(r"writeInterval\s+\d+;", "writeInterval   50;", t, count=1)
t = re.sub(r"purgeWrite\s+\d+;", "purgeWrite      2;", t, count=1)
open(p, "w").write(t)
PY
"$OF" -c "cd '$W' && . \$WM_PROJECT_DIR/bin/tools/RunFunctions && rm -rf processor* log.decomposePar log.simpleFoam log.reconstructPar && \
    runApplication decomposePar -latestTime && runParallel simpleFoam && runApplication reconstructPar -latestTime && rm -rf processor*" \
  || { echo "iteration check FAILED (see $W/log.*)"; exit 1; }
NEW=$(ls "$W" | grep -E '^[0-9]+$' | sort -n | tail -1)
echo "continued from $LAST to $NEW"
"$PV" DLR/atlas/sample_grid.py "$W" "$PWD/DLR/atlas/grid.noindex/iter_${RUN}.npz" > "$W/log.sample" 2>&1 || exit 1
.venv/bin/python DLR/atlas/iteration_check.py "$RUN" "$NEW" "$LAST"
