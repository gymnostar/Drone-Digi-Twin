#!/bin/sh
# Continue every atlas run whose probe drift (last 50 iterations, probes >= 10 m, / Uref) is above
# the criterion in extend_runs.py (default 0.5 %), by +200 iterations, resample its grid, and record
# how much the field changed. Then refresh validation, figures, planner, deck and the report once.
# Waits for run_atlas.out to say "ALL DONE" so it never competes with the queue.
cd "$(dirname "$0")/../.." || exit 1
OF=/Applications/OpenFOAM-v2306.app/Contents/Resources/etc/openfoam
PV=/Applications/ParaView-6.2.0.app/Contents/bin/pvpython
PY=.venv/bin/python
until grep -q "ALL DONE" DLR/atlas/run_atlas.out 2>/dev/null; do sleep 60; done
$PY DLR/atlas/validate_atlas.py > /dev/null
for run in $($PY DLR/atlas/extend_runs.py --list); do
    C=DLR/atlas/cases.noindex/$run; W=DLR/atlas/cases.noindex/ext_$run
    [ -d "$C" ] || continue
    echo "== extend $run start $(date '+%H:%M')"
    rm -rf "$W"; cp -r "$C" "$W"; rm -f "$W/DONE"
    $PY - "$W/system/controlDict" <<'PY'
import re, sys
p = sys.argv[1]; t = open(p).read()
old = int(re.findall(r"(?m)^endTime\s+(\d+);", t)[0])
t = re.sub(r"startFrom\s+\w+;", "startFrom       latestTime;", t)
t = re.sub(r"endTime\s+\d+;", f"endTime         {old + 200};", t)
t = re.sub(r"writeControl\s+\w+;", "writeControl    runTime;", t, count=1)
t = re.sub(r"writeInterval\s+\d+;", "writeInterval   50;", t, count=1)
t = re.sub(r"purgeWrite\s+\d+;", "purgeWrite      2;", t, count=1)
open(p, "w").write(t)
PY
    "$OF" -c "cd '$W' && . \$WM_PROJECT_DIR/bin/tools/RunFunctions && rm -rf processor* log.decomposePar log.simpleFoam log.reconstructPar && \
        runApplication decomposePar -latestTime && runParallel simpleFoam && runApplication reconstructPar -latestTime && rm -rf processor*" \
        || { echo "!!! extend $run FAILED (see $W/log.*)"; continue; }
    cp DLR/atlas/grid.noindex/$run.npz DLR/atlas/grid.noindex/prev_$run.npz        # keep the shorter run's field
    "$PV" DLR/atlas/sample_grid.py "$W" "$PWD/DLR/atlas/grid.noindex/$run.npz" > "$W/log.sample" 2>&1 || echo "!!! sampling $run failed"
    $PY DLR/atlas/extend_runs.py --compare "$run"
    echo "== extend $run end $(date '+%H:%M')"
done
$PY DLR/atlas/validate_atlas.py > /dev/null && $PY DLR/atlas/make_figures.py > /dev/null
$PY DLR/planner/run_experiments.py > DLR/atlas/extend_planner.log 2>&1
sh DLR/slides/2026-11-26/build.sh; $PY docs/overnight-2026-11-25/build_report.py
echo "EXTEND DONE $(date)"
