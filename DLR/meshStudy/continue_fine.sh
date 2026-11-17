#!/bin/sh
# Fine mesh only: continue from iteration 1000 to 2000, because its residuals were still falling at 1000.
# Run inside the OpenFOAM shell. New probe/residual lines are appended to the original files so
# make_table.py and the report read one continuous history.
cd "$(dirname "$0")/fine" || exit 1
sed -i '' -E 's/^endTime +[0-9]+;/endTime         2000;/' system/controlDict
rm -f log.decomposePar log.reconstructPar
decomposePar -latestTime > log.decomposePar 2>&1 || exit 1
mpirun -np 8 simpleFoam -parallel >> log.simpleFoam 2>&1 || exit 1
reconstructPar -latestTime > log.reconstructPar 2>&1 && rm -rf processor*
for d in probes residuals; do
  for f in postProcessing/$d/1000/*; do
    grep -v '^#' "$f" >> "postProcessing/$d/0/$(basename "$f" | sed -E 's/_1000//')"
  done
done
echo "FINE CONTINUED $(date)"
