#!/bin/sh
# Wait for the first atlas pass to finish, then run it again: only the runs without DONE (255) execute.
cd "$(dirname "$0")/../.." || exit 1
until grep -q "ALL DONE" DLR/atlas/run_atlas.out 2>/dev/null; do sleep 60; done
echo "== chain: second pass $(date)" >> DLR/atlas/run_atlas.out
sh DLR/atlas/run_atlas.sh >> DLR/atlas/run_atlas.out 2>&1
