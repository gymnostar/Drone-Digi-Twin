#!/bin/sh
# Double-click in Finder: opens ParaView with the DLR wind view ready.
cd -- "$(dirname -- "$0")"
exec /Applications/ParaView-6.2.0.app/Contents/MacOS/paraview --script="$PWD/scripts/wind_view.py"
