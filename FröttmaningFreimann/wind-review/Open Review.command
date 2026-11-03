#!/bin/sh
set -eu
cd -- "$(dirname -- "$0")"
exec /Applications/ParaView-6.2.0.app/Contents/MacOS/paraview --state="$PWD/review.pvsm"
