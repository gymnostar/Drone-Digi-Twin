#!/bin/zsh
set -e
cd -- "${0:A:h}"
scene="${1:-joint-baseline}"
case "$scene" in
  joint-baseline|turbine-cheap|turbine-reference) ;;
  *) print 'Choose joint-baseline, turbine-cheap or turbine-reference'; exit 1 ;;
esac
apps=(/Applications/ParaView*.app(N))
if (( ${#apps} == 0 )); then
  print 'Install ParaView, then run this launcher again.'; exit 1
fi
export DRONE_VIEW="$scene"
export DRONE_VIEW_DIR="$PWD"
exec "${apps[-1]}/Contents/MacOS/paraview" --script="$PWD/view.py"
