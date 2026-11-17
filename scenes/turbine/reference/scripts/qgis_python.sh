#!/usr/bin/env bash
# macOS QGIS 4.2.3 bundled Python/GDAL; does not open a desktop renderer.
set -euo pipefail
QGIS_APP="${QGIS_APP:-/Applications/QGIS-final-4_2_3.app}"
[[ -x "$QGIS_APP/Contents/MacOS/python3.12" ]] || { echo 'Set QGIS_APP to the QGIS 4.2.3 application path.' >&2; exit 2; }
export PYTHONHOME="$QGIS_APP/Contents/Frameworks"
export PYTHONPATH="$QGIS_APP/Contents/Resources/qgis/python:$QGIS_APP/Contents/Resources/python3.12/site-packages"
export PROJ_DATA="$QGIS_APP/Contents/Resources/qgis/proj"
export GDAL_DATA="$QGIS_APP/Contents/Resources/qgis/gdal"
export QT_QPA_PLATFORM=offscreen
exec "$QGIS_APP/Contents/MacOS/python3.12" "$@"
