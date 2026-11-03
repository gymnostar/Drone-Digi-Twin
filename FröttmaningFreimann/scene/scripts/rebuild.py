"""Stage a scene-only rebuild from downloaded tiles; no CFD setup or solver.

Use a new workspace. Then run the printed QGIS/Python commands.
Turbine IDs/dimensions in prepare_scene.py are specific to this Munich site.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil


def stage(raw, workspace):
    raw, workspace = raw.resolve(), workspace.resolve()
    if workspace.exists():
        raise FileExistsError(f'Refusing to overwrite {workspace}')
    tiles = sorted(p for p in raw.iterdir() if re.fullmatch(r'\d+_\d+\.(tif|gml)', p.name))
    footprints = raw / '091_Oberbayern_Hausumringe.zip'
    if not any(p.suffix == '.tif' for p in tiles) or not any(p.suffix == '.gml' for p in tiles):
        raise ValueError('Expected terrain TIFFs and LoD2 GML tiles')
    if not footprints.is_file():
        raise FileNotFoundError(footprints)
    scripts = Path(__file__).resolve().parent
    manifest = json.loads((scripts.parent / 'qgis/scene_manifest.json').read_text())
    expected = {Path(s['path']).name: s['sha256'] for s in manifest['sources']}
    expected[footprints.name] = manifest['official_footprints']['sha256']
    changes = []
    for path in [*tiles, footprints]:
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        if digest != expected.get(path.name):
            changes.append(path.name)
    scene = workspace / 'Scene'
    scene.mkdir(parents=True)
    for path in [*tiles, footprints]:
        shutil.copy2(path, workspace / path.name)
    for name in ['prepare_scene.py', 'prepare_footprints.py', 'check_height_conflicts.py', 'setup_scene.py']:
        shutil.copy2(scripts / name, scene / name)
    print(f'Staged {len(tiles)} tiles and the footprint archive in {workspace}')
    print('New or changed source editions/tiles:', ', '.join(changes) if changes else 'none; all original hashes match')
    print('Using scripts/qgis_python.sh (or a compatible QGIS/GDAL Python), run:')
    for name in ['prepare_scene.py', 'prepare_footprints.py', 'check_height_conflicts.py']:
        print(scene / name)
    print('Then run setup_scene.py in a new QGIS project via its Python Console.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--workspace', type=Path, required=True)
    args = parser.parse_args()
    stage(args.raw, args.workspace)
