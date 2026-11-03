"""Package an existing prepared scene for QGIS, without raw-download dependencies.

Run with scripts/qgis_python.sh. Refuses to replace an existing output folder.
Only packages the geographic scene; does not prepare or run CFD.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZIP_DEFLATED
import argparse
import hashlib
import json
import shutil
import sqlite3
import xml.etree.ElementTree as ET

import numpy as np
from osgeo import gdal

gdal.UseExceptions()
WORK = Path(__file__).resolve().parents[1]
PROJECT = 'Frottmaning_Freimann_3D.qgz'
OPEN_GUIDE = '''# Open the Fröttmaning–Freimann scene

Open **Frottmaning_Freimann_3D.qgz** in QGIS 4.2.3. It opens the full 4 km × 3 km area with buildings, bridges and both turbines.

On this Mac, from this folder:

```sh
open -a /Applications/QGIS-final-4_2_3.app Frottmaning_Freimann_3D.qgz
```

Keep this folder together. No raw downloads or wind simulation are needed to view it. There is one full-area 3D view; if closed, reopen **Fröttmaning & Freimann 3D** from **View > 3D Map Views**.

The [site reference image](../crop-reference.png) is an optional illustration for future work; it does not define active QGIS scenes or validated CFD domains. The [Project 3 notes](../../NOTES.md) explain the research plan. Rebuild/verification scripts are in `../scripts/`.

Geobasisdaten: [Bayerische Vermessungsverwaltung](https://www.geodaten.bayern.de/) (Daten verändert), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Data was selected, merged and styled; rotors are illustrative.
'''


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_project(path):
    with ZipFile(path) as archive:
        name = next(n for n in archive.namelist() if n.endswith('.qgs'))
        return name, ET.fromstring(archive.read(name))


def geopackage_content_digest(path):
    """Hash feature data and CRS while ignoring SQLite's mutable journal header."""
    value = hashlib.sha256()
    with sqlite3.connect(path.resolve().as_uri()+'?mode=ro&immutable=1', uri=True) as db:
        for query in ['SELECT table_name,data_type,srs_id FROM gpkg_contents ORDER BY table_name',
                      'SELECT * FROM gpkg_geometry_columns ORDER BY table_name']:
            value.update(json.dumps(db.execute(query).fetchall()).encode())
        for (table,) in db.execute('SELECT table_name FROM gpkg_contents ORDER BY table_name'):
            quoted = '"'+table.replace('"','""')+'"'
            value.update(json.dumps(db.execute('PRAGMA table_info('+quoted+')').fetchall()).encode())
            for row in db.execute('SELECT * FROM '+quoted+' ORDER BY fid'):
                normalized = [v.hex() if isinstance(v,bytes) else v for v in row]
                value.update(json.dumps(normalized,ensure_ascii=True).encode())
    return value.hexdigest()


def verify(folder):
    folder = folder.resolve()
    manifest = json.loads((folder / 'scene_manifest.json').read_text())
    _, root = read_project(folder / PROJECT)
    layers = root.findall('./projectlayers/maplayer')
    for layer in layers:
        source = layer.findtext('datasource', '').split('|')[0]
        path = (folder / source).resolve()
        if Path(source).is_absolute() or not path.is_relative_to(folder) or not path.is_file():
            raise ValueError(f'External or missing QGIS source: {source}')
    for name, expected in manifest['packaged_files_sha256'].items():
        path = folder / name
        if sha256(path) != expected:
            if path.suffix == '.gpkg':
                logical = manifest.get('geopackage_content_sha256',{}).get(name)
                unchanged = logical and geopackage_content_digest(path)==logical
            else:
                unchanged = False
            if not unchanged:
                raise ValueError(f'Changed packaged data: {name}')
        if path.suffix == '.gpkg':
            with sqlite3.connect(path.as_uri() + '?mode=ro&immutable=1', uri=True) as db:
                if db.execute('PRAGMA quick_check').fetchall() != [('ok',)]:
                    raise ValueError(f'Invalid database: {path}')
    terrain = gdal.Open(str(folder / manifest['files']['dem']))
    assert (terrain.RasterXSize, terrain.RasterYSize) == tuple(manifest['dem_size_pixels'])
    assert terrain.GetSpatialRef().GetAuthorityCode(None) == '25832'
    assert terrain.GetGeoTransform() == (manifest['extent_utm32'][0], 1., 0.,
                                        manifest['extent_utm32'][3], 0., -1.)
    views = root.findall('./mapViewDocks3D/view')
    assert any(view.get('isOpen') == '1' for view in views), 'No saved open 3D view'
    print(f'PASS: {len(layers)} local QGIS layers; data hashes, databases, 1 m terrain and saved 3D view.')


def package(source, output):
    source = source.resolve()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    manifest = json.loads((source / 'scene_manifest.json').read_text())
    qgs_name, root = read_project(source / PROJECT)
    files = {
        manifest['files']['dem']: 'terrain.tif',
        manifest['files']['hillshade']: 'hillshade.tif',
        manifest['files']['geometry']: 'geometry.gpkg',
        manifest['official_footprints']['output']: 'footprints.gpkg',
        manifest['height_review']['output']: 'height_review.gpkg',
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='.scene-package-', dir=output.parent) as temp:
        stage = Path(temp)
        dem = gdal.Open(str(source / manifest['files']['dem']))
        packed = gdal.Translate(str(stage / 'terrain.tif'), dem, format='GTiff',
                                creationOptions=['COMPRESS=DEFLATE', 'PREDICTOR=3', 'TILED=YES'])
        assert packed.GetGeoTransform() == dem.GetGeoTransform()
        assert packed.GetProjection() == dem.GetProjection()
        for y in range(0, dem.RasterYSize, 256):
            height = min(256, dem.RasterYSize - y)
            assert np.array_equal(dem.ReadAsArray(0, y, dem.RasterXSize, height),
                                  packed.ReadAsArray(0, y, dem.RasterXSize, height))
        packed = None
        for old, new in files.items():
            if new == 'terrain.tif':
                continue
            if new.endswith('.gpkg'):
                with sqlite3.connect((source / old).as_uri() + '?mode=ro', uri=True) as src:
                    with sqlite3.connect(stage / new) as dst:
                        src.backup(dst)
            else:
                shutil.copy2(source / old, stage / new)

        # Keep the eight prepared layers. The eighteen raw tile layers are redundant.
        raw = root.find("./layer-tree-group/layer-tree-group[@name='Original downloaded tiles']")
        if raw is None:
            raise ValueError('Expected the original prepared-scene layer group')
        removed = {node.get('id') for node in raw.findall('.//layer-tree-layer')}
        for parent in root.iter():
            for child in list(parent):
                if child is raw or child.get('id') in removed or child.text in removed:
                    parent.remove(child)
                elif child.tag == 'maplayer' and child.findtext('id') in removed:
                    parent.remove(child)
        for layer in root.findall('./projectlayers/maplayer'):
            ds = layer.find('datasource')
            old, *options = ds.text.split('|')
            ds.text = './' + files[Path(old).name] + ('|' + '|'.join(options) if options else '')
            if layer.findtext('layername') == 'LoD2 | original roofs and walls':
                for edges in layer.findall('./renderer-3d//edges'):
                    edges.set('enabled', '1')
                    edges.set('color', '69,77,79,255')
                    edges.set('width', '0.65')
        for node in root.findall('.//layer-tree-layer'):
            node.set('source', next(layer.findtext('datasource') for layer in
                     root.findall('./projectlayers/maplayer') if layer.findtext('id') == node.get('id')))
        group = root.find("./layer-tree-group/layer-tree-group[@name='Scene | existing DGM1 coverage']")
        group.set('name', 'Scene | terrain, buildings and turbines')
        # Use a relative manifest path, including in project variables.
        for node in root.iter('value'):
            if (node.text or '').endswith('/scene_manifest.json'):
                node.text = './scene_manifest.json'
        view = root.find('./mapViewDocks3D/view')
        view.set('isOpen', '1')
        view.set('name', 'Fröttmaning & Freimann 3D')
        view.set('width', '1200')
        view.set('height', '800')
        view.set('x', '80')
        view.set('y', '80')
        settings = view.find('qgis3d')
        settings.find('origin').attrib.update(x='695500', y='5343500', z='500')
        settings.find('color').set('background', '201,220,232,255')
        lights = settings.find('lights')
        lights.clear()
        ET.SubElement(lights, 'directional-light', {
            'id': '064b9a44-5be5-455d-add5-41b48d4f0d3d', 'color': '255,244,229,255',
            'intensity': '0.85', 'x': '-0.6', 'y': '-0.4', 'z': '-0.8'})
        # Exclude diagnostic markers and 2D footprints from the 3D geometry list.
        layers3d = settings.find('./terrain/layers')
        layers3d.clear()
        for layer in root.findall('./projectlayers/maplayer'):
            if layer.findtext('layername') in ['LoD2 | original roofs and walls',
                    'Turbine rotors | illustrative, fixed pose', 'Terrain hillshade',
                    'DGM1 | 1 m terrain, metres NH']:
                ET.SubElement(layers3d, 'layer', {'id': layer.findtext('id')})
        view.find('camera').attrib.update(dist='5200', pitch='48', yaw='15',
                                         xMap='696000', yMap='5343500', zMap='525')
        root.attrib.pop('saveDateTime', None)
        with ZipFile(stage / PROJECT, 'w', ZIP_DEFLATED) as archive:
            archive.writestr(qgs_name, ET.tostring(root, encoding='utf-8', xml_declaration=True))
        manifest['files'] = {key: files[name] for key, name in manifest['files'].items()}
        manifest['official_footprints']['output'] = 'footprints.gpkg'
        manifest['height_review']['output'] = 'height_review.gpkg'
        manifest['source_snapshot'] = {'project_sha256': sha256(source / PROJECT),
            'note': 'Source paths are relative to the original prepared-scene directory. Download hashes are provenance, not preview dependencies. Redownloads may be a different edition.'}
        manifest['packaged_files_sha256'] = {name: sha256(stage / name) for name in files.values()}
        manifest['geopackage_content_sha256'] = {name: geopackage_content_digest(stage/name)
            for name in files.values() if name.endswith('.gpkg')}
        manifest['packaging'] = {'terrain': 'All 12,000,000 elevation samples, CRS and transform checked equal to the source VRT; no resampling.',
            'geometry': 'SQLite backups of prepared GeoPackages; no geometry corrections.',
            'display': 'Structure edges enabled; 3D overview and directional light saved; raw tile layers omitted.',
            'wind': 'Historical source scene, wind inputs, exports and results untouched.'}
        (stage / 'scene_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
        (stage / 'README.md').write_text(OPEN_GUIDE)
        verify(stage)
        stage.rename(output)
    print(f'Open {output / PROJECT}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--output', type=Path, default=WORK / 'qgis')
    parser.add_argument('--verify', type=Path, metavar='SCENE_FOLDER')
    args = parser.parse_args()
    if args.verify:
        verify(args.verify)
    elif args.source:
        package(args.source, args.output)
    else:
        parser.error('Supply --source PREPARED_SCENE or --verify PACKAGED_SCENE')
