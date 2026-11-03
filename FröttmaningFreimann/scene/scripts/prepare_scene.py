"""Prepare the existing Bavarian tiles; run with QGIS's bundled Python/GDAL."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()
OUT = Path(__file__).resolve().parent
DATA = OUT.parent
EXTENT = None
NS = {
    'b': 'http://www.opengis.net/citygml/building/1.0',
    'g': 'http://www.opengis.net/gml',
    'gen': 'http://www.opengis.net/citygml/generics/1.0',
}
GML_ID = '{http://www.opengis.net/gml}id'
TURBINES = {
    'DEBY_LOD2_5149447': ('Fröttmaning', 66.8, 66.0, 'Enercon E-66'),
    'DEBY_LOD2_108248982': ('Freimann', 80.0, 138.0, 'Enercon E-138'),
}


def field(layer, name, kind=ogr.OFTString):
    layer.CreateField(ogr.FieldDefn(name, kind))


def feature(layer, geometry, values):
    f = ogr.Feature(layer.GetLayerDefn())
    for name, value in values.items():
        if value is not None:
            f.SetField(name, value)
    f.SetGeometry(geometry)
    layer.CreateFeature(f)


def polygon(rings, z=True):
    p = ogr.Geometry(ogr.wkbPolygon25D if z else ogr.wkbPolygon)
    for points in rings:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for point in points:
            ring.AddPoint(*point) if z else ring.AddPoint_2D(*point[:2])
        if points[0] != points[-1]:
            ring.AddPoint(*points[0]) if z else ring.AddPoint_2D(*points[0][:2])
        p.AddGeometry(ring)
    return p


def read_polygon(element):
    rings = []
    for boundary in ['exterior', 'interior']:
        for element_ring in element.findall('g:' + boundary + '/g:LinearRing', NS):
            pos = element_ring.find('g:posList', NS)
            if pos is None:
                raise ValueError('Unsupported ring without posList')
            values = list(map(float, pos.text.split()))
            if len(values) % 3:
                raise ValueError('Expected XYZ coordinates')
            rings.append([tuple(values[i:i + 3]) for i in range(0, len(values), 3)])
    if not rings:
        raise ValueError('Polygon has no rings')
    return rings


def main():
    global EXTENT
    tifs = sorted(p for p in DATA.glob('*.tif') if re.fullmatch(r'\d+_\d+\.tif',p.name))
    if not tifs:
        raise ValueError('No DGM1 tiles found')
    bounds = []
    for path in tifs:
        raster = gdal.Open(str(path))
        transform = raster.GetGeoTransform()
        if transform[1:3] != (1.,0.) or transform[4:] != (0.,-1.):
            raise ValueError('Expected aligned native 1 m DGM1 tiles: '+str(path))
        sr = osr.SpatialReference(wkt=raster.GetProjection())
        if sr.GetAuthorityCode(None) != '25832':
            raise ValueError('Expected EPSG:25832: '+str(path))
        bounds.append((transform[0],transform[3]-raster.RasterYSize,
                       transform[0]+raster.RasterXSize,transform[3]))
    EXTENT = (min(b[0] for b in bounds),min(b[1] for b in bounds),
              max(b[2] for b in bounds),max(b[3] for b in bounds))
    area_id = '_'.join(str(int(value/1000)) for value in EXTENT)
    files = {'dem':f'dgm1_{area_id}.vrt','hillshade':f'hillshade_{area_id}.tif',
             'geometry':f'scene_geometry_{area_id}.gpkg'}
    gpkg = OUT / files['geometry']
    if gpkg.exists():
        raise FileExistsError('Prepared geometry already exists: '+str(gpkg))
    vrt_path = OUT / files['dem']
    vrt = gdal.BuildVRT(str(vrt_path), list(map(str, tifs)),
                        options=gdal.BuildVRTOptions(resolution='highest'))
    vrt = None
    xml = ET.parse(vrt_path)
    for source in xml.findall('.//SourceFilename'):
        source.text = '../' + Path(source.text).name
        source.set('relativeToVRT', '1')
    xml.write(vrt_path, encoding='UTF-8', xml_declaration=True)
    dem = gdal.Open(str(vrt_path))
    if (dem.ReadAsArray() == -9999).any():
        raise ValueError('The available DGM1 tiles do not form a complete rectangle')
    stats = dem.GetRasterBand(1).ComputeStatistics(False)
    gdal.DEMProcessing(str(OUT / files['hillshade']), dem, 'hillshade',
                       options=gdal.DEMProcessingOptions(
                           multiDirectional=True, computeEdges=True,
                           creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']))
    ds = ogr.GetDriverByName('GPKG').CreateDataSource(str(gpkg))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(25832)
    surfaces = ds.CreateLayer('lod2_surfaces', srs, ogr.wkbPolygon25D)
    for name in ['building_id', 'surface', 'source_tile', 'function']:
        field(surfaces, name)
    for name in ['ground_m', 'roof_m', 'height_m']:
        field(surfaces, name, ogr.OFTReal)
    footprints = ds.CreateLayer('building_footprints', srs, ogr.wkbMultiPolygon)
    for name in ['building_id', 'source_tile']:
        field(footprints, name)
    field(footprints, 'height_m', ogr.OFTReal)
    rotors = ds.CreateLayer('illustrative_rotors', srs, ogr.wkbPolygon25D)
    for name in ['name', 'part', 'colour', 'geometry_note']:
        field(rotors, name)
    points = ds.CreateLayer('turbine_locations', srs, ogr.wkbPoint25D)
    for name in ['name', 'model', 'lod2_id', 'dimensions_source']:
        field(points, name)
    for name in ['hub_m', 'rotor_diameter_m', 'tip_height_m', 'base_elevation_m']:
        field(points, name, ogr.OFTReal)
    ds.StartTransaction()
    counts = Counter()
    turbine_data = []
    gmls = sorted(DATA.glob('*.gml'))
    for path in gmls:
        root = ET.parse(path).getroot()
        for b in root.findall('.//b:Building', NS):
            counts['source_buildings'] += 1
            bid = b.attrib[GML_ID]
            coords = [float(v) for pos in b.findall('.//g:posList', NS)
                      for v in pos.text.split()]
            if not coords:
                raise ValueError('Building without coordinates: ' + bid)
            xmin, xmax = min(coords[0::3]), max(coords[0::3])
            ymin, ymax = min(coords[1::3]), max(coords[1::3])
            if xmax < EXTENT[0] or xmin > EXTENT[2] or ymax < EXTENT[1] or ymin > EXTENT[3]:
                continue
            counts['scene_buildings'] += 1
            attrs = {a.attrib['name']: a.findtext('gen:value', namespaces=NS)
                     for a in b.findall('gen:stringAttribute', NS)}
            ground = float(attrs['HoeheGrund'])
            roof = float(attrs['HoeheDach'])
            height = float(b.findtext('b:measuredHeight', namespaces=NS))
            values = dict(building_id=bid, source_tile=path.name,
                          function=b.findtext('b:function', namespaces=NS),
                          ground_m=ground, roof_m=roof, height_m=height)
            fp = ogr.Geometry(ogr.wkbMultiPolygon)
            for surface_type in ['GroundSurface', 'WallSurface', 'RoofSurface', 'ClosureSurface']:
                for surface in b.findall('.//b:boundedBy/b:' + surface_type, NS):
                    for element in surface.findall('.//g:Polygon', NS):
                        rings = read_polygon(element)
                        feature(surfaces, polygon(rings), dict(values, surface=surface_type))
                        counts[surface_type] += 1
                        if surface_type == 'GroundSurface':
                            fp.AddGeometry(polygon(rings, z=False))
            if fp.GetGeometryCount():
                feature(footprints, fp, dict(building_id=bid, source_tile=path.name,
                                             height_m=height))
                counts['footprints'] += 1
            if bid in TURBINES:
                name, hub, diameter, model = TURBINES[bid]
                turbine_data.append(dict(name=name, x=(xmin+xmax)/2, y=(ymin+ymax)/2,
                                         ground=ground, hub=hub, diameter=diameter,
                                         model=model, lod2_id=bid))
        print(path.name, 'processed', flush=True)
    if len(turbine_data) != 2:
        raise ValueError('Both source turbine towers must be present')
    for t in turbine_data:
        x, y, base, hub, radius = t['x'], t['y'], t['ground'], t['hub'], t['diameter']/2
        point = ogr.Geometry(ogr.wkbPoint25D)
        point.AddPoint(x, y, base)
        feature(points, point, dict(name=t['name'], model=t['model'], lod2_id=t['lod2_id'],
                                    hub_m=hub, rotor_diameter_m=radius*2,
                                    tip_height_m=hub+radius, base_elevation_m=base,
                                    dimensions_source='https://www.swm.de/unternehmen/magazin/energie/windrad-muenchen'))
        # Fixed north-facing display pose. Blade profile and hub are illustrative.
        def face(vertices, part, colour='white'):
            feature(rotors, polygon([vertices]), dict(
                name=t['name'], part=part, colour=colour,
                geometry_note='Illustrative; fixed north-facing pose; dimensions from SWM; centre/base from LoD2'))
        cy = y - 3.0
        cz = base + hub
        for blade in range(3):
            angle = blade * math.tau / 3
            def xyz(r, width, depth):
                return (x + r*math.sin(angle) + width*math.cos(angle), cy+depth,
                        cz + r*math.cos(angle) - width*math.sin(angle))
            sections = [(radius*f, radius*w) for f,w in
                        [(0.035,0.009),(0.15,0.029),(0.35,0.023),
                         (0.65,0.014),(0.88,0.006),(1.0,0.0007)]]
            for i, ((r0,w0),(r1,w1)) in enumerate(zip(sections, sections[1:])):
                vertices = [xyz(r0,-w0,-.22),xyz(r0,w0,-.22),
                            xyz(r1,w1,-.22),xyz(r1,-w1,-.22),
                            xyz(r0,-w0,.22),xyz(r0,w0,.22),
                            xyz(r1,w1,.22),xyz(r1,-w1,.22)]
                colour = 'red' if t['name']=='Freimann' and i==4 else 'white'
                for ids in [(0,1,2,3),(7,6,5,4),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)]:
                    face([vertices[j] for j in ids],f'blade_{blade+1}',colour)
        r = 1.7 if t['name']=='Freimann' else 1.1
        for j in range(8):
            def sphere(i, k):
                phi=math.pi*(j+i)/8
                theta=math.tau*k/16
                return (x+r*math.sin(phi)*math.cos(theta),
                        cy+r*math.sin(phi)*math.sin(theta),cz+r*math.cos(phi))
            for k in range(16):
                face([sphere(0,k),sphere(0,k+1),sphere(1,k+1),sphere(1,k)],'hub')
        t['terrain_sample_m'] = float(dem.ReadAsArray(
            int(x-EXTENT[0]), int(EXTENT[3]-y), 1, 1)[0,0])
    ds.CommitTransaction()
    ds = None
    manifest = {
        'crs': 'EPSG:25832', 'vertical_datum': 'DHHN2016 NH (metres)',
            'extent_utm32': EXTENT,
            'area_km2': (EXTENT[2]-EXTENT[0])*(EXTENT[3]-EXTENT[1])/1000000,
            'dem_resolution_m': 1,
        'dem_size_pixels': [dem.RasterXSize,dem.RasterYSize],
        'dem_min_max_m': stats[:2], 'counts': dict(counts),
        'selection': 'Complete buildings whose XY bounding boxes intersect existing DGM1 coverage; no source XYZ modified.',
        'turbines': turbine_data,
        'rotor_limitations': 'Illustrative blade/hub profiles; fixed north-facing orientation, not surveyed or animated. Existing LoD2 towers retained.',
        'sources': [{'path': '../'+p.name,'bytes':p.stat().st_size,
                     'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in tifs+gmls],
        'missing_dgm1_tiles_for_full_lod2_grid': [f'{e}_{n}.tif'
            for e in range(min(int(p.stem.split('_')[0]) for p in gmls),
                           max(int(p.stem.split('_')[0]) for p in gmls)+2)
            for n in range(min(int(p.stem.split('_')[1]) for p in gmls),
                           max(int(p.stem.split('_')[1]) for p in gmls)+2)
            if not (DATA/f'{e}_{n}.tif').exists()],
        'files': files,
        'dlr_included': False,
    }
    (OUT/'scene_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(manifest,indent=2,ensure_ascii=False))


if __name__ == '__main__':
    main()
