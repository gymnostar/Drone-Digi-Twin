"""Flag sampled LoD2/terrain intersections without modifying either source."""
from pathlib import Path
from collections import defaultdict
import json
import numpy as np
from scipy.ndimage import map_coordinates
from osgeo import gdal,ogr,osr

gdal.UseExceptions();ogr.UseExceptions()
folder=Path(__file__).resolve().parent
m=json.loads((folder/'scene_manifest.json').read_text())
ds=ogr.Open(str(folder/m['files']['geometry']))
layer=ds.GetLayerByName('lod2_surfaces')
xmin,ymin,xmax,ymax=m['extent_utm32'];boundary_ids=set()
for f in ds.GetLayerByName('building_footprints'):
 lo,hi,south,north=f.GetGeometryRef().GetEnvelope()
 if lo<xmin or hi>xmax or south<ymin or north>ymax:boundary_ids.add(f['building_id'])
raster=gdal.Open(str(folder/m['files']['dem']))
a=raster.ReadAsArray();transform=raster.GetGeoTransform()
roofs={};sample_count=0
for f in layer:
 if f['surface']!='RoofSurface' or f['building_id'] in boundary_ids:continue
 points=[]
 geometry=f.GetGeometryRef()
 for ring in geometry:
  vertices=ring.GetPoints()
  for left,right in zip(vertices,vertices[1:]):
   n=max(1,int(np.ceil(np.hypot(right[0]-left[0],right[1]-left[1]))))
   for t in np.linspace(0,1,n,endpoint=False):
    points.append(tuple(left[j]+(right[j]-left[j])*t for j in range(3)))
  if vertices:points.append(vertices[-1])
 if not points:continue
 points=np.array(points)
 terrain=map_coordinates(a,[(transform[3]-points[:,1])-.5,
                            (points[:,0]-transform[0])-.5],order=1,mode='nearest')
 gaps=terrain-points[:,2];sample_count+=len(points)
 i=int(np.argmax(gaps));bid=f['building_id']
 record={'building_id':bid,'function':f['function'],'overlap_m':float(gaps[i]),
         'x':float(points[i,0]),'y':float(points[i,1]),'source_z_m':float(points[i,2]),
         'terrain_z_m':float(terrain[i]),'source_tile':f['source_tile']}
 if bid not in roofs or record['overlap_m']>roofs[bid]['overlap_m']:roofs[bid]=record
flagged=sorted((r for r in roofs.values() if r['overlap_m']>.3),key=lambda r:r['overlap_m'],reverse=True)
path=folder/f"height_conflicts_{int(m['area_km2'])}km2.gpkg"
if path.exists():raise FileExistsError(path)
ds=None
ds=ogr.GetDriverByName('GPKG').CreateDataSource(str(path));srs=osr.SpatialReference();srs.ImportFromEPSG(25832)
out=ds.CreateLayer('height_review',srs,ogr.wkbPoint)
for name in ['building_id','object_type','source_tile','review_note']:out.CreateField(ogr.FieldDefn(name,ogr.OFTString))
for name in ['overlap_m','source_z_m','terrain_z_m']:out.CreateField(ogr.FieldDefn(name,ogr.OFTReal))
for r in flagged:
 f=ogr.Feature(out.GetLayerDefn());f['building_id']=r['building_id'];f['object_type']='Bridge' if r['function']=='53001_1800' else 'Building/structure'
 for name in ['source_tile','overlap_m','source_z_m','terrain_z_m']:f[name]=r[name]
 f['review_note']='Sampled roof/deck boundary below DGM1 by >0.30 m; review source compatibility. Models and DEM unchanged.'
 point=ogr.Geometry(ogr.wkbPoint);point.AddPoint_2D(r['x'],r['y']);f.SetGeometry(point);out.CreateFeature(f)
ds=None
m['height_review']={'policy':'Preserve exact LoD2 Z and native DGM1 samples. No geometry lifting, flattening or offsets.',
 'method':'Roof/deck boundary vertices and interpolated edge points spaced at most 1 m; bilinear sampling of native DGM1; threshold terrain minus source Z >0.30 m. This is a screening check, not full continuous surface-intersection validation.',
 'objects_checked':len(roofs),'sampled_boundary_points':sample_count,'flagged_objects':len(flagged),
 'flagged_bridges':sum(r['function']=='53001_1800' for r in flagged),
 'max_sampled_overlap_m':max(r['overlap_m'] for r in flagged) if flagged else 0,
 'output':path.name,'objects':flagged}
(folder/'scene_manifest.json').write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:v for k,v in m['height_review'].items() if k!='objects'},indent=2))
