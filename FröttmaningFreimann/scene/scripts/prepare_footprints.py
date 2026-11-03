"""Extract original official footprint polygons intersecting the study area."""
from pathlib import Path
import hashlib,json
from osgeo import gdal,ogr

gdal.UseExceptions();ogr.UseExceptions()
folder=Path(__file__).resolve().parent
archive=folder.parent/'091_Oberbayern_Hausumringe.zip'
manifest_path=folder/'scene_manifest.json'
m=json.loads(manifest_path.read_text())
output=folder/('area_footprints_'+'_'.join(str(int(v/1000)) for v in m['extent_utm32'])+'.gpkg')
source='/vsizip/'+str(archive)+'/hausumringe.shp'
if output.exists():raise FileExistsError(output)
ds=gdal.VectorTranslate(str(output),source,options=gdal.VectorTranslateOptions(
 format='GPKG',layerName='area_footprints',spatFilter=m['extent_utm32'],
 layerCreationOptions=['SPATIAL_INDEX=YES']))
ds=None
city_ds=ogr.Open(str(output));layer=city_ds.GetLayer(0)
assert layer.GetSpatialRef().GetAuthorityCode(None)=='25832'
count=layer.GetFeatureCount();assert count>0
m['official_footprints']={'provider':'Bayerische Vermessungsverwaltung',
 'download':archive.name,'regional_feature_count':2703736,
 'selection':'Complete original polygons intersecting the study rectangle; no citywide layer loaded.',
 'extent_filter_utm32':m['extent_utm32'],'feature_count':count,
 'output':output.name,'layer':'area_footprints',
 'geometry':'Original 2D polygons; no height or extrusion added.',
 'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}
manifest_path.write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n')
print('STUDY AREA FOOTPRINTS',count,output)
