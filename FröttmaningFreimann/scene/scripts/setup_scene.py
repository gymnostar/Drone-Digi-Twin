"""Run inside QGIS Python Console to load and style the prepared scene."""
from pathlib import Path
import json
import qgis.utils
from qgis.core import *
from qgis._3d import *
from qgis.PyQt.QtGui import QColor, QFont, QPainter
from qgis.PyQt.QtWidgets import QDockWidget
from qgis.PyQt import sip

FOLDER = Path(__file__).resolve().parent
MANIFEST = json.loads((FOLDER/'scene_manifest.json').read_text())
PROJECT_PATH = str(FOLDER/'Frottmaning_Freimann_3D.qgz')
project = QgsProject.instance()
root = project.layerTreeRoot()
Qgs3D.initialize()
project.setCrs(QgsCoordinateReferenceSystem('EPSG:25832'))
vertical = QgsCoordinateReferenceSystem('EPSG:7837')
assert vertical.description() == 'DHHN2016 height', vertical.description()
project.setVerticalCrs(vertical)
project.setTitle('Fröttmaning & Freimann | DGM1 + LoD2')
project.setFileName(PROJECT_PATH)
project.setFilePathStorage(Qgis.FilePathType.Relative)
extent = QgsRectangle(*MANIFEST['extent_utm32'])

# Keep every original layer available in a collapsed, hidden group.
raw = root.findGroup('Original downloaded tiles')
if raw is None:
    original_nodes = list(root.children())
    raw = root.addGroup('Original downloaded tiles')
    for node in original_nodes:
        raw.addChildNode(node.clone())
        root.removeChildNode(node)
raw.setItemVisibilityChecked(False)
raw.setExpanded(False)
# New downloads remain available at their complete original extent.
known_sources={layer.source().split('|',1)[0] for layer in project.mapLayers().values()}
for source in MANIFEST['sources']:
    path=(FOLDER/source['path']).resolve()
    if str(path) not in known_sources:
        layer=(QgsRasterLayer(str(path),path.stem) if path.suffix=='.tif' else
               QgsVectorLayer(str(path)+'|layername=Building',path.stem+' — Building','ogr'))
        assert layer.isValid(),str(path)
        project.addMapLayer(layer,False)
        raw.addLayer(layer).setItemVisibilityChecked(False)
scene = root.findGroup('Scene | existing DGM1 coverage')
if scene is None:
    scene = root.insertGroup(0,'Scene | existing DGM1 coverage')


def add(layer):
    assert layer.isValid(), layer.source()
    base_source = layer.source().split('|subset=',1)[0]
    matches = [existing for existing in project.mapLayersByName(layer.name())
               if existing.source().split('|subset=',1)[0] == base_source]
    if matches:
        project.removeMapLayers([existing.id() for existing in matches[1:]])
        return matches[0]
    same_name=[existing for existing in project.mapLayersByName(layer.name())
               if scene.findLayer(existing.id())]
    if same_name:
        existing=same_name[0]
        existing.setDataSource(layer.source(),layer.name(),layer.providerType())
        assert existing.isValid(),layer.source()
        return existing
    project.addMapLayer(layer,False)
    scene.addLayer(layer)
    root.findLayer(layer.id()).setExpanded(False)
    return layer


def load_vector(table,name):
    return add(QgsVectorLayer(str(FOLDER/MANIFEST['files']['geometry'])+'|layername='+table,name,'ogr'))


dem = add(QgsRasterLayer(str(FOLDER/MANIFEST['files']['dem']),'DGM1 | 1 m terrain, metres NH'))
shade = add(QgsRasterLayer(str(FOLDER/MANIFEST['files']['hillshade']),'Terrain hillshade'))
lo,hi = MANIFEST['dem_min_max_m']
ramp = QgsColorRampShader(lo,hi)
ramp.setColorRampType(Qgis.ShaderInterpolationMethod.Linear)
ramp.setColorRampItemList([
    QgsColorRampShader.ColorRampItem(lo,QColor('#e4ead9'),f'{lo:.1f} m'),
    QgsColorRampShader.ColorRampItem(500,QColor('#c5d5b1'),'500 m'),
    QgsColorRampShader.ColorRampItem(535,QColor('#a6bd87'),'535 m'),
    QgsColorRampShader.ColorRampItem(hi,QColor('#e1d6b4'),f'{hi:.1f} m'),
])
shader = QgsRasterShader(lo,hi)
shader.setRasterShaderFunction(ramp)
dem.setRenderer(QgsSingleBandPseudoColorRenderer(dem.dataProvider(),1,shader))
shade_renderer = QgsSingleBandGrayRenderer(shade.dataProvider(),1)
enhancement = QgsContrastEnhancement(shade.dataProvider().dataType(1))
enhancement.setContrastEnhancementAlgorithm(QgsContrastEnhancement.ContrastEnhancementAlgorithm.StretchToMinimumMaximum)
enhancement.setMinimumValue(0)
enhancement.setMaximumValue(255)
shade_renderer.setContrastEnhancement(enhancement)
shade_renderer.setOpacity(.42)
shade.setRenderer(shade_renderer)
shade.setBlendMode(QPainter.CompositionMode.CompositionMode_Multiply)
provider = QgsRasterDemTerrainProvider()
provider.setLayer(dem)
project.elevationProperties().setTerrainProvider(provider)
dem.elevationProperties().setEnabled(True)

surfaces = load_vector('lod2_surfaces','LoD2 | original roofs and walls')
rotors = load_vector('illustrative_rotors','Turbine rotors | illustrative, fixed pose')
footprints = load_vector('building_footprints','Building footprints | derived from LoD2')
sites = load_vector('turbine_locations','Wind turbines | locations and dimensions')
for previous in project.mapLayersByName('Official Hausumringe | Munich building footprints'):
    previous.setName('Official Hausumringe | study area footprints')
official_info=MANIFEST['official_footprints']
official = add(QgsVectorLayer(str(FOLDER/official_info['output'])+'|layername='+official_info['layer'],
                             'Official Hausumringe | study area footprints','ogr'))
official.renderer().setSymbol(QgsFillSymbol.createSimple({
    'style':'no','outline_color':'#1a7591','outline_width':'0.13'}))
official.setCustomProperty('scene/source','Bavarian Hausumringe; study rectangle only; 2D only')
review=None
if 'height_review' in MANIFEST:
    review=add(QgsVectorLayer(str(FOLDER/MANIFEST['height_review']['output'])+'|layername=height_review',
                             'Height review | source terrain intersections','ogr'))
    review.renderer().setSymbol(QgsMarkerSymbol.createSimple({
        'name':'triangle','color':'#ffbf47','outline_color':'#5e3b0c',
        'outline_width':'0.25','size':'3.2'}))
    metadata=review.metadata()
    metadata.setAbstract(MANIFEST['height_review']['method']+' '+MANIFEST['height_review']['policy'])
    review.setMetadata(metadata)

# Hide buildings crossing the terrain boundary rather than placing their outer
# sections over empty space. All candidates remain in the GeoPackage.
footprints.setSubsetString('')
boundary_ids = []
for f in footprints.getFeatures():
    if not extent.contains(f.geometry().boundingBox()):
        boundary_ids.append(f['building_id'])
quoted_ids = ','.join("'"+bid.replace("'","''")+"'" for bid in boundary_ids)
subset = '"building_id" NOT IN ('+quoted_ids+')' if boundary_ids else ''
assert surfaces.setSubsetString(subset)
assert footprints.setSubsetString(subset)
count_request = QgsFeatureRequest().setFlags(Qgis.FeatureRequestFlag.NoGeometry).setNoAttributes()
displayed_buildings = sum(1 for _ in footprints.getFeatures(count_request))
displayed_surfaces = sum(1 for _ in surfaces.getFeatures(count_request))
MANIFEST['counts']['boundary_buildings_hidden'] = len(boundary_ids)
MANIFEST['counts']['displayed_buildings'] = displayed_buildings
MANIFEST['counts']['displayed_surfaces'] = displayed_surfaces
MANIFEST['boundary_policy'] = 'Display complete buildings whose ground footprints fit within the available terrain; boundary candidates remain in the GeoPackage and original tiles.'
(FOLDER/'scene_manifest.json').write_text(json.dumps(MANIFEST,indent=2,ensure_ascii=False)+'\n')


def symbol3d(colour):
    symbol = QgsPolygon3DSymbol()
    symbol.setAltitudeClamping(Qgis.AltitudeClamping.Absolute)
    symbol.setAltitudeBinding(Qgis.AltitudeBinding.Vertex)
    symbol.setHeight(0)
    symbol.setExtrusionHeight(0)
    symbol.setCullingMode(Qgs3DTypes.CullingMode.NoCulling)
    symbol.setAddBackFaces(True)
    symbol.setEdgesEnabled(False)
    material = QgsPhongMaterialSettings()
    material.setDiffuse(QColor(colour))
    material.setAmbient(QColor(colour).lighter(110))
    material.setSpecular(QColor('#333333'))
    material.setShininess(12)
    material.setAmbientCoefficient(.25)
    material.setDiffuseCoefficient(.75)
    symbol.setMaterialSettings(material)
    return symbol


def style_surfaces(layer,key,colours):
    categories = []
    rule_root = QgsRuleBased3DRenderer.Rule(None)
    for value,colour in colours.items():
        fill = QgsFillSymbol.createSimple({'color':colour,'outline_style':'no'})
        categories.append(QgsRendererCategory(value,fill,value))
        rule_root.appendChild(QgsRuleBased3DRenderer.Rule(
            symbol3d(colour),f'"{key}" = \'{value}\'',value))
    layer.setRenderer(QgsCategorizedSymbolRenderer(key,categories))
    renderer = QgsRuleBased3DRenderer(rule_root)
    renderer.setLayer(layer)
    layer.setRenderer3D(renderer)
    layer.elevationProperties().setClamping(Qgis.AltitudeClamping.Absolute)


style_surfaces(surfaces,'surface',{
    'GroundSurface':'#b7ad9c', 'WallSurface':'#e3ded3',
    'RoofSurface':'#b67355', 'ClosureSurface':'#c7bda9',
})
style_surfaces(rotors,'colour',{'white':'#f6f5ef','red':'#b84439'})
# Bridges are road structures, not terracotta building roofs. Preserve their Z.
bridge_symbol=QgsFillSymbol.createSimple({'color':'#7d898e','outline_style':'no'})
bridge_2d=QgsRuleBasedRenderer(bridge_symbol)
bridge_2d.rootRule().children()[0].setFilterExpression('"function" = \'53001_1800\'')
building_rules=QgsRuleBasedRenderer.convertFromRenderer(surfaces.renderer())
for rule in building_rules.rootRule().children():
    rule.setFilterExpression('('+rule.filterExpression()+') AND "function" <> \'53001_1800\'')
    bridge_2d.rootRule().appendChild(rule.clone())
surfaces.setRenderer(bridge_2d)
three_d=surfaces.renderer3D()
for rule in three_d.rootRule().children():
    rule.setFilterExpression('('+rule.filterExpression()+') AND "function" <> \'53001_1800\'')
three_d.rootRule().appendChild(QgsRuleBased3DRenderer.Rule(
    symbol3d('#7d898e'),'"function" = \'53001_1800\'','Bridges — original deck elevations'))
footprints.renderer().setSymbol(QgsFillSymbol.createSimple({
    'color':'#f2c987','outline_color':'#79583e','outline_width':'0.15'}))
sites.renderer().setSymbol(QgsMarkerSymbol.createSimple({
    'name':'circle','color':'#254d65','outline_color':'white',
    'outline_width':'0.4','size':'3.0'}))
label = QgsPalLayerSettings()
label.fieldName = 'name'
fmt = QgsTextFormat()
fmt.setFont(QFont('Arial',10))
fmt.setSize(10)
fmt.setColor(QColor('#244a5c'))
buffer = QgsTextBufferSettings()
buffer.setEnabled(True)
buffer.setColor(QColor('white'))
buffer.setSize(1.1)
fmt.setBuffer(buffer)
label.setFormat(fmt)
sites.setLabeling(QgsVectorLayerSimpleLabeling(label))
sites.setLabelsEnabled(True)

for layer,description in [
    (dem,'Native 1 m mosaic of all downloaded DGM1 tiles. No resampling, flattening or vertical scaling. Elevations in DHHN2016 NH metres.'),
    (surfaces,f"{displayed_buildings} complete buildings within the existing terrain strip. Original XYZ roof/wall/ground rings preserved. Boundary and outside buildings remain available in the GeoPackage/original tiles."),
    (footprints,'XY GroundSurface outlines derived from Bavarian LoD2; not a separate official Hausumringe download.'),
    (rotors,MANIFEST['rotor_limitations']),
    (sites,'Locations from the centre of each original LoD2 tower; hub heights and rotor diameters from SWM. Rotor pose is illustrative.'),
    (official,'Official Bayerische Vermessungsverwaltung Hausumringe, spatially filtered from the downloaded Oberbayern archive to the study area. Original complete 2D polygons, no invented height or extrusion. Visible in the 2D map; no citywide layer loaded.'),
]:
    metadata=layer.metadata()
    metadata.setAbstract(description)
    layer.setMetadata(metadata)

# Use only terrain layers for the ground texture, avoiding building roofs on the ground.
theme = QgsMapThemeCollection.MapThemeRecord()
for layer in [shade,dem]:
    record = QgsMapThemeCollection.MapThemeLayerRecord(layer)
    record.isVisible=True
    theme.addLayerRecord(record)
project.mapThemeCollection().insert('Terrain texture',theme)
for node in scene.findLayers():
    node.setItemVisibilityChecked(node.layer()!=footprints)
scene.setItemVisibilityChecked(True)
root.setHasCustomLayerOrder(True)
root.setCustomLayerOrder(([review] if review else [])+[sites,rotors,surfaces,official,footprints,shade,dem]+
    [n.layer() for n in raw.findLayers()])
project.viewSettings().setDefaultViewExtent(QgsReferencedRectangle(extent,project.crs()))
project.viewSettings().setPresetFullExtent(QgsReferencedRectangle(extent,project.crs()))
iface.mapCanvas().setExtent(extent)
iface.mapCanvas().setCanvasColor(QColor('#edf0ef'))
iface.mapCanvas().refresh()

canvases = iface.mapCanvases3D()
for existing_canvas in canvases:
    sip.transferto(existing_canvas,None)
matching = [c for c in canvases if c.objectName() == 'frottmaning_freimann_scene']
if not matching:
    matching = [c for c in canvases if c.mapSettings() and
                c.mapSettings().terrainSettings().type() == 'dem' and
                c.mapSettings().terrainSettings().layer() == dem]
canvas = matching[0] if matching else iface.createNewMapCanvas3D('Fröttmaning & Freimann 3D')
assert canvas is not None,'Could not create 3D view'
sip.transferto(canvas,None)
qgis.utils._frottmaning_scene_canvas = canvas
canvas.setObjectName('frottmaning_freimann_scene')
settings = canvas.mapSettings()
settings.setCrs(project.crs())
settings.setExtent(extent)
settings.setOrigin(QgsVector3D(695500,5343500,500))
settings.setTransformContext(project.transformContext())
settings.setPathResolver(project.pathResolver())
settings.setMapThemeCollection(project.mapThemeCollection())
terrain = QgsDemTerrainSettings()
terrain.setLayer(dem)
terrain.setVerticalScale(1)
terrain.setElevationOffset(0)
terrain.setResolution(256)
terrain.setMapTileResolution(1024)
terrain.setMaximumGroundError(1)
terrain.setMaximumScreenError(.5)
terrain.setSkirtHeight(10)
settings.setTerrainSettings(terrain)
settings.setTerrainMapTheme('Terrain texture')
settings.setTerrainRenderingEnabled(True)
settings.setTerrainShadingEnabled(True)
settings.setLayers([rotors,surfaces,shade,dem])
settings.setBackgroundColor(QColor('#c9dce8'))
settings.setShowLabels(False)
sun = QgsDirectionalLightSettings()
sun.setDirection(QgsVector3D(-0.6,-0.4,-0.8))
sun.setIntensity(.85)
sun.setColor(QColor('#fff4e5'))
settings.setLightSources([sun])
settings.setFieldOfView(42)
pose = QgsCameraPose()
pose.setCenterPoint(settings.mapToWorldCoordinates(QgsVector3D(extent.center().x(),extent.center().y(),525)))
pose.setDistanceFromCenterPoint(4200)
pose.setPitchAngle(48)
pose.setHeadingAngle(15)
canvas.cameraController().setCameraPose(pose)

metadata = project.metadata()
metadata.setTitle('Fröttmaning & Freimann — Bavarian DGM1 and LoD2')
metadata.setAbstract(f"{MANIFEST['area_km2']:g} km² scene using {len([s for s in MANIFEST['sources'] if s['path'].endswith('.tif')])} existing 1 m DGM1 tiles and {displayed_buildings} of {MANIFEST['counts']['source_buildings']} LoD2 objects, including bridges. "+
    'Horizontal CRS EPSG:25832; height datum DHHN2016. No vertical exaggeration. '+
    MANIFEST['rotor_limitations']+' Original layers retained. Official Hausumringe added for the study area only. No DLR data or additional terrain downloaded by the agent. Actual terrain/structure intersections remain at original heights for simulation review.')
project.setMetadata(metadata)
QgsExpressionContextUtils.setProjectVariable(project,'scene_manifest',str(FOLDER/'scene_manifest.json'))
QgsExpressionContextUtils.setProjectVariable(project,'rotor_geometry','illustrative_fixed_pose')
for dock in iface.mainWindow().findChildren(QDockWidget):
    if dock.objectName()=='PythonConsole':
        dock.hide()
iface.setActiveLayer(sites)
assert project.write(PROJECT_PATH)
print('Saved',PROJECT_PATH)
print('Scene features:',displayed_surfaces,'LoD2 faces;',displayed_buildings,'objects;',sites.featureCount(),'turbines')
print('Terrain:',dem.width(),dem.height(),'pixels;',settings.terrainSettings().type(),'; scale',settings.terrainSettings().verticalScale())
