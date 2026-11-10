"""Build the standard DLR wind view in ParaView: hub-height slice coloured by wind speed,
arrows for direction, buildings/terrain as solid walls.
GUI:   double-click "Open Wind View.command"   (runs: paraview --script=scripts/wind_view.py)
Batch: pvpython scripts/wind_view.py out.png   (saves a screenshot instead)"""
import os, sys
from paraview.simple import *

# The ParaView GUI runs --script without __file__; the .command cd's into DLR/, so fall back to cwd.
try:
    here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    here = os.path.join(os.getcwd(), 'scripts')
case = os.path.join(here, '..', 'windCase', 'dlr.foam')
SLICE_Z = 35.0          # slice height in metres (local frame; ground is ~12-25 m)

r = OpenFOAMReader(registrationName='DLR wind (dlr.foam)', FileName=case)
r.MeshRegions = ['internalMesh', 'patch/buildings', 'patch/terrain']
r.CellArrays = ['U', 'p', 'k']
r.UpdatePipelineInformation()
last = r.TimestepValues[-1]
scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps(); scene.AnimationTime = last

v = GetActiveViewOrCreate('RenderView')
air   = ExtractBlock(registrationName='air', Input=r);   air.Selectors = ['/Root/internalMesh']
walls = ExtractBlock(registrationName='buildings + terrain', Input=r); walls.Selectors = ['/Root/patch/buildings', '/Root/patch/terrain']
sl = Slice(registrationName='slice at %d m' % SLICE_Z, Input=air)
sl.SliceType.Origin = [0, 0, SLICE_Z]; sl.SliceType.Normal = [0, 0, 1]
arrows = Glyph(registrationName='wind arrows', Input=sl, GlyphType='Arrow')
arrows.OrientationArray = ['CELLS', 'U']; arrows.ScaleArray = ['POINTS', 'No scale array']
arrows.ScaleFactor = 30
# a fixed number of evenly spread arrows, independent of how many cells the mesh has
arrows.GlyphMode = 'Uniform Spatial Distribution (Bounds Based)'; arrows.MaximumNumberOfSamplePoints = 1500

w = Show(walls, v); w.ColorArrayName = ['POINTS', '']; w.DiffuseColor = [0.85, 0.55, 0.4]
d = Show(sl, v); ColorBy(d, ('CELLS', 'U', 'Magnitude'))
lut = GetColorTransferFunction('U')
try: lut.ApplyPreset('Viridis', True)
except Exception: pass
lut.RescaleTransferFunction(0.0, 7.0)
d.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lut, v); bar.Title = 'Wind speed (m/s)'; bar.ComponentTitle = ''
a = Show(arrows, v); a.ColorArrayName = ['POINTS', '']; a.DiffuseColor = [0.1, 0.1, 0.1]

v.ViewTime = last
Render()   # first render resets the camera, so set ours afterwards
# look down at the campus from the south-west (domain centre = campus centre)
v.CameraFocalPoint = [60, 60, 25]; v.CameraPosition = [-300, -800, 900]; v.CameraViewUp = [0, 0, 1]
v.CameraViewAngle = 30
SetActiveSource(sl)
try: Hide3DWidgets(proxy=sl.SliceType)   # hide the red plane handle; untick/tick 'Show Plane' to get it back
except Exception: pass
Render()

if len(sys.argv) > 1 and sys.argv[1].endswith('.png'):
    try: v.UseColorPaletteForBackground = 0
    except Exception: pass
    v.Background = [1, 1, 1]
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]; bar.RangeLabelFormat = '%.0f'
    SaveScreenshot(sys.argv[1], v, ImageResolution=[1600, 900]); print('saved', sys.argv[1])
