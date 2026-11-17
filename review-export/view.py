#!/usr/bin/env pvpython
"""Export saved CFD volumes and open portable ParaView views; no solver execution."""
import json
import os
from pathlib import Path
import sys
from paraview.simple import *

# The ParaView GUI executes --script without defining __file__.
ROOT=(Path(os.environ['DRONE_VIEW_DIR']).resolve() if 'DRONE_VIEW_DIR' in os.environ
      else Path(globals().get('__file__', 'view.py')).resolve().parent)
PROJECT=ROOT.parent
SCENES={
    'joint-baseline':('scenes/joint/runs/poc-N-baseline-common',True),
    'turbine-cheap':('scenes/turbine/runs/poc-studio-development-strict-cheap',False),
    'turbine-reference':('scenes/turbine/runs/poc-studio-development-reference',False),
}

def build(name):
    ResetSession()
    relative,joint=SCENES[name]
    case=PROJECT/relative
    marker=case/'view.foam';marker.touch(exist_ok=True)
    reader=OpenFOAMReader(FileName=str(marker))
    reader.MeshRegions=['internalMesh'];reader.CellArrays=['U','k','p']
    reader.UpdatePipelineInformation();iteration=max(reader.TimestepValues)
    reader.UpdatePipeline(time=iteration)
    merged=MergeBlocks(Input=reader);merged.UpdatePipeline(time=iteration)
    count=merged.GetDataInformation().GetNumberOfCells()
    velocity=merged.CellData['U'].GetRange(-1)
    volume=ROOT/f'{name}.vtu'
    SaveData(str(volume),proxy=merged,ChooseArraysToWrite=1,CellDataArrays=['U','k','p'],FieldDataArrays=[],CompressorType='ZLib',DataMode='Appended')
    Delete(merged);Delete(reader)
    flow=XMLUnstructuredGridReader(registrationName='Saved CFD volume',FileName=[str(volume)])
    flow.UpdatePipeline()
    points=CellDatatoPointData(registrationName='Velocity interpolation for display',Input=flow)
    points.PassCellData=1
    clip=Clip(registrationName='Display region',Input=points,ClipType='Box')
    clip.ClipType.Position=[-220,-470,0] if joint else [-55,-40,0]
    clip.ClipType.Length=[440,690,110] if joint else [210,80,80]
    clip.Invert=1
    view=CreateView('RenderView');view.ViewSize=[1440,1000]
    view.Background=[.96,.97,.99];view.UseColorPaletteForBackground=0
    view.OrientationAxesVisibility=1
    horizontal=Slice(registrationName='Velocity at rotor height',Input=clip)
    horizontal.SliceType.Origin=[0,0,32.1];horizontal.SliceType.Normal=[0,0,1]
    vertical=Slice(registrationName='Vertical wake section',Input=clip)
    vertical.SliceType.Origin=[0,0,0];vertical.SliceType.Normal=[1,0,0] if joint else [0,1,0]
    for item,opacity in ((horizontal,.30),(vertical,.40)):
        display=Show(item,view);ColorBy(display,('POINTS','U','Magnitude'));display.Opacity=opacity
        display.SetScalarBarVisibility(view,True)
    lut=GetColorTransferFunction('U');lut.ApplyPreset('Cool to Warm',True)
    lut.RescaleTransferFunction(0,12)
    bar=GetScalarBar(lut,view);bar.Title='Mean speed';bar.ComponentTitle='m/s'
    bar.TitleColor=[.08,.10,.14];bar.LabelColor=[.08,.10,.14]
    stream=StreamTracer(registrationName='Steady mean-flow streamlines',Input=clip,SeedType='Line')
    stream.Vectors=['POINTS','U'];stream.IntegrationDirection='FORWARD'
    stream.SeedType.Point1=[-175,185,40] if joint else [-45,-28,32.1]
    stream.SeedType.Point2=[175,185,40] if joint else [-45,28,32.1]
    stream.SeedType.Resolution=48
    stream.MaximumStreamlineLength=1100 if joint else 230
    tubes=Tube(registrationName='Flow paths',Input=stream);tubes.Radius=.65 if joint else .18
    tube_display=Show(tubes,view);ColorBy(tube_display,('POINTS','U','Magnitude'))
    if joint:
        surface=ROOT/'buildings.stl'
        original=STLReader(FileNames=[str(case/'constant/triSurface/buildings.stl')])
        SaveData(str(surface),proxy=original,FileType='Binary');Delete(original)
        buildings=STLReader(registrationName='Simulated building surfaces',FileNames=[str(surface)])
        display=Show(buildings,view);display.ColorArrayName=['POINTS',''];display.DiffuseColor=[.78,.79,.82]
        ground=Plane(registrationName='Ground')
        ground.Origin=[-220,-470,0];ground.Point1=[220,-470,0];ground.Point2=[-220,220,0]
    else:
        ground=Plane(registrationName='Ground')
        ground.Origin=[-55,-40,0];ground.Point1=[155,-40,0];ground.Point2=[-55,40,0]
    display=Show(ground,view);display.ColorArrayName=['POINTS',''];display.DiffuseColor=[.88,.88,.87]
    rotor=Disk(registrationName='Actuator disk location')
    rotor.InnerRadius=13.1;rotor.OuterRadius=13.5;rotor.CircumferentialResolution=64
    transform=Transform(registrationName='Rotor source marker',Input=rotor)
    transform.Transform.Rotate=[90,0,0] if joint else [0,90,0]
    transform.Transform.Translate=[0,-280,32.1] if joint else [0,0,32.1]
    display=Show(transform,view);display.ColorArrayName=['POINTS',''];display.DiffuseColor=[.92,.35,.12];display.Opacity=.95
    label=Text(registrationName='Simulation description')
    label.Text='Coupled city + turbine | saved steady CFD' if joint else f'Isolated turbine | {name.split("-")[1]} mesh | saved steady CFD'
    display=Show(label,view);display.Color=[.08,.10,.14];display.FontSize=16
    view.CameraPosition=[650,-900,700] if joint else [210,-230,150]
    view.CameraFocalPoint=[0,-110,25] if joint else [50,0,28]
    view.CameraViewUp=[0,0,1]
    view.CameraParallelProjection=1;view.CameraParallelScale=410 if joint else 105
    SetActiveSource(flow);SetActiveView(view)
    lut.RescaleTransferFunction(0,12);lut.AutomaticRescaleRangeMode='Never'
    Render(view)
    state=ROOT/f'{name}.pvsm'
    SaveState(str(state))
    state.write_text(state.read_text().replace(str(ROOT)+'/', ''))
    if '--qa' in sys.argv:
        SaveScreenshot(f'/tmp/{name}-3d-qa.png',view,ImageResolution=[1440,1000])
    return {'source_case':relative,'saved_SIMPLE_iteration':iteration,'cells':count,'U_magnitude_range_m_s':list(velocity),'volume_file':volume.name,'volume_bytes':volume.stat().st_size,'display':'Actual saved velocity fields; slices and streamlines interpolate the steady mean flow. Source disk is a location marker, not blade geometry.'}

if '--build' in sys.argv:
    checks={name:build(name) for name in SCENES}
    (ROOT/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))
else:
    name=os.environ.get('DRONE_VIEW','joint-baseline')
    if name not in SCENES:raise ValueError('Choose joint-baseline, turbine-cheap or turbine-reference')
    LoadState(str(ROOT/f'{name}.pvsm'),data_directory=str(ROOT),restrict_to_data_directory=True)
    RenderAllViews()
    print(f"Loaded 3D simulation: {name}", flush=True)
