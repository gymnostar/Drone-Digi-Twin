"""Make a portable inspection copy of the retired pilot; never run a solver.

Run with ParaView 6.2 pvpython, --source PATH_TO_OLD_CFD_PILOT.
Writes display data and a two-view state beside this script.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from paraview.simple import *


def export_data(source, output):
    data = output / 'data'
    data.mkdir(exist_ok=True)
    mesh = np.load(source / 'environment_mesh.npz')
    counts = {}
    for region in ['terrain', 'buildings', 'bridges', 'turbine_structure']:
        faces = mesh['faces'][mesh['regions'] == region]
        indices, inverse = np.unique(faces, return_inverse=True)
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(mesh['vertices'][indices], deep=True))
        cells = vtk.vtkCellArray()
        offsets = np.arange(0, 3 * (len(faces) + 1), 3, dtype=np.int64)
        cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True),
                      numpy_to_vtkIdTypeArray(inverse, deep=True))
        poly = vtk.vtkPolyData()
        poly.SetPoints(points)
        poly.SetPolys(cells)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(str(data / (region + '.vtp')))
        writer.SetInputData(poly)
        writer.SetDataModeToAppended()
        writer.SetCompressorTypeToZLib()
        assert writer.Write() == 1
        counts[region] = poly.GetNumberOfCells()
    disk = vtk.vtkSTLReader()
    disk.SetFileName(str(source / 'case/constant/triSurface/actuator_disk_reference.stl'))
    disk.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(disk.GetOutput())
    writer.SetFileName(str(data / 'disk.vtp'))
    assert writer.Write() == 1
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(source / 'case/VTK/case_750/internal.vtu'))
    reader.Update()
    grid = reader.GetOutput()
    for attrs in [grid.GetCellData(), grid.GetPointData()]:
        for i in reversed(range(attrs.GetNumberOfArrays())):
            if attrs.GetArrayName(i) not in ['U', 'p', 'k', 'omega']:
                attrs.RemoveArray(i)
    assert grid.GetNumberOfCells() == 196568
    for name in ['U', 'p', 'k', 'omega']:
        assert np.isfinite(vtk_to_numpy(grid.GetCellData().GetArray(name))).all()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(data / 'flow.vtu'))
    writer.SetInputData(grid)
    writer.SetDataModeToAppended()
    writer.SetCompressorTypeToZLib()
    assert writer.Write() == 1
    config = json.loads((source / 'pilot_config.json').read_text())
    # Only the assumptions needed to interpret this review, with no old paths.
    metadata = {name: config[name] for name in ['qgis_extent_xy_m', 'local_extent_xy_m',
                'block_z_local_m', 'disk', 'numerical_test', 'runtime_limits']}
    metadata.update(status='Retired numerical pilot; scope review only',
        origin_utm32_dhhn2016_m=[696000, 5343500, 500],
        wind_from_degrees=270, inflow_m_s=5.4,
        model='steady simpleFoam / kOmegaSST / actuator disk',
        candidates_in_source_area=229, volume_cells=grid.GetNumberOfCells(),
        surface_triangles=counts,
        limitations=['54 underdetermined and 6704 concave air cells in historical extended check.',
                    'No measurement/benchmark, grid or domain independence validation.',
                    'Uniform inflow, smooth walls and Ct/Cp are provisional.',
                    'No drone, route planner, live yaw, gusts or shelter classification.'],
        data_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data.iterdir())})
    (output / 'assumptions.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return config


def make_view(output, config):
    layout = CreateLayout('Retired pilot review')
    geometry = CreateView('RenderView')
    flow_view = CreateView('RenderView')
    AssignViewToLayout(view=geometry, layout=layout, hint=0)
    layout.SplitHorizontal(0, .5)
    AssignViewToLayout(view=flow_view, layout=layout, hint=2)
    layout.SetSize(1800, 950)
    for view in [geometry, flow_view]:
        view.UseColorPaletteForBackground = 0
        view.Background = [.96, .97, .98]
        view.CameraViewUp = [0, 0, 1]
        view.CameraParallelProjection = 1
    colors = {'terrain': [.62, .70, .57], 'buildings': [.78, .40, .24],
              'bridges': [.30, .39, .45], 'turbine_structure': [.90, .90, .84]}
    for region, color in colors.items():
        reader = XMLPolyDataReader(registrationName=region, FileName=[str(output / 'data' / (region + '.vtp'))])
        for view in [geometry, flow_view]:
            display = Show(reader, view)
            display.ColorArrayName = ['POINTS', '']
            display.DiffuseColor = color
            display.Ambient = .25
            if region != 'terrain':
                display.Representation = 'Surface With Edges'
                display.EdgeColor = [.20, .23, .25]
                display.LineWidth = .5
            if view == flow_view:
                display.Opacity = .32
    disk = XMLPolyDataReader(registrationName='Rotor disk: reference only', FileName=[str(output / 'data/disk.vtp')])
    for view in [geometry, flow_view]:
        display = Show(disk, view)
        display.ColorArrayName = ['POINTS', '']
        display.DiffuseColor = [.90, .08, .12]
        display.Opacity = .55
    flow = XMLUnstructuredGridReader(registrationName='Retired iteration 750: 3D air field', FileName=[str(output / 'data/flow.vtu')])
    outline = Outline(Input=flow)
    display = Show(outline, geometry)
    display.ColorArrayName = ['POINTS', '']
    display.DiffuseColor = [.45, .48, .52]
    hx, hy, hz = config['disk']['centre_local_m']
    section = Slice(registrationName='Vertical wind slice through hub', Input=flow)
    section.SliceType = 'Plane'
    section.SliceType.Origin = [hx, hy, hz]
    section.SliceType.Normal = [0, 1, 0]
    display = Show(section, flow_view)
    ColorBy(display, ('CELLS', 'U', 'Magnitude'))
    lut = GetColorTransferFunction('U')
    lut.RGBPoints = [0, .267, .005, .329, 2, .23, .322, .546,
                     4, .128, .567, .551, 6, .369, .789, .383, 8, .993, .906, .144]
    lut.RescaleTransferFunction(0, 8)
    display.SetScalarBarVisibility(flow_view, True)
    bar = GetScalarBar(lut, flow_view)
    bar.Orientation = 'Horizontal'
    bar.WindowLocation = 'Lower Center'
    bar.ScalarBarLength = .65
    bar.TitleFontSize = 12
    bar.LabelFontSize = 11
    bar.AutomaticLabelFormat = 0
    bar.LabelFormat = '{:.1f}'
    bar.RangeLabelFormat = '{:.1f}'
    bar.Title = 'Mean speed (m/s)'
    bar.ComponentTitle = ''
    bar.TitleColor = [.1, .1, .1]
    bar.LabelColor = [.1, .1, .1]
    centers = CellCenters(Input=section)
    arrows = Glyph(registrationName='Velocity directions, not flight paths', Input=centers, GlyphType='Arrow')
    arrows.OrientationArray = ['POINTS', 'U']
    arrows.ScaleArray = ['POINTS', 'No scale array']
    arrows.ScaleFactor = 18
    arrows.GlyphMode = 'Every Nth Point'
    arrows.Stride = 18
    display = Show(arrows, flow_view)
    display.ColorArrayName = ['POINTS', '']
    display.DiffuseColor = [.12, .12, .15]
    labels = [(geometry, 'SIMPLIFIED GEOMETRY | E-66 only\nBuildings orange; bridges grey; rotor model red'),
              (flow_view, 'RETIRED WIND PILOT | west to east, 5.4 m/s\nOne vertical slice; arrows show airflow, not a route')]
    for view, title in labels:
        text = Text()
        text.Text = title
        display = Show(text, view)
        display.Color = [.10, .13, .17]
        display.FontSize = 12
    geometry.CameraPosition = [-1100, -1600, 1700]
    geometry.CameraFocalPoint = [-650, -45, 80]
    geometry.CameraParallelScale = 700
    flow_view.CameraPosition = [-650, -1800, 180]
    flow_view.CameraFocalPoint = [-650, -45, 180]
    flow_view.CameraParallelScale = 650
    RenderAllViews()
    SaveScreenshot(str(output / 'review.png'), layout, ImageResolution=[1800, 950])
    state = output / 'review.pvsm'
    SaveState(str(state))
    state.write_text(state.read_text().replace(str(output) + '/', './'))
    print('Saved portable review:', state)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    output = Path(__file__).resolve().parent
    config = export_data(args.source.resolve(), output)
    make_view(output, config)
