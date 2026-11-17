"""Build the retained SWiFT neutral, steady actuator-disk pilot (stdlib only)."""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def foam(kind, obj, body):
    return (f'FoamFile {{ version 2.0; format ascii; class {kind}; object {obj}; }}\n'
            + body.strip() + '\n')


def rotor_coefficients(raw, uref):
    table = (raw / 'upstream/wt_data/v3/Simulated_data/RotorPerformances_OF2.csv').read_text().splitlines()
    rows = [[float(v) for v in line.split()] for line in table[1:] if line.strip()]
    if not rows[0][0] <= uref <= rows[-1][0]:
        raise ValueError(f'Hub speed {uref} is outside the locked v3 performance table.')
    low = max(r for r in rows if r[0] <= uref)
    high = min(r for r in rows if r[0] >= uref)
    fraction = 0. if high[0] == low[0] else (uref-low[0])/(high[0]-low[0])
    return tuple(low[i] + fraction*(high[i]-low[i]) for i in [4, 5])


def main(name='baseline-v2', level=2, iterations=600, empty=False, base=13.5, scheme='upwind', closure='kEpsilon', source_level=None, poc_role=None, uhub=None, roughness=None, condition_role='development', enlarged_domain=False, strict_stopping=False):
    started = time.perf_counter()
    if Path(name).name != name or name in {'.', '..'}:
        raise ValueError('Run name must be one directory name.')
    protocol = None
    if poc_role is not None:
        protocol_path = ROOT.parents[1] / 'scene-selection/poc-protocol.json'
        protocol = json.loads(protocol_path.read_text())
        if closure != protocol['common_toy_closure']['model']:
            raise ValueError('POC role must use the frozen common toy closure.')
        if source_level != 4 or base != 13.5:
            raise ValueError('POC comparisons require the common13.5m base and fixed D/32 source box.')
        if not 1 <= iterations <= 2000:
            raise ValueError('POC solve limit is1..2000 iterations within the stage time cap.')
        if uhub is None or roughness is None:
            raise ValueError('POC cases require explicit hub speed and roughness.')
    if uhub is not None and (not math.isfinite(uhub) or uhub <= 0):
        raise ValueError('Hub speed must be finite and positive.')
    if roughness is not None and (not math.isfinite(roughness) or roughness <= 0):
        raise ValueError('Roughness must be finite and positive.')
    if source_level is not None and source_level < level:
        raise ValueError('Fixed source refinement must be at least the wake level.')
    raw = ROOT / 'data/raw'
    lock = json.loads((ROOT / 'reference/data-sources.json').read_text())
    for item in lock['files']:
        path = raw / item['path']
        assert path.stat().st_size == item['bytes'], path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], path
    case = ROOT / 'runs' / name
    if case.exists():
        raise FileExistsError(f'Retained run exists: {case}; choose another --name')
    D, hub, uref, ustar, cmu, kappa = 27., 32.1, 8.7 if uhub is None else uhub, .45, .09, .41
    z0 = hub / math.expm1(uref * kappa / ustar) if roughness is None else roughness
    if roughness is not None:
        ustar = kappa * uref / math.log1p(hub/z0)
    k = ustar**2 / math.sqrt(cmu)
    eps = ustar**3 / (kappa*(hub + z0))
    cp, ct = rotor_coefficients(raw, uref)
    for sub in ['0', 'constant', 'system']:
        (case / sub).mkdir(parents=True, exist_ok=True)
    area = math.pi * (D/2)**2
    params = dict(name=name, source='SWiFT neutral ensemble, v3', rotor_enabled=not empty, rotor_D_m=D,
                  hub_m=hub, Uhub_m_s=uref, ustar_m_s=ustar, z0_m=z0,
                  k_equilibrium_m2_s2=k, measured_three_component_k_m2_s2=.822,
                  published_horizontal_TI_percent=10.7, Cp=cp, Ct=ct,
                  area_m2=area, disk_thickness_m=D/4, refinement_level=level,
                  base_cell_m=base, smallest_wake_cell_m=base/2**level, max_iterations=iterations,
                  source_refinement_level=source_level if source_level is not None else level,
                  source_cell_m=base/2**(source_level if source_level is not None else level),
                  velocity_convection_scheme=scheme,
                  density_kg_m3=1.02, expected_T_at_Uhub_N=.5*1.02*area*ct*uref**2,
                  closure=closure, force_variant='variableScaling',
                  monitor_points_m=[[-65.,0.,hub]],
                  model_limitations=['Tower/nacelle drag omitted in this first pilot; retained display geometry is not a force model.',
                    'Neutral equilibrium k differs from measured component variance and published horizontal TI.',
                    'Uniform initial field; equilibrium inlet is transported during convergence.',
                    'Flat ground, zero tilt/yaw and no swirl; averaged near-neutral comparison.',
                    'No mesh/domain independence established by this pilot.'],
                  data_split=dict(training_dates=['2017-01-23'], validation_dates=['2017-06-22'],
                                  final_test_dates=['2017-07-03'],
                                  initial_unfitted_ensemble_uses_all_dates=True,
                                  restriction='No calibration/generalization claim; reserve date groups for any later fitting.'))
    params['domain_bounds_m'] = [[-270,-270,0],[540,270,270]] if enlarged_domain else [[-135,-135,0],[405,135,162]]
    params['data_split'].update(current_population_role='DEVELOPMENT', restriction='All six measured periods have entered method selection; historical date groups are exploratory, not a confirmatory holdout.')
    if protocol is not None:
        params['poc_study'] = dict(role=poc_role, protocol_revision=protocol['revision'],
            protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
            condition_role=condition_role,
            max_cells=protocol['resources']['cheap_target_max_cells' if poc_role == 'cheap' else 'reference_max_cells'],
            max_mesh_seconds=60*protocol['resources']['max_mesh_minutes'],
            max_solve_seconds=60*protocol['resources']['max_solve_minutes'],
            residual_control=protocol['convergence'][('reference' if strict_stopping else poc_role)+'_residual_control'],
            stopping_recipe='reference-strict' if strict_stopping else poc_role,
            observation_comparison_eligible=abs(uref-8.7)<1e-10 and abs(z0-hub/math.expm1(8.7*kappa/.45))<1e-5,
            reference_kind='Checked numerical RANS reference; not physical truth')
        support = dict(base_cell_m=base, grid_origin_m=[-135,-135,0],
            box_min_m=[-6.75,-20.25,6.75], box_max_m=[6.75,20.25,60.75],
            refinement_level=source_level, cell_m=base/2**source_level,
            cylinder_p1_m=[-3.375,0,hub], cylinder_p2_m=[3.375,0,hub], radius_m=D/2)
        params['force_support'] = dict(**support, configuration_sha256=hashlib.sha256(json.dumps(support,sort_keys=True).encode()).hexdigest(),
            proof_scope='Same axis-aligned base lattice and source refinement box; selected source count and uniform refinement level must also be checked after meshing.')
    (case / 'parameters.json').write_text(json.dumps(params, indent=2)+'\n')
    def put(rel, kind, obj, text):
        (case / rel).write_text(foam(kind, obj, text))
    put('system/blockMeshDict', 'dictionary', 'blockMeshDict', '''
convertToMeters 1;
vertices ((-135 -135 0) (405 -135 0) (405 135 0) (-135 135 0)
          (-135 -135 162) (405 -135 162) (405 135 162) (-135 135 162));
blocks (hex (0 1 2 3 4 5 6 7) (NX NY NZ) simpleGrading (1 1 1));
edges ();
boundary (
 inlet {type patch; faces ((0 4 7 3));}
 outlet {type patch; faces ((1 2 6 5));}
 sides {type patch; faces ((0 1 5 4) (3 7 6 2));}
 ground {type wall; faces ((0 3 2 1));}
 top {type patch; faces ((4 5 6 7));}
);
mergePatchPairs ();
'''.replace('NX',str(round(540/base))).replace('NY',str(round(270/base))).replace('NZ',str(round(162/base))))
    if enlarged_domain:
        # Integer base-cell offsets retain the exact local mesh lattice.
        mesh_path = case/'system/blockMeshDict'
        text = mesh_path.read_text()
        text = text.replace('(-135 -135 0) (405 -135 0) (405 135 0) (-135 135 0)',
                            '(-270 -270 0) (540 -270 0) (540 270 0) (-270 270 0)')
        text = text.replace('(-135 -135 162) (405 -135 162) (405 135 162) (-135 135 162)',
                            '(-270 -270 270) (540 -270 270) (540 270 270) (-270 270 270)')
        text = text.replace(f'({round(540/base)} {round(270/base)} {round(162/base)}) simpleGrading',
                            f'({round(810/base)} {round(540/base)} {round(270/base)}) simpleGrading')
        mesh_path.write_text(text)
    support_geometry = ('sourceSupport {type searchableBox; min (-6.75 -20.25 6.75); max (6.75 20.25 60.75);}'
                        if source_level is not None else '')
    support_refinement = (f'sourceSupport {{mode inside; levels ((1e15 {source_level}));}}'
                          if source_level is not None else '')
    if protocol is not None:
        # A common near-ground band controls the already identified cellPoint/wall artifact.
        support_geometry += ' groundSupport {type searchableBox; min (-40.5 -54 0); max (189 54 6.75);}'
        support_refinement += ' groundSupport {mode inside; levels ((1e15 3));}'
        params['poc_study']['fixed_ground_band'] = dict(min_m=[-40.5,-54,0],max_m=[189,54,6.75],refinement_level=3,
            reason='Prespecified common near-ground resolution, independent of surrounding wake level; not fitted to observations.')
        (case/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    put('system/snappyHexMeshDict', 'dictionary', 'snappyHexMeshDict', f'''
castellatedMesh true; snap false; addLayers false;
geometry {{ wake {{type searchableBox; min (-27 -40.5 0); max (175.5 40.5 81);}} {support_geometry} }}
castellatedMeshControls
{{
 maxLocalCells 1000000; maxGlobalCells 1000000; minRefinementCells 0;
 nCellsBetweenLevels 3; maxLoadUnbalance 0.1; features ();
 refinementSurfaces {{}} resolveFeatureAngle 30;
 refinementRegions {{wake {{mode inside; levels ((1e15 {level}));}} {support_refinement}}}
 locationInMesh (-100 0 100); allowFreeStandingZoneFaces true;
}}
snapControls {{nSmoothPatch 3; tolerance 2; nSolveIter 30; nRelaxIter 5;}}
addLayersControls {{relativeSizes true; layers {{}} expansionRatio 1.2;
 finalLayerThickness .3; minThickness .1; nGrow 0; featureAngle 60;
 nRelaxIter 5; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
 maxFaceThicknessRatio .5; maxThicknessToMedialRatio .3; minMedianAxisAngle 90;
 nBufferCellsNoExtrude 0; nLayerIter 50;}}
meshQualityControls {{maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4;
 maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist .02;
 minDeterminant .001; minFaceWeight .02; minVolRatio .01; minTriangleTwist -1;
 nSmoothScale 4; errorReduction .75;}}
writeFlags (); mergeTolerance 1e-6;
''')
    put('system/topoSetDict', 'dictionary', 'topoSetDict', '''
actions (
 {name rotor; type cellSet; action new; source cylinderToCell;
  p1 (-3.375 0 32.1); p2 (3.375 0 32.1); radius 13.5;}
);
''')
    put('constant/fvOptions', 'dictionary', 'fvOptions', f'''
rotor
{{
 type actuationDiskSource; variant variableScaling; selectionMode cellSet;
 cellSet rotor; diskArea {area:.12g}; diskDir (1 0 0); sink true;
 Cp {cp:.12g}; Ct {ct:.12g};
 monitorMethod points; monitorCoeffs {{points ((-65 0 {hub}));}}
 writeToFile true; writeFileStart 0; writeFileEnd 1000000;
}}
''')
    if empty:
        (case / 'constant/fvOptions').unlink()
    put('constant/transportProperties', 'dictionary', 'transportProperties',
        'transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 1.5e-5;')
    put('constant/turbulenceProperties', 'dictionary', 'turbulenceProperties',
        f'simulationType RAS;\nRAS {{RASModel {closure}; turbulence on; printCoeffs on;}}')
    if closure != 'kEpsilon' and protocol is None:
        params['sensitivity'] = dict(
            kind='unfitted stock closure sensitivity', changed_property='RASModel',
            baseline_model='kEpsilon', tested_model=closure,
            restriction='Stock installed model, not an atmospheric-adapted paper model. Same imposed inlet/ground parameters; matched empty-ABL transport must be assessed before wake interpretation.')
        params['data_split']['restriction'] = ('All six periods are now a DEVELOPMENT reference used to choose methods. '
            'Previously nominated date groups are exploratory only; no confirmatory holdout remains for this closure selection.')
        (case/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    abl = f'flowDir (1 0 0); zDir (0 0 1); Uref {uref}; Zref {hub}; z0 uniform {z0:.12g}; zGround uniform 0; d uniform 0;'
    patches = {
        'U': ('volVectorField', '[0 1 -1 0 0 0 0]', f'({uref} 0 0)',
              f'inlet {{type atmBoundaryLayerInletVelocity; {abl} value uniform ({uref} 0 0);}}\n'
              f'outlet {{type inletOutlet; inletValue uniform ({uref} 0 0); value uniform ({uref} 0 0);}}\n'
              'ground {type noSlip;}\nsides {type slip;} top {type slip;}'),
        'p': ('volScalarField', '[0 2 -2 0 0 0 0]', '0',
              'inlet {type zeroGradient;} outlet {type fixedValue; value uniform 0;}\n'
              'ground {type zeroGradient;} sides {type zeroGradient;} top {type zeroGradient;}'),
        'k': ('volScalarField', '[0 2 -2 0 0 0 0]', str(k),
              f'inlet {{type atmBoundaryLayerInletK; {abl} value uniform {k};}}\n'
              f'outlet {{type inletOutlet; inletValue uniform {k}; value uniform {k};}}\n'
              f'ground {{type kqRWallFunction; value uniform {k};}} sides {{type zeroGradient;}} top {{type zeroGradient;}}'),
        'epsilon': ('volScalarField', '[0 2 -3 0 0 0 0]', str(eps),
              f'inlet {{type atmBoundaryLayerInletEpsilon; {abl} value uniform {eps};}}\n'
              f'outlet {{type inletOutlet; inletValue uniform {eps}; value uniform {eps};}}\n'
              f'ground {{type epsilonWallFunction; value uniform {eps};}} sides {{type zeroGradient;}} top {{type zeroGradient;}}'),
        'nut': ('volScalarField', '[0 2 -1 0 0 0 0]', '0',
              'inlet {type calculated; value uniform 0;} outlet {type calculated; value uniform 0;}\n'
              f'ground {{type atmNutkWallFunction; z0 uniform {z0:.12g}; boundNut false; value uniform 0;}}\n'
              'sides {type calculated; value uniform 0;} top {type calculated; value uniform 0;}'),
    }
    for field, (kind, dims, value, boundary) in patches.items():
        put(f'0/{field}', kind, field, f'dimensions {dims};\ninternalField uniform {value};\nboundaryField {{\n{boundary}\n}}')
    put('system/fvSchemes', 'dictionary', 'fvSchemes', '''
ddtSchemes {default steadyState;}
gradSchemes {default Gauss linear;}
divSchemes {default none; div(phi,U) bounded Gauss upwind;
 div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;}
snGradSchemes {default corrected;}
wallDist {method meshWave;}
''')
    if scheme == 'linearUpwind':
        path = case / 'system/fvSchemes'
        path.write_text(path.read_text().replace('div(phi,U) bounded Gauss upwind;',
                                                'div(phi,U) bounded Gauss linearUpwind grad(U);'))
    put('system/fvSolution', 'dictionary', 'fvSolution', '''
solvers {
 p {solver GAMG; tolerance 1e-7; relTol .05; smoother GaussSeidel;}
 "(U|k|epsilon)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol .05;}
}
SIMPLE {nNonOrthogonalCorrectors 1; residualControl {p 1e-5; U 1e-5; k 1e-5; epsilon 1e-5;}}
relaxationFactors {fields {p .3;} equations {U .7; k .7; epsilon .7;}}
''')
    if protocol is not None:
        path = case / 'system/fvSolution'
        old = 'residualControl {p 1e-5; U 1e-5; k 1e-5; epsilon 1e-5;}'
        limits = params['poc_study']['residual_control']
        new = 'residualControl {'+' '.join(f'{key} {value:.12g};' for key,value in limits.items())+'}'
        assert old in path.read_text()
        path.write_text(path.read_text().replace(old,new))
    put('system/decomposeParDict', 'dictionary', 'decomposeParDict',
        'numberOfSubdomains 3;\nmethod scotch;')
    cloud = []
    for i in range(39):
        for j in range(49):
            for iz in range(33):
                cloud.append(f'({20.25+i*3.375:g} {-81+j*3.375:g} {1.6875+iz*3.375:g})')
    put('system/sampleGrid', 'dictionary', 'sampleGrid', '''
type sets; libs (sampling); interpolationScheme cellPoint; setFormat raw;
fields (U k); sets {grid {type cloud; axis xyz; points (
'''+'\n'.join(cloud)+');}}')
    write_every = 2000 if protocol is not None else 200
    put('system/controlDict', 'dictionary', 'controlDict', f'''
application simpleFoam; startFrom startTime; startTime 0; stopAt endTime;
endTime {iterations}; deltaT 1; writeControl timeStep; writeInterval {write_every};
purgeWrite 0; writeFormat ascii; writePrecision 10; timeFormat general; timePrecision 8;
runTimeModifiable true;
functions
{{
 probes {{type probes; libs (sampling); fields (U k); writeControl timeStep; writeInterval 1;
 probeLocations ((-65 0 32.1) (0 0 32.1) (54 0 32.1) (81 0 32.1) (108 0 32.1) (135 0 32.1));}}
 inletFlux {{type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name inlet;
 operation sum; fields (phi); writeControl timeStep; writeInterval 1; writeFields false;}}
 outletFlux {{type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name outlet;
 operation sum; fields (phi); writeControl timeStep; writeInterval 1; writeFields false;}}
}}
''')
    if protocol is not None:
        regions = {
            'wakeROI': [(27+3.375*i,-40.5+3.375*j,1.6875+3.375*iz) for i in range(33) for j in range(25) for iz in range(24)],
            'background': [(-128.25+13.5*i,-128.25+13.5*j,6.75+13.5*iz) for i in range(40) for j in range(20) for iz in range(12)],
        }
        body = 'type sets; libs (sampling); interpolationScheme cellPoint; setFormat raw; fields (U k); sets {\n'
        for key,points in regions.items():
            body += key+' {type cloud; axis xyz; points (\n'
            body += '\n'.join(f'({x:g} {y:g} {z:g})' for x,y,z in points)+');}\n'
        body += '}\n'
        put('system/sampleStudy','dictionary','sampleStudy',body)
        params['poc_study']['query_regions'] = {key:dict(points=len(points),definition='Fixed common native cellPoint query coordinates; no extrapolation or moved points') for key,points in regions.items()}
        (case/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    (case/'cost.prepare.json').write_text(json.dumps(dict(stage='prepare',elapsed_wall_seconds=time.perf_counter()-started,
        timing_scope='Pinned source verification and complete case/query dictionary preparation; source acquisition excluded.'),indent=2)+'\n')
    (case / 'swift.foam').touch()
    print(json.dumps(dict(case=str(case), source_files_verified=len(lock['files']), parameters=params), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', default='baseline-v2')
    p.add_argument('--level', type=int, default=2, choices=[1,2,3])
    p.add_argument('--iterations', type=int, default=600)
    p.add_argument('--empty', action='store_true', help='Undisturbed-inflow transport check, no rotor force')
    p.add_argument('--base', type=float, default=13.5, choices=[13.5,9.])
    p.add_argument('--scheme', default='upwind', choices=['upwind','linearUpwind'])
    p.add_argument('--closure', default='kEpsilon', choices=['kEpsilon','realizableKE'])
    p.add_argument('--source-level', type=int, choices=[2,3,4], help='Fixed rotor/source refinement independent of surrounding wake level')
    p.add_argument('--poc-role', choices=['cheap','reference'])
    p.add_argument('--strict-stopping', action='store_true')
    p.add_argument('--enlarged-domain', action='store_true')
    p.add_argument('--uhub', type=float)
    p.add_argument('--roughness', type=float)
    p.add_argument('--condition-role', default='development', choices=['development','numerical-holdout'])
    args = p.parse_args()
    main(args.name, args.level, args.iterations, args.empty, args.base, args.scheme, args.closure, args.source_level,
         args.poc_role, args.uhub, args.roughness, args.condition_role, args.enlarged_domain, args.strict_stopping)
