#!/usr/bin/env python3
"""Prepare AIJ E after-construction, at the measured 1:250 model scale."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import urllib.request

ROOT = Path(__file__).resolve().parent
IMAGE = 'sha256:7df8c49ff2f966deb0d8afba54a5deb19befc93ea4346e974f0f13869572c289'

def header(obj, cls='dictionary'):
    return f'FoamFile\n{{ version 2.0; format ascii; class {cls}; object {obj}; }}\n'

def write(path, body, cls='dictionary'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header(path.name, cls) + body + '\n')

def read_csv(name):
    with (ROOT/'data/raw'/name).open(encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))

def build():
    parser=argparse.ArgumentParser()
    parser.add_argument('--case', default='baseline')
    parser.add_argument('--empty', action='store_true', help='Prepare no-obstacle inflow-transport verification')
    parser.add_argument('--download', action='store_true', help='Fetch missing pinned inputs; never overwrite differing bytes')
    parser.add_argument('--epsilon-mode', choices=['mixing-length','shear-balance'],default='mixing-length',
                        help='shear-balance is an independent inlet-only sensitivity variant, not measured epsilon')
    parser.add_argument('--background-m',type=float,choices=[20.,30.,40.],help='Declared bounded numerical mesh variant')
    parser.add_argument('--consistent-simple',action='store_true',help='Final bounded SIMPLEC solver variant; numerical POC only')
    parser.add_argument('--solver-repair',action='store_true',help='Tighter inner pressure solve and two nonorthogonal corrections; numerical POC only')
    parser.add_argument('--mesh-repair',action='store_true',help='Uniform surface level 3 and conservative snapping; numerical POC only')
    parser.add_argument('--scenario',choices=['experimental','numerical-abl'],default='experimental')
    parser.add_argument('--closure',choices=['kEpsilon','realizableKE'])
    parser.add_argument('--resolution',choices=['cheap','reference'],default='reference')
    parser.add_argument('--direction',choices=['N','NE','NW','E','S','W'],default='N')
    args=parser.parse_args()
    if args.scenario=='numerical-abl':
        from poc import build_poc
        build_poc(args);return
    if args.closure is not None or args.direction!='N':
        raise SystemExit('Closure/direction extensions apply only to the explicit numerical-ABL scenario.')
    case=ROOT/'runs'/args.case
    if (case/'constant/polyMesh').exists():
        raise SystemExit('Refusing to overwrite a meshed case; choose a fresh case name.')
    lock=json.loads((ROOT/'inputs.json').read_text())
    for entry in lock['files']:
        target=ROOT/'data/raw'/entry['file']
        if not target.exists() and args.download:
            body=urllib.request.urlopen(entry['url'],timeout=60).read()
            assert hashlib.sha256(body).hexdigest()==entry['sha256']
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(body)
        assert hashlib.sha256(target.read_bytes()).hexdigest()==entry['sha256']
    af=[(float(r['z(m)'])/250,float(r['U(m/s)']),float(r['k(m2/s2)'])) for r in read_csv('AF_caseE.csv')]
    records=read_csv('RS_caseE.csv')
    ac=[r for r in records if r['case']=='ac']
    points=[r for r in ac if r['Wind_direction']=='N']
    assert len(records)==2560 and len(ac)==1280 and len(points)==80
    assert len({r['No.'] for r in points})==80
    probes=[tuple(float(r[k])/250 for k in ('x(m)','y(m)','z(m)')) for r in points]
    zref=15.9/250
    def interp(z, col):
        for left,right in zip(af,af[1:]):
            if left[0]<=z<=right[0]:
                f=(z-left[0])/(right[0]-left[0]);return left[col]*(1-f)+right[col]*f
        raise ValueError(z)
    uref=interp(zref,1)
    z0=0.0001
    epsilons=[]
    for i,(z,u,k) in enumerate(af):
        left=af[max(0,i-1)];right=af[min(len(af)-1,i+1)]
        derivative=(right[1]-left[1])/(right[0]-left[0])
        epsilons.append(.09**.5*k*abs(derivative))
    def epsilon(z,k):
        if args.epsilon_mode=='mixing-length':return .09**.75*k**1.5/(.41*(z+z0))
        if z<=af[0][0]:return epsilons[0]
        if z>=af[-1][0]:return epsilons[-1]
        for i,(left,right) in enumerate(zip(af,af[1:])):
            if left[0]<=z<=right[0]:
                fraction=(z-left[0])/(right[0]-left[0]);return epsilons[i]*(1-fraction)+epsilons[i+1]*fraction
        raise ValueError(z)
    params={'dataset_record':'15429018','case':'ac','direction':'N','wind_from_deg':0,
            'empty_inflow_verification':args.empty,
            'coordinate_frame':'x east, y north, z up; wind N -> (0,-U,0)',
            'coordinate_scale_csv_to_model':1/250,'stl_already_model_scale':True,
            'Uref_height_model_m':zref,'Uref_m_s':uref,'measurement_height_model_m':0.008,
            'held_out_directions':['E','S','W'],'fit_performed':False,
            'mesh_background_m':0.08,'surface_levels':[2,3],
            'ground_roughness_model_m':z0,'roughness_status':'assumed, not measured or fitted',
            'epsilon_model':'Cmu^0.75 k^1.5 / [kappa*(z+z0)]; Cmu=.09,kappa=.41; assumed mixing length',
            'epsilon_mode':args.epsilon_mode,
            'inlet_U_k':'piecewise linear measured profile, clamp above last point; zero U at floor',
            'boundary_assumption':'open-domain sides/top slip; exact facility confinement not supplied',
            'nu_m2_s':1.5e-5,'closure':'standard kEpsilon','iterations_cap':600,
            'Re_source_discrepancy':'Source README quotes Re approximately17000; supplied inletU at H=.24m with nu1.5e-5 gives approximately8e4. Preserve suppliedU and statedairnu; unresolved definition/metadata discrepancy.',
            'measurement_operator_limitation':'magnitude of RANS mean vector compared to thermistor mean scalar speed; these are not exactly identical statistics',
            'image':IMAGE,'source_hashes':{e['file']:e['sha256'] for e in lock['files']}}
    if args.epsilon_mode=='shear-balance':
        params['epsilon_model']='sqrt(Cmu)*k*abs(dU/dz), Cmu=.09; local P=epsilon equilibrium assumption using sourceAF finite differences; NOT measured/fitted to city response'
    (case/'constant/triSurface').mkdir(parents=True,exist_ok=True)
    shutil.copyfile(ROOT/'data/raw/BD_caseE.stl',case/'constant/triSurface/buildings.stl')
    profile=[(0.0,0.0,af[0][2])]+af+[(1.44,af[-1][1],af[-1][2])]
    def fixed_profile(index, vector=False):
        rows=[]
        for z,u,k in profile:
            value=(f'(0 {-u:.9g} 0)' if vector else
                   f'{k:.9g}' if index==2 else f'{epsilon(z,k):.9g}')
            rows.append(f'({z:.9g} {value})')
        return 'type fixedProfile; direction (0 0 1); origin 0; profile table\n(\n'+'\n'.join(rows)+'\n);'
    walls={'U':'type noSlip;', 'p':'type zeroGradient;',
           'k':'type kqRWallFunction; value uniform 0.5;',
           'epsilon':'type epsilonWallFunction; value uniform 1;',
           'nut':'type nutkWallFunction; value uniform 0;'}
    values={'U':f'(0 {-uref} 0)','p':'0','k':str(interp(zref,2)),
            'epsilon':str(epsilon(zref,interp(zref,2))), 'nut':'0'}
    dimensions={'U':'[0 1 -1 0 0 0 0]','p':'[0 2 -2 0 0 0 0]',
                'k':'[0 2 -2 0 0 0 0]','epsilon':'[0 2 -3 0 0 0 0]','nut':'[0 2 -1 0 0 0 0]'}
    for name in values:
        inlet=(fixed_profile(1,True) if name=='U' else fixed_profile(2) if name=='k' else
               fixed_profile(3) if name=='epsilon' else 'type zeroGradient;')
        if name=='U':outlet='type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);'
        elif name=='p':outlet='type fixedValue; value uniform 0;'
        elif name in ('k','epsilon'):outlet=f'type inletOutlet; inletValue uniform {values[name]}; value uniform {values[name]};'
        else:outlet='type calculated; value uniform 0;'
        if name=='nut':
            inlet='type calculated; value uniform 0;'
            ground=f'type atmNutkWallFunction; z0 uniform {z0}; value uniform 0;'
            side='type calculated; value uniform 0;'
        else:ground=walls[name];side='type slip;' if name=='U' else 'type zeroGradient;'
        body=f'dimensions {dimensions[name]};\ninternalField uniform {values[name]};\nboundaryField\n{{\n'
        body+=f'inlet {{ {inlet} }}\noutlet {{ {outlet} }}\nground {{ {ground} }}\n'
        body+=f'"(sides|top)" {{ {side} }}\nbuildings {{ {walls[name]} }}\n}}'
        write(case/'0'/name,body,'volVectorField' if name=='U' else 'volScalarField')
    write(case/'constant/transportProperties','transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 1.5e-5;')
    write(case/'constant/turbulenceProperties','simulationType RAS;\nRAS { RASModel kEpsilon; turbulence on; printCoeffs on; }')
    write(case/'system/blockMeshDict','''convertToMeters 1;
vertices ((-2 -4.4 0) (2 -4.4 0) (2 2 0) (-2 2 0)
          (-2 -4.4 1.44) (2 -4.4 1.44) (2 2 1.44) (-2 2 1.44));
blocks (hex (0 1 2 3 4 5 6 7) (50 80 18) simpleGrading (1 1 1));
edges ();
boundary (
outlet {type patch; faces ((0 1 5 4));}
inlet {type patch; faces ((3 7 6 2));}
sides {type patch; faces ((0 4 7 3) (1 2 6 5));}
ground {type wall; faces ((0 3 2 1));}
top {type patch; faces ((4 5 6 7));}
);
mergePatchPairs ();''')
    write(case/'system/snappyHexMeshDict','''castellatedMesh true; snap true; addLayers false;
geometry {
 buildings.stl { type triSurfaceMesh; name buildings; }
 pedestrian { type searchableBox; min (-0.6 -0.45 0); max (0.6 0.5 0.025); }
 district { type searchableBox; min (-0.85 -1.3 0); max (0.85 0.85 0.3); }
 badCellRefinement { type searchableBox; min (0.505 -0.58 0); max (0.545 -0.54 0.03); }
}
castellatedMeshControls {
 maxLocalCells 1000000; maxGlobalCells 1000000; minRefinementCells 5;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
 features ();
 refinementSurfaces { buildings {level (2 3); patchInfo {type wall;}} }
 resolveFeatureAngle 30;
 refinementRegions {pedestrian {mode inside; levels ((1E15 3));} district {mode inside; levels ((1E15 1));} badCellRefinement {mode inside; levels ((1E15 4));}}
 locationInMesh (1.5 1.5 0.4); allowFreeStandingZoneFaces true;
}
snapControls {nSmoothPatch 5; tolerance 1.5; nSolveIter 100; nRelaxIter 10;
 nFeatureSnapIter 10; implicitFeatureSnap true; explicitFeatureSnap false; multiRegionFeatureSnap false;}
addLayersControls {relativeSizes true; layers {}; expansionRatio 1; finalLayerThickness 0.3;
 minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1;
 nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
 minMedianAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 10;}
meshQualityControls {
 maxNonOrtho 65; maxBoundarySkewness 4; maxInternalSkewness 2; maxConcave 10;
 minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02;
 minDeterminant 0.01; minFaceWeight 0.02; minVolRatio 0.01;
 minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
}
debug 0; mergeTolerance 1e-6;''')
    if args.empty:
        path=case/'system/snappyHexMeshDict';text=path.read_text()
        text=text.replace('snap true;', 'snap false;')
        text=text.replace(' buildings.stl { type triSurfaceMesh; name buildings; }\n','')
        text=text.replace('pedestrian { type searchableBox; min (-0.6 -0.45 0); max (0.6 0.5 0.025); }',
                          'pedestrian { type searchableBox; min (-0.08 -0.8 0); max (0.08 2 0.08); }')
        text='\n'.join(line for line in text.splitlines() if not line.strip().startswith(('district {','badCellRefinement {')))
        text=text.replace('refinementSurfaces { buildings {level (2 3); patchInfo {type wall;}} }','refinementSurfaces {}')
        text='\n'.join(' refinementRegions {pedestrian {mode inside; levels ((1E15 3));}}' if line.strip().startswith('refinementRegions') else line for line in text.splitlines())
        path.write_text(text+'\n')
    quality=(case/'system/snappyHexMeshDict').read_text().split('meshQualityControls {',1)[1].split('\n}',1)[0]
    write(case/'system/meshQualityDict',quality)
    write(case/'system/fvSchemes','''ddtSchemes {default steadyState;}
gradSchemes {default cellLimited Gauss linear 1;}
divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear limited 0.5;}
interpolationSchemes {default linear;}
snGradSchemes {default limited 0.5;}
wallDist {method meshWave;}''')
    write(case/'system/fvSolution','''solvers {
 p {solver GAMG; tolerance 1e-7; relTol 0.05; smoother GaussSeidel;}
 "(U|k|epsilon)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.05;}
}
SIMPLE {nNonOrthogonalCorrectors 1; consistent no;
 residualControl {p 1e-5; U 1e-5; k 1e-5; epsilon 1e-5;}}
relaxationFactors {fields {p 0.3;} equations {U 0.6; k 0.6; epsilon 0.6;}}''')
    write(case/'system/decomposeParDict','numberOfSubdomains 3; method scotch;')
    point_text='\n'.join('('+ ' '.join(f'{v:.9g}' for v in p)+')' for p in probes)
    approach_points='\n'.join(f'(0 1.2 {z:.9g})' for z,_,_ in af)
    functions=f'''functions {{
 sample {{type sets; libs (sampling); writeControl timeStep; writeInterval 50;
  interpolationScheme cellPoint; setFormat raw; fields (U k);
  sets {{measurements {{type cloud; axis xyz; points ({point_text});}}
         approach {{type cloud; axis xyz; points ({approach_points});}}}}
 }}
 flux {{type surfaceFieldValue; libs (fieldFunctionObjects); writeFields false; writeControl timeStep;
  writeInterval 10; regionType patch; name inlet; operation sum; fields (phi);}}
 outflux {{type surfaceFieldValue; libs (fieldFunctionObjects); writeFields false; writeControl timeStep;
  writeInterval 10; regionType patch; name outlet; operation sum; fields (phi);}}
}}'''
    write(case/'system/controlDict',f'''application simpleFoam; startFrom startTime; startTime 0;
stopAt endTime; endTime 600; deltaT 1; writeControl timeStep; writeInterval 100;
purgeWrite 2; writeFormat ascii; writePrecision 10; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable false;
{functions}''')
    (case/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    (case/'measurements.json').write_text(json.dumps([{'point_id':int(r['No.']),'xyz_model_m':p,'reference_ratio':float(r['Velocity_Ratio'])} for r,p in zip(points,probes)],indent=2)+'\n')
    print(json.dumps({'case':str(case),'geometry_model_bounds_m':[[-.792,.790],[-.788,.793],[0,.24]],'N_measurements':len(points),'Uref_m_s':uref,'data_gates':'pass_for_scoped_scalar_speed_baseline'},indent=2))

if __name__=='__main__':build()
