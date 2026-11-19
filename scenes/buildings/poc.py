#!/usr/bin/env python3
"""Explicit numerical-only full-size Niigata ABL case and matched field comparison."""
import hashlib
import json
import math
import platform
import re
import socket
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT.parent/'poc-protocol.json'
IMAGE = 'sha256:7df8c49ff2f966deb0d8afba54a5deb19befc93ea4346e974f0f13869572c289'
DIRECTIONS = {'N':0., 'NbE':11.25, 'NNE':22.5, 'NE':45., 'NW':315., 'E':90., 'S':180., 'W':270.}

def write(path, body, cls='dictionary'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'FoamFile {{version 2.0; format ascii; class {cls}; object {path.name};}}\n'+body+'\n')

def transform(point, degrees, inverse=False):
    c=math.cos(math.radians(degrees));s=math.sin(math.radians(degrees))
    x,y,z=point
    return (c*x+s*y,-s*x+c*y,z) if inverse else (c*x-s*y,s*x+c*y,z)

def box_string(a,b):
    return 'min ('+' '.join(map(str,a))+'); max ('+' '.join(map(str,b))+');'

def bbox_world_to_solver(a,b,degrees):
    ps=[transform((x,y,z),degrees) for x in (a[0],b[0]) for y in (a[1],b[1]) for z in (a[2],b[2])]
    return tuple(min(p[i] for p in ps) for i in range(3)),tuple(max(p[i] for p in ps) for i in range(3))

def world_refinement_box(low,high,degrees,oriented=False):
    if not (oriented and degrees):
        return 'type searchableBox; '+box_string(*bbox_world_to_solver(low,high,degrees))
    # OpenFOAM v2306 maps [0, span] from its local box frame using this corner and axes.
    origin=transform(low,degrees)
    span=tuple(b-a for a,b in zip(low,high))
    axis=transform((1,0,0),degrees)
    def vector(values):
        return '('+' '.join(map(str,values))+')'
    return (f'type searchableRotatedBox; origin {vector(origin)}; span {vector(span)}; '
            f'e1 {vector(axis)}; e3 (0 0 1);')


def build_poc(args, case_root=None):
    prepare_start=time.monotonic();prepare_cpu=time.process_time()
    protocol=json.loads(PROTOCOL.read_text())
    if protocol['revision']!=3:raise SystemExit('This case builder implements frozen protocol revision 3 only.')
    if Path(args.case).name != args.case or args.case in {'.', '..'}:
        raise ValueError('Case name must be one directory name.')
    case=(Path(case_root) if case_root is not None else ROOT/'runs')/args.case
    if case.exists():raise SystemExit('Numerical cases require a fresh name; existing evidence is never overwritten.')
    lock=json.loads((ROOT/'inputs.json').read_text())
    for entry in lock['files']:
        p=ROOT/'data/raw'/entry['file']
        if hashlib.sha256(p.read_bytes()).hexdigest()!=entry['sha256']:raise SystemExit('Input hash mismatch: '+str(p))
    closure=args.closure or protocol['common_toy_closure']['model']
    direction=args.direction
    theta=DIRECTIONS[direction]
    cheap=args.resolution=='cheap'
    mesh_repair=getattr(args,'mesh_repair',False)
    solver_repair=getattr(args,'solver_repair',False)
    surface_levels=[3,3] if mesh_repair else [2,3]
    snap_tolerance=.5 if mesh_repair else 1.5
    dx=getattr(args,'background_m',None) or (40. if cheap else 20.)
    uref=8.7;zref=32.1;z0=.01158;kappa=.41;cmu=.09
    ustar=uref*kappa/math.log((zref+z0)/z0)
    k=ustar**2/math.sqrt(cmu)
    eps=ustar**3/(kappa*(zref+z0))
    maxcells=protocol['resources']['cheap_target_max_cells' if cheap else 'reference_max_cells']
    residual=protocol['convergence']['cheap_residual_control' if cheap else 'reference_residual_control']
    params={
        'scenario':'numerical-abl','scope':'Numerical-only full source Niigata geometry; does not reproduce the wind-tunnel inflow or validate physical city flow.',
        'protocol_revision':3,'protocol_sha256':hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        'empty_inflow_verification':args.empty,'direction':direction,'wind_from_deg':theta,
        'coordinate_frame':'World x east/y north/z up; solver frame rotates geometry and queries so incoming wind is -y.',
        'raw_stl_to_world_scale':250.,'closure':closure,'closure_status':'Frozen shared-toy assumption; not evidence of superior physical city accuracy.',
        'resolution':args.resolution,'mesh_background_m':dx,'surface_levels':surface_levels,'snap_tolerance':snap_tolerance,
        'repair_recipe':{'mesh':mesh_repair,'solver':solver_repair,'consistent_simple':getattr(args,'consistent_simple',False)},
        'Uref_m_s':uref,'Uref_height_m':zref,'z0_m':z0,'ustar_m_s':ustar,'k_m2_s2':k,
        'max_cells':maxcells,'iterations_cap':2000,'mesh_timeout_s':600,'solve_timeout_s':600,
        'residual_thresholds':residual,'flux_threshold':1e-4,'probe_change_RMSE_over_Uref_threshold':.002,
        'image':IMAGE,'host':socket.gethostname(),'host_platform':platform.platform(),
        'source_hashes':{e['file']:e['sha256'] for e in lock['files']},
        'measured_RMSE_gate':'not applicable: new incoming conditions differ from original experiment',
        'boundary_assumption':'Open-domain sides/top slip; common coherent neutral ABL ground/inlet; domain sensitivity not yet established.',
        'comparison_sampling':'Exact fixed world grid; common valid fluid samples only; no moved, filled or extrapolated samples.',
    }
    case.mkdir(parents=True)
    source=ROOT/'data/raw/BD_caseE.stl';dest=case/'constant/triSurface/buildings.stl'
    dest.parent.mkdir(parents=True)
    vertex=re.compile(r'(\s*vertex\s+)([-+\deE.]+)\s+([-+\deE.]+)\s+([-+\deE.]+)')
    normal=re.compile(r'(\s*facet\s+normal\s+)([-+\deE.]+)\s+([-+\deE.]+)\s+([-+\deE.]+)')
    with source.open() as incoming,dest.open('w') as outgoing:
        for line in incoming:
            match=vertex.fullmatch(line.rstrip())
            norm=normal.fullmatch(line.rstrip())
            if match:
                p=transform(tuple(250*float(match[i]) for i in (2,3,4)),theta)
                outgoing.write(match[1]+' '.join(f'{x:.12g}' for x in p)+'\n')
            elif norm:
                n=transform(tuple(float(norm[i]) for i in (2,3,4)),theta)
                outgoing.write(norm[1]+' '.join(f'{x:.12g}' for x in n)+'\n')
            else:outgoing.write(line)
    params['scaled_stl_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest()
    scalar={'p':0.,'k':k,'epsilon':eps,'nut':0.}
    dimensions={'U':'[0 1 -1 0 0 0 0]','p':'[0 2 -2 0 0 0 0]',
                'k':'[0 2 -2 0 0 0 0]','epsilon':'[0 2 -3 0 0 0 0]','nut':'[0 2 -1 0 0 0 0]'}
    common=f'flowDir (0 -1 0); zDir (0 0 1); Uref {uref}; Zref {zref}; z0 uniform {z0}; d uniform 0;'
    for name in dimensions:
        value=f'(0 {-uref} 0)' if name=='U' else str(scalar[name])
        if name=='U':
            inlet=f'type atmBoundaryLayerInletVelocity; {common} value uniform {value};'
            outlet=f'type inletOutlet; inletValue uniform (0 0 0); value uniform {value};'
            ground='type noSlip;';building='type noSlip;';side='type slip;'
        elif name=='p':
            inlet='type zeroGradient;';outlet='type fixedValue; value uniform 0;'
            ground=building=side='type zeroGradient;'
        elif name in ('k','epsilon'):
            kind='atmBoundaryLayerInletK' if name=='k' else 'atmBoundaryLayerInletEpsilon'
            inlet=f'type {kind}; {common} value uniform {value};'
            outlet=f'type inletOutlet; inletValue uniform {value}; value uniform {value};'
            ground=building=(f'type kqRWallFunction; value uniform {value};' if name=='k' else f'type epsilonWallFunction; value uniform {value};')
            side='type zeroGradient;'
        else:
            inlet=outlet=side='type calculated; value uniform 0;'
            ground=f'type atmNutkWallFunction; z0 uniform {z0}; value uniform 0;'
            building='type nutkWallFunction; value uniform 0;'
        write(case/'0'/name,f'dimensions {dimensions[name]};\ninternalField uniform {value};\nboundaryField {{ inlet {{{inlet}}} outlet {{{outlet}}} ground {{{ground}}} buildings {{{building}}} "(sides|top)" {{{side}}} }}',
              'volVectorField' if name=='U' else 'volScalarField')
    write(case/'constant/transportProperties','transportModel Newtonian; nu [0 2 -1 0 0 0 0] 1.5e-5;')
    write(case/'constant/turbulenceProperties',f'simulationType RAS; RAS {{RASModel {closure}; turbulence on; printCoeffs on;}}')
    write(case/'system/blockMeshDict',f'''convertToMeters 1;
vertices ((-600 -1200 0) (600 -1200 0) (600 600 0) (-600 600 0)
          (-600 -1200 360) (600 -1200 360) (600 600 360) (-600 600 360));
blocks (hex (0 1 2 3 4 5 6 7) ({int(1200/dx)} {int(1800/dx)} {int(360/dx)}) simpleGrading (1 1 1));
edges ();
boundary (
outlet {{type patch; faces ((0 1 5 4));}} inlet {{type patch; faces ((3 7 6 2));}}
sides {{type patch; faces ((0 4 7 3) (1 2 6 5));}}
ground {{type wall; faces ((0 3 2 1));}} top {{type patch; faces ((4 5 6 7));}}
); mergePatchPairs ();''')
    oriented=getattr(args,'oriented_world_boxes',False)
    floor=world_refinement_box((-205,-205,0),(205,205,10),theta,oriented)
    district=world_refinement_box((-230,-400,0),(230,230,100),theta,oriented)
    repair=world_refinement_box((126.25,-145,0),(136.25,-135,7.5),theta,oriented)
    geometry='' if args.empty else ' buildings.stl {type triSurfaceMesh; name buildings;}\n'
    surfaces='' if args.empty else f'buildings {{level ({surface_levels[0]} {surface_levels[1]}); patchInfo {{type wall;}}}}'
    regions=' floorRegion {mode inside; levels ((1E15 3));} inletStrip {mode inside; levels ((1E15 3));}'
    if not args.empty:regions+=' district {mode inside; levels ((1E15 1));} repair {mode inside; levels ((1E15 4));}'
    write(case/'system/snappyHexMeshDict',f'''castellatedMesh true; snap {'false' if args.empty else 'true'}; addLayers false;
geometry {{
{geometry} floorRegion {{{floor}}}
 inletStrip {{type searchableBox; min (-10 200 0); max (10 600 80);}}
 district {{{district}}}
 repair {{{repair}}}
}}
castellatedMeshControls {{
 maxLocalCells {maxcells}; maxGlobalCells {maxcells}; minRefinementCells 5;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3; features ();
 refinementSurfaces {{{surfaces}}} resolveFeatureAngle 30;
 refinementRegions {{{regions}}}
 locationInMesh (500 500 200); allowFreeStandingZoneFaces true;
}}
snapControls {{nSmoothPatch 5; tolerance {snap_tolerance}; nSolveIter 100; nRelaxIter 10;
 nFeatureSnapIter 10; implicitFeatureSnap true; explicitFeatureSnap false; multiRegionFeatureSnap false;}}
addLayersControls {{relativeSizes true; layers {{}}; expansionRatio 1; finalLayerThickness 0.3;
 minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1;
 nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
 minMedianAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 10;}}
meshQualityControls {{
 maxNonOrtho 65; maxBoundarySkewness 4; maxInternalSkewness 2; maxConcave 10;
 minVol 1e-8; minTetQuality 1e-15; minArea -1; minTwist 0.02;
 minDeterminant 0.01; minFaceWeight 0.02; minVolRatio 0.01;
 minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
}}
debug 0; mergeTolerance 1e-6;''')
    quality=(case/'system/snappyHexMeshDict').read_text().split('meshQualityControls {',1)[1].split('\n}',1)[0]
    write(case/'system/meshQualityDict',quality)
    write(case/'system/fvSchemes', '''ddtSchemes {default steadyState;} gradSchemes {default cellLimited Gauss linear 1;}
divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear limited 0.5;}
interpolationSchemes {default linear;} snGradSchemes {default limited 0.5;}
wallDist {method meshWave;}''')
    write(case/'system/fvSolution',f'''solvers {{
 p {{solver GAMG; tolerance 1e-7; relTol 0.05; smoother GaussSeidel;}}
 "(U|k|epsilon)" {{solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.05;}}
}}
SIMPLE {{nNonOrthogonalCorrectors 1; consistent no; residualControl {{p {residual['p']}; U {residual['U']}; k {residual['k']}; epsilon {residual['epsilon']};}}}}
relaxationFactors {{fields {{p 0.3;}} equations {{U 0.6; k 0.6; epsilon 0.6;}}}}''')
    if solver_repair:
        path=case/'system/fvSolution'
        text=path.read_text().replace('tolerance 1e-7; relTol 0.05; smoother GaussSeidel;', 'tolerance 1e-8; relTol 0.01; smoother GaussSeidel;').replace('nNonOrthogonalCorrectors 1;', 'nNonOrthogonalCorrectors 2;')
        path.write_text(text)
    if getattr(args,'consistent_simple',False):
        path=case/'system/fvSolution';path.write_text(path.read_text().replace('consistent no;', 'consistent yes;'))
    write(case/'system/decomposeParDict','numberOfSubdomains 3; method scotch;')
    monitor_world=[(x,y,z) for x,y in ((0,150),(0,0),(0,-150)) for z in (10,40,80)]
    monitor_solver=[transform(p,theta) for p in monitor_world]
    monitor=' '.join('('+' '.join(f'{v:.12g}' for v in p)+')' for p in monitor_solver)
    functions=f'''functions {{
 monitor {{type sets; libs (sampling); writeControl timeStep; writeInterval 25;
 interpolationScheme cellPoint; setFormat raw; fields (U);
 sets {{ probes {{type cloud; axis xyz; points ({monitor});}} }} }}
 flux {{type surfaceFieldValue; libs (fieldFunctionObjects); writeFields false; writeControl timeStep;
 writeInterval 10; regionType patch; name inlet; operation sum; fields (phi);}}
 outflux {{type surfaceFieldValue; libs (fieldFunctionObjects); writeFields false; writeControl timeStep;
 writeInterval 10; regionType patch; name outlet; operation sum; fields (phi);}}
}}'''
    write(case/'system/controlDict',f'''application simpleFoam; startFrom startTime; startTime 0;
stopAt endTime; endTime 2000; deltaT 1; writeControl timeStep; writeInterval 100;
purgeWrite 2; writeFormat ascii; writePrecision 10; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable false;
{functions}''')
    queries=[]
    for z in (5,10,20,40,80):
        for y in range(-200,201,5):
            for x in range(-200,201,5):queries.append(('city',(float(x),float(y),float(z))))
    for z in (5,10,20,40,80):
        for y in range(-1000,501,40):
            for x in range(-500,501,40):queries.append(('domain',(float(x),float(y),float(z))))
    for z in (5,10,20,32.1,40,80,160):
        queries.append(('inflow',(0.,350.,float(z))))
    sets=[]
    query_records=[]
    for group in ('city','domain','inflow'):
        group_points=[p for name,p in queries if name==group]
        text=' '.join('('+' '.join(f'{v:.12g}' for v in transform(p,theta))+')' for p in group_points)
        sets.append(f'{group} {{type cloud; axis xyz; points ({text});}}')
        for index,p in enumerate(group_points):query_records.append({'region':group,'id':index,'world_xyz_m':p,'solver_xyz_m':transform(p,theta)})
    write(case/'system/roiSampling','functions {roi {type sets; libs (sampling); writeControl timeStep; writeInterval 1; interpolationScheme cellPoint; setFormat raw; fields (U k); sets {'+'\n'.join(sets)+'}}}')
    (case/'queries.json').write_text(json.dumps(query_records,separators=(',',':'))+'\n')
    params['query_sha256']=hashlib.sha256((case/'queries.json').read_bytes()).hexdigest()
    (case/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    (case/'cost.prepare.json').write_text(json.dumps({'stage':'prepare','elapsed_s':time.monotonic()-prepare_start,'cpu_s':time.process_time()-prepare_cpu,'peak_cgroup_memory_bytes_including_cache':None,'scope':'Host case preparation, geometry transformation and fixed-query generation; source acquisition excluded','host':params['host']},indent=2)+'\n')
    print(json.dumps({'case':str(case),'scenario':params['scenario'],'closure':closure,'resolution':args.resolution,'Uref':uref,'query_points':len(queries),'max_cells':maxcells}))

def raw_rows(path):
    return [list(map(float,line.split())) for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]

def numeric_times(folder):
    return sorted((p for p in folder.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)),key=lambda p:float(p.name))

def sampled_world(case):
    params=json.loads((case/'parameters.json').read_text())
    queries=json.loads((case/'queries.json').read_text())
    folders=numeric_times(case/'postProcessing/roi')
    last=folders[-1]
    result={}
    for region in ('city','domain','inflow'):
        matching=list(last.glob(region+'_*.xy'))
        path=next((p for p in matching if p.stem==region+'_k_U'),None)
        if path:
            values=[r[:3]+r[4:7] for r in raw_rows(path)]
        else:values=raw_rows(last/f'{region}_U.xy')
        lookup={tuple(round(v,7) for v in row[:3]):row[3:6] for row in values}
        for q in (q for q in queries if q['region']==region):
            key=tuple(round(v,7) for v in q['solver_xyz_m'])
            u=lookup.get(key)
            if u and all(math.isfinite(v) and abs(v)<1e10 for v in u):
                result[(region,q['id'])]={'world_xyz_m':q['world_xyz_m'],'U_world_m_s':transform(u,params['wind_from_deg'],inverse=True)}
    return result,last.name

def diagnostics(case):
    params=json.loads((case/'parameters.json').read_text())
    log=(case/'log.simpleFoam').read_text()
    residual={}
    terminal=re.split(r'^Time = ',log,flags=re.M)[-1]
    for field,value in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+)',terminal):
        residual[field]=max(residual.get(field,0.),float(value))
    checks={}
    for field,bound in params['residual_thresholds'].items():
        vals=[v for name,v in residual.items() if name==field or (field=='U' and name in ('Ux','Uy','Uz'))]
        checks[field]=bool(vals) and max(vals)<=bound
    flux={}
    for name in ('flux','outflux'):
        paths=list((case/'postProcessing'/name).glob('*/surfaceFieldValue.dat'))
        flux[name]=raw_rows(paths[-1])[-1][-1]
    imbalance=abs(flux['flux']+flux['outflux'])/abs(flux['flux'])
    monitor=case/'postProcessing/monitor'
    changes=[]
    if monitor.exists():
        folders=numeric_times(monitor)
        if len(folders)>=2:
            a,b=folders[-2:]
            av=raw_rows(a/'probes_U.xy');bv=raw_rows(b/'probes_U.xy')
            aa={tuple(round(v,7) for v in r[:3]):r[3:6] for r in av if all(abs(v)<1e10 for v in r[3:6])}
            bb={tuple(round(v,7) for v in r[:3]):r[3:6] for r in bv if all(abs(v)<1e10 for v in r[3:6])}
            keys=aa.keys()&bb.keys()
            if keys:changes=[{'iterations':[a.name,b.name],'samples':len(keys),'vector_change_RMSE_over_Uref':math.sqrt(sum(math.dist(aa[k],bb[k])**2 for k in keys)/len(keys))/params['Uref_m_s']}]
    settling=bool(changes) and changes[-1]['vector_change_RMSE_over_Uref']<=params['probe_change_RMSE_over_Uref_threshold']
    meshlog=(case/'log.checkMesh.standard').read_text()
    cellmatch=re.search(r'^\s*cells:\s+(\d+)',meshlog,re.M)
    cost={stage:json.loads((case/f'cost.{stage}.json').read_text()) for stage in ('prepare','mesh','solve','sample') if (case/f'cost.{stage}.json').exists()}
    mesh_checks={name:'Mesh OK.' in (case/f'log.checkMesh.{name}').read_text() for name in ('standard','quality')}
    cell_count=int(cellmatch[1]) if cellmatch else None
    mesh_checks['cell_cap']=cell_count is not None and cell_count<=params['max_cells']
    all_geometry=(case/'log.checkMesh').read_text()
    mesh_checks['allGeometry_or_concavity_only']='Mesh OK.' in all_geometry or ('Failed 1 mesh checks.' in all_geometry and 'Concave cells (using face planes) found' in all_geometry)
    stage_checks={stage:stage in cost and cost[stage].get('exit_code',0)==0 and not cost[stage].get('timed_out',False) for stage in ('prepare','mesh','solve','sample')}
    valid=all(stage_checks.values()) and all(mesh_checks.values()) and all(checks.values()) and imbalance<=params['flux_threshold'] and settling and 'SIMPLE solution converged' in log
    return {'case':case.name,'scenario':'numerical-abl','scope':params['scope'],'closure':params['closure'],'resolution':params['resolution'],
            'host':params['host'],'cells':cell_count,'mesh_checks':mesh_checks,
            'preparation_cost_included':'prepare' in cost,'stage_checks':stage_checks,
            'solver_final_iteration':re.findall(r'^Time = ([\d.]+)',log,re.M)[-1],
            'last_initial_residuals':residual,'frozen_residual_checks':checks,'solver_reported_convergence':'SIMPLE solution converged' in log,
            'boundary_flux_m3_s':flux,'normalized_flux_imbalance':imbalance,'probe_settling':changes,
            'frozen_numerical_checks_pass':valid,'status':'checked numerical result' if valid else 'provisional: frozen convergence/flux/output criteria unmet',
            'cost':cost,'complete_map_generation_elapsed_s':sum(v['elapsed_s'] for v in cost.values()),
            'complete_map_generation_cpu_s':sum(v['cpu_s'] for v in cost.values()) if all(v['cpu_s'] is not None for v in cost.values()) else None,
            'accounted_cpu_s_partial':sum(v['cpu_s'] or 0 for v in cost.values()),
            'peak_memory_accounting_complete':all(v['peak_cgroup_memory_bytes_including_cache'] is not None for name,v in cost.items() if name!='prepare'),
            'peak_cgroup_memory_bytes_including_cache':max((v['peak_cgroup_memory_bytes_including_cache'] or 0 for v in cost.values()),default=0),
            'physical_validation':'not performed: numerical-only ABL toy; experimental conditions intentionally differ'}

def evaluate_poc(case,reference=None):
    params=json.loads((case/'parameters.json').read_text());data,iteration=sampled_world(case)
    result=diagnostics(case);result['sample_iteration']=iteration
    queries=json.loads((case/'queries.json').read_text())
    result['sampling_counts']={region:{'total_requested':sum(q['region']==region for q in queries),'valid_fluid':sum(k[0]==region for k in data)} for region in ('city','domain','inflow')}
    profile=[]
    for q in (q for q in queries if q['region']=='inflow'):
        row=data.get(('inflow',q['id']))
        if row:
            z=q['world_xyz_m'][2]
            expected=params['ustar_m_s']/.41*math.log((z+params['z0_m'])/params['z0_m'])
            profile.append({'z_m':z,'mean_speed_m_s':math.sqrt(sum(v*v for v in row['U_world_m_s'])),'prescribed_speed_m_s':expected,
                            'speed_error_over_Uref':(math.sqrt(sum(v*v for v in row['U_world_m_s']))-expected)/params['Uref_m_s']})
    result['upstream_profile']=profile
    if reference:
        rp=json.loads((reference/'parameters.json').read_text())
        for key in ('closure','direction','Uref_m_s','z0_m','query_sha256','scaled_stl_sha256'):
            if params[key]!=rp[key]:raise SystemExit('Unmatched comparison parameter: '+key)
        ref,refiter=sampled_world(reference);rd=diagnostics(reference)
        comparison={'reference':reference.name,'reference_status':rd['status'],'reference_checked':rd['frozen_numerical_checks_pass'],
                    'scope':'Same-physics numerical field agreement, not physical validation','regions':{},
                    'region_definitions':{'city_low':'Fixed city ROI at z <=40 m AGL, including obstacle vicinity','city_downstream':'Fixed city ROI with solver-frame downstream y <=-100 m and z <=40 m AGL; geometric wake proxy, not identified physical wake'},
                    'complete_elapsed_ratio_descriptive':rd['complete_map_generation_elapsed_s']/result['complete_map_generation_elapsed_s'],
                    'complete_map_generation_speedup':rd['complete_map_generation_elapsed_s']/result['complete_map_generation_elapsed_s'] if rd['frozen_numerical_checks_pass'] else None}
        for region in ('city','city_low','city_downstream','domain'):
            base_region='city' if region.startswith('city') else region
            def in_region(key):
                if key[0]!=base_region:return False
                x,y,z=transform(ref[key]['world_xyz_m'],params['wind_from_deg'])
                return (z<=40 if region=='city_low' else y<=-100 and z<=40 if region=='city_downstream' else True)
            rk={k for k in ref if in_region(k)};keys=sorted(rk&data.keys())
            errors=sorted(math.dist(data[k]['U_world_m_s'],ref[k]['U_world_m_s']) for k in keys)
            p=.95*(len(errors)-1);lo=int(p);p95=errors[lo]+(errors[min(lo+1,len(errors)-1)]-errors[lo])*(p-lo) if errors else None
            rmse=math.sqrt(sum(e*e for e in errors)/len(errors)) if errors else None
            coverage=len(keys)/len(rk) if rk else 0.
            comparison['regions'][region]={'reference_valid':len(rk),'common_valid':len(keys),'common_fraction_of_reference_valid':coverage,
                'vector_RMSE_m_s':rmse,'vector_RMSE_over_Uref':rmse/params['Uref_m_s'] if rmse is not None else None,
                'pointwise_vector_error_P95_m_s':p95,'P95_over_Uref':p95/params['Uref_m_s'] if p95 is not None else None,
                'maximum_vector_error_m_s':max(errors) if errors else None}
        city=comparison['regions']['city']
        comparison['frozen_accuracy_and_coverage_targets_pass']=bool(city['vector_RMSE_over_Uref'] is not None and city['vector_RMSE_over_Uref']<=.05 and city['P95_over_Uref']<=.15 and city['common_fraction_of_reference_valid']>=.95)
        comparison['speedup_target_pass']=comparison['complete_map_generation_speedup'] is not None and comparison['complete_map_generation_speedup']>=2.
        comparison['overall_checked_target_pass']=bool(rd['frozen_numerical_checks_pass'] and result['frozen_numerical_checks_pass'] and comparison['frozen_accuracy_and_coverage_targets_pass'] and comparison['speedup_target_pass'])
        result['comparison']=comparison
    out=ROOT/'results'/case.name;out.mkdir(parents=True,exist_ok=True)
    (out/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    (out/'samples.json').write_text(json.dumps([{'region':k[0],'id':k[1],**v} for k,v in sorted(data.items())],separators=(',',':'))+'\n')
    print(json.dumps(result,indent=2))
    return result
