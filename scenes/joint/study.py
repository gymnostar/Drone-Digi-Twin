#!/usr/bin/env python3
"""Bounded coupled numerical transfer; reuses city builder and turbine force audit."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import re
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parent
SCENES=ROOT.parent
sys.path.insert(0,str(SCENES/'turbine/scripts'))
import prepare as turbine_prepare

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result
city=module('city_poc',SCENES/'buildings/poc.py')
runner=module('city_run',SCENES/'buildings/run.py')
source=module('turbine_poc',SCENES/'turbine/scripts/poc.py')
TRANSPORT_TOLERANCE=re.compile(
    r'("\(U\|k\|epsilon\)"\s*\{\s*solver smoothSolver;\s*'
    r'smoother symGaussSeidel;\s*tolerance )([0-9.eE+-]+)(;\s*relTol 0\.05;\s*\})'
)
RETAINED_RELAXATION='relaxationFactors {fields {p 0.3;} equations {U 0.6; k 0.6; epsilon 0.6;}}'
# OpenFOAM v2306 tutorials/incompressible/simpleFoam/pitzDaily/system/fvSolution
# (sha256 c562b517547fd9445cb070fb1826f9eaffad42ec2f34694fc727bea99f6b7dd1) relaxes
# only the equations, U and ".*" by 0.9, with no pressure field relaxation. 0.7 is
# the one pre-registered fallback in lookup/campaign.json simplec_relaxation_trials.
SIMPLEC_RELAXATION_FACTORS=(0.9,0.7)

def simplec_relaxation(factor):
    return f'relaxationFactors {{equations {{U {factor}; ".*" {factor};}}}}'

def simplec_relaxation_allowed(p):
    """The strict NNE reference and the NE and NW production baselines only."""
    strict_nne=(p.get('role')=='reference' and p.get('tight_transport',False)
                and p.get('direction')=='NNE')
    production=(p.get('role')=='baseline' and p.get('direction') in ('NE','NW')
                and bool(p.get('post_snap_agglomeration')))
    return (strict_nne or production) and not p.get('standard_simple_recovery',False)

def support_audit(case,p):
    result=source.support_audit(case,p)
    result['nominal_selected_volume_m3']=result.pop('selected_volume_m3')
    result['nominal_volume_over_geometric_cylinder']=result['nominal_selected_volume_m3']/(math.pi*13.5**2*6.75)
    result['volume_note']='Count times nominal1.25m cubed; snapped cells have small deformations. Exact selected-cell spatial geometry hash, not nominal volume, proves equal support across recipes. Cylinder center selection remains a finite support discretization.'
    return result

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')

def validate_solver(case,p):
    text=(case/'system/fvSolution').read_text()
    if p.get('reference_variant') is not None:
        raise ValueError('Historical reference variant is no longer supported')
    tight_transport=p.get('tight_transport',False)
    if tight_transport:
        if p['role']!='reference' or not p['common_reference_mesh_recipe']:
            raise ValueError('Tight transport requires the common-mesh reference')
        protocol_bytes=city.PROTOCOL.read_bytes()
        convergence=json.loads(protocol_bytes)['convergence']
        if (p['protocol_sha256']!=hashlib.sha256(protocol_bytes).hexdigest()
                or p['residual_thresholds']!=convergence['reference_residual_control']
                or p['flux_threshold']!=convergence['relative_flux_imbalance_max']
                or p['probe_change_RMSE_over_Uref_threshold']!=
                    convergence['terminal_probe_vector_change_RMSE_over_Uref_max']):
            raise ValueError('Tight transport requires the frozen reference protocol')
    transport=TRANSPORT_TOLERANCE.findall(text)
    expected_tolerance=1e-10 if tight_transport else 1e-7
    if len(transport)!=1 or float(transport[0][1])!=expected_tolerance:
        raise ValueError('Unexpected inner transport solver settings')
    block=re.search(r'residualControl \{([^}]+)\}',text)[1]
    actual={key:float(value) for key,value in re.findall(r'(p|U|k|epsilon) ([0-9.eE+-]+);',block)}
    if actual!=p['residual_thresholds']:raise ValueError('Frozen outer residual controls changed')
    simplec=p.get('simplec_conventional_relaxation')
    if simplec is not None and (type(simplec) is not float or simplec not in SIMPLEC_RELAXATION_FACTORS):
        raise ValueError('Invalid simplec_conventional_relaxation setting')
    if simplec is not None and not simplec_relaxation_allowed(p):
        raise ValueError('Conventional SIMPLEC relaxation requires the strict NNE reference or the NE or NW production baseline')
    if text.count('relaxationFactors')!=1:raise ValueError('Retained relaxation factors changed')
    relaxation='relaxationFactors'+text.split('relaxationFactors',1)[1]
    if simplec is not None:
        if relaxation.strip()!=simplec_relaxation(simplec):
            raise ValueError('Conventional SIMPLEC relaxation factors changed')
    else:
        expected={'p':.3,'U':.6,'k':.6,'epsilon':.6}
        actual={key:float(value) for key,value in re.findall(r'(p|U|k|epsilon) ([0-9.eE+-]+);',relaxation)}
        if actual!=expected or '"' in relaxation:raise ValueError('Retained relaxation factors changed')
    standard_simple = p.get('standard_simple_recovery', False)
    if not isinstance(standard_simple, bool):
        raise ValueError('Invalid standard_simple_recovery setting')
    if standard_simple and not (tight_transport and p['direction'] == 'NNE'):
        raise ValueError('Standard SIMPLE recovery requires the strict NNE reference')
    expected = 'yes' if (p['role'] != 'baseline' or simplec is not None) and not standard_simple else 'no'
    if re.findall(r'\bconsistent\s+(\w+)\s*;', text) != [expected]:
        raise ValueError('Unmatched SIMPLE convention')

HUB_WORLD=(0,-280,32.1)
DEFAULT_RECIPES={'N':'production','NbE':'strict','NNE':'strict','NE':'production','NW':'strict'}
JOINT_DIRECTIONS=tuple(DEFAULT_RECIPES)
RECIPES={
    'strict': ('reference', True, True, False),
    'production': ('baseline', True, False, True),
}


def rotor_point(north_point, degrees):
    """Keep the world hub fixed and yaw north-case offsets with the inflow."""
    if degrees==0:
        return north_point
    hub=city.transform(HUB_WORLD,degrees)
    x,y,z=north_point
    return (hub[0]+x-HUB_WORLD[0],hub[1]+y-HUB_WORLD[1],z)


def foam_point(point):
    return '('+' '.join(map(str,point))+')'


def build(name,role,common_mesh=False,tight_transport=False,agglomerate=False,direction='N',standard_simple=False,simplec_relaxation_factor=None):
    started=time.monotonic();cpu=time.process_time()
    if Path(name).name!=name:raise ValueError('One fresh case name required')
    if tight_transport and (role!='reference' or not common_mesh):
        raise ValueError('Tight transport requires --role reference --common-mesh')
    if agglomerate and (role!='baseline' or not common_mesh):
        raise ValueError('Agglomeration requires --role baseline --common-mesh')
    if standard_simple and not (direction == 'NNE' and tight_transport):
        raise ValueError('Standard SIMPLE recovery requires the strict NNE reference')
    if simplec_relaxation_factor is not None:
        simplec_relaxation_factor=float(simplec_relaxation_factor)
        if not simplec_relaxation_allowed({'role':role,'direction':direction,'tight_transport':tight_transport,
                                           'post_snap_agglomeration':agglomerate,'standard_simple_recovery':standard_simple}):
            raise ValueError('Conventional SIMPLEC relaxation requires the strict NNE reference or the NE or NW production baseline')
    case=ROOT/'runs'/name
    if direction not in JOINT_DIRECTIONS:
        raise ValueError('Unsupported coupled direction '+direction)
    degrees=city.DIRECTIONS[direction]
    args=SimpleNamespace(case=name,closure='realizableKE',direction=direction,oriented_world_boxes=True,resolution='reference',empty=False,background_m=20.,mesh_repair=False,solver_repair=role!='baseline',consistent_simple=role!='baseline' and not standard_simple)
    city.build_poc(args, case_root=ROOT/'runs')
    p=json.loads((case/'parameters.json').read_text())
    p.update(joint=True,role=role,scope='Coupled full Niigata geometry and hypothetical full-size SWiFT rotor placement; numerical transfer only, no interaction observations.',max_cells=500000)
    if standard_simple:
        p['standard_simple_recovery'] = True
    p['city_refinement_frame']='Fixed world-space floor, district and repair boxes.'
    protocol=json.loads(city.PROTOCOL.read_text())
    thresholds=protocol['convergence']['reference_residual_control' if role=='reference' else 'cheap_residual_control']
    p['residual_thresholds']=thresholds
    fv=case/'system/fvSolution';text=fv.read_text()
    residual='residualControl {'+' '.join(f'{field} {val};' for field,val in thresholds.items())+'}'
    text=re.sub(r'residualControl \{[^}]*\}',residual,text)
    if tight_transport:
        transport=TRANSPORT_TOLERANCE.findall(text)
        if len(transport)!=1 or float(transport[0][1])!=1e-7:
            raise ValueError('Expected one original inner transport tolerance')
        text=TRANSPORT_TOLERANCE.sub(r'\g<1>1e-10\g<3>',text)
        p['tight_transport']=True
    if simplec_relaxation_factor is not None:
        if text.count(RETAINED_RELAXATION)!=1 or len(re.findall(r'\bconsistent (?:yes|no);',text))!=1:
            raise ValueError('Expected the retained relaxation and one SIMPLE convention')
        text=text.replace(RETAINED_RELAXATION,simplec_relaxation(simplec_relaxation_factor))
        text=re.sub(r'\bconsistent (?:yes|no);','consistent yes;',text)
        p['simplec_conventional_relaxation']=simplec_relaxation_factor
    fv.write_text(text)
    wake_level=3 if common_mesh else {'reference':3,'baseline':2,'candidate':1}[role]
    support_box=((-20.25,-286.75,6.75),(20.25,-273.25,60.75))
    wake_box=((-40.5,-455.5,0),(40.5,-253,81))
    support_box=tuple(rotor_point(p,degrees) for p in support_box)
    wake_box=tuple(rotor_point(p,degrees) for p in wake_box)
    snappy=case/'system/snappyHexMeshDict';text=snappy.read_text()
    text=text.replace('geometry {\n','geometry {\n sourceSupport {type searchableBox; '+city.box_string(*support_box)+'}\n rotorWake {type searchableBox; '+city.box_string(*wake_box)+'}\n',1)
    text=text.replace('refinementRegions {','refinementRegions { sourceSupport {mode inside; levels ((1E15 4));} rotorWake {mode inside; levels ((1E15 '+str(wake_level)+'));}',1)
    # Extended mesh checks locate these near-ground building corners in world metres.
    repairs_by_direction={
        'NNE': (((80,142.5,0),(92.5,155,7.5)), ((-125,-147.5,0),(-110,-135,7.5))),
        'NbE': (((10,17.5,0),(25,30,7.5)),),
    }
    if direction in repairs_by_direction:
        repairs=repairs_by_direction[direction]
        p['directional_mesh_repairs_world_m']=repairs
        for i,(low,high) in enumerate(repairs):
            name=f'buildingCorner{i}'
            geometry=city.world_refinement_box(low,high,degrees,oriented=True)
            text=text.replace('geometry {\n',f'geometry {{\n {name} {{{geometry}}}\n',1)
            text=text.replace('refinementRegions {',
                              f'refinementRegions {{ {name} {{mode inside; levels ((1E15 4));}}',1)
    if agglomerate:
        protected=tuple(rotor_point(p,degrees) for p in ((-26,-292,1.5),(26,-268,66)))
        city.write(case/'system/agglomerationBox',city.box_string(*protected))
        shutil.copytree(ROOT/'agglomeration',case/'system/agglomeration')
        p['post_snap_agglomeration']='complete interior snapped sibling octets'
    snappy.write_text(text)
    lock=json.loads((SCENES/'turbine/reference/data-sources.json').read_text())
    for item in lock['files']:
        path=SCENES/'turbine/data/raw'/item['path']
        if path.stat().st_size!=item['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('Turbine input hash mismatch')
    cp,ct=turbine_prepare.rotor_coefficients(SCENES/'turbine/data/raw',8.7)
    area=math.pi*13.5**2
    p1=rotor_point((0,-283.375,32.1),degrees)
    p2=rotor_point((0,-276.625,32.1),degrees)
    monitor=rotor_point((0,-215,32.1),degrees)
    city.write(case/'system/topoSetDict',f'actions ({{name rotor; type cellSet; action new; source cylinderToCell; p1 {foam_point(p1)}; p2 {foam_point(p2)}; radius 13.5;}});')
    city.write(case/'constant/fvOptions',f'rotor {{type actuationDiskSource; variant variableScaling; selectionMode cellSet; cellSet rotor; diskArea {area:.12g}; diskDir (0 -1 0); sink true; Cp {cp:.12g}; Ct {ct:.12g}; monitorMethod points; monitorCoeffs {{points ({foam_point(monitor)});}} writeToFile true; writeFileStart 0; writeFileEnd 1000000;}}')
    settling_points=[rotor_point((0.,y,32.1),degrees) for y in (-215.,-334.,-361.,-388.,-415.)]
    monitor_text=' '.join(foam_point(point) for point in settling_points)
    control=case/'system/controlDict';control.write_text(control.read_text().replace('points (','points ('+monitor_text+' ',1))
    p['additional_terminal_settling_points_m']=settling_points
    support_config={'base_m':20,'origin_m':[-600,-1200,0],'level':4,'cell_m':1.25,'box':support_box,'cylinder':{'p1':p1,'p2':p2,'radius':13.5}}
    p['common_reference_mesh_recipe']=(common_mesh or role=='reference') and not agglomerate
    p.update(
        rotor_D_m=27., rotor_center_m=HUB_WORLD, rotor_center_solver_m=city.transform(HUB_WORLD,degrees),
        monitor_m=city.transform(monitor,degrees,inverse=True), monitor_solver_m=monitor,
        diskDir=[0,-1,0], disk_direction_world=city.transform((0,-1,0),degrees,inverse=True),
        Cp=cp, Ct=ct, area_m2=area, density_kg_m3=1.02,
        source_cell_m=1.25, source_refinement_level=4, wake_level=wake_level,
        force_support={
            'configuration_sha256':hashlib.sha256(json.dumps(support_config,sort_keys=True).encode()).hexdigest(),
            'configuration':support_config,
        },
        rotor_coefficient_convention='Locked v3 Cp/Ct interpolated at prescribed8.7m/s; native variableScaling uses solved local monitor vector; no joint fitting.',
        transfer_limitation='Separate30m cheap recipe cannot share the20m source lattice. This adapted20m candidate is not an exact replay or established optimized method.',
        cell_limit_adjustment={
            'original_cheap_target':150000,'joint_cap':500000,
            'reason':'One allowed adjustment: full city geometry accepted20m background plus common1.25m force support. Cheap150k target remains reported unmet if exceeded.',
        },
        physical_validation='No joint observations; numerical agreement is not physical or flight validation.',
    )
    queries=json.loads((case/'queries.json').read_text())
    for region,points in (
        ('rotor_wake',[(x/2,y/2,z/2) for z in range(10,131,5) for y in range(-830,-613,5) for x in range(-54,55,5)]),
        ('flight',[(float(x),float(y),float(z)) for z in (20,40,60) for y in range(-450,-219,5) for x in range(-50,51,5)]),
        ('rotor_monitor',[(0.,-215.,32.1)])):
        for i,point in enumerate(points):
            if not degrees:
                world=solver=point
            elif region=='flight':
                world=point
                solver=city.transform(world,degrees)
            else:
                solver=rotor_point(point,degrees)
                world=city.transform(solver,degrees,inverse=True)
            queries.append({'region':region,'id':i,'world_xyz_m':world,'solver_xyz_m':solver})
    sets=[]
    for region in dict.fromkeys(q['region'] for q in queries):
        pts=' '.join('('+' '.join(map(str,q['solver_xyz_m']))+')' for q in queries if q['region']==region)
        sets.append(f'{region} {{type cloud; axis xyz; points ({pts});}}')
    city.write(case/'system/roiSampling','functions {roi {type sets; libs (sampling); writeControl timeStep; writeInterval 1; interpolationScheme cellPoint; setFormat raw; fields (U k); sets {'+'\n'.join(sets)+'}}}')
    (case/'queries.json').write_text(json.dumps(queries,separators=(',',':'))+'\n')
    p['query_sha256']=hashlib.sha256((case/'queries.json').read_bytes()).hexdigest()
    validate_solver(case,p)
    dump(case/'parameters.json',p)
    preparation_cost={'stage':'prepare','elapsed_s':time.monotonic()-started,'cpu_s':time.process_time()-cpu,'peak_cgroup_memory_bytes_including_cache':None,'scope':'Full source hash checks, case dictionaries, geometry transform and frozen query preparation; acquisition excluded'}
    if tight_transport:
        preparation_cost['tight_transport']=True
    dump(case/'cost.prepare.json',preparation_cost)
    print(json.dumps({'case':str(case),'role':role,'queries':len(queries),'max_cells':p['max_cells']}))

def matched_reference(parameters, reference):
    if reference is None:
        suffix='tight-transport' if parameters.get('post_snap_agglomeration') else 'fixed'
        reference=ROOT/'runs'/f"poc-{parameters['direction']}-reference-{suffix}"
    rp=json.loads((reference/'parameters.json').read_text())
    if (rp['role']!='reference' or rp['direction']!=parameters['direction']
            or rp['wind_from_deg']!=parameters['wind_from_deg']):
        raise ValueError('Joint comparison requires a reference with the same wind direction')
    return reference,rp


def mesh(case,reference=None):
    outer_start=time.monotonic();outer_cpu=time.process_time()
    p=json.loads((case/'parameters.json').read_text())
    validate_solver(case,p)
    if p.get('post_snap_agglomeration'):
        reference,_=matched_reference(p,reference)
    commands='blockMesh > log.blockMesh 2>&1 && snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 && '
    if p.get('post_snap_agglomeration'):
        commands+=(
            '(cd system/agglomeration && time wmake -j3) > log.agglomeration-build 2>&1 && '
            '/tmp/agglomerate > log.agglomeration 2>&1 && '
            'rm -f constant/polyMesh/{surfaceIndex,refinementHistory} && '
            'rm -rf constant/polyMesh/sets && '
        )
    commands+='checkMesh -allTopology -allGeometry > log.checkMesh 2>&1 && checkMesh > log.checkMesh.standard 2>&1 && checkMesh -meshQuality > log.checkMesh.quality 2>&1 && topoSet > log.topoSet 2>&1'
    runner.run_stage(case,'mesh',commands,600)
    runner.accept_mesh(case);runner.check_budget(case,p)
    audit=support_audit(case,p)
    dump(case/'support-audit.json',audit)
    if p.get('post_snap_agglomeration'):
        expected=json.loads((reference/'support-audit.json').read_text())
        keys=('configuration_sha256','selected_cells','spatial_cell_vertices_sha256')
        if not all(audit[key]==expected[key] for key in keys):
            raise ValueError('Agglomeration changed exact rotor support')
    if not audit['all_selected_cells_at_fixed_level']:raise ValueError('Source support not uniformly fixed level')
    mesh_cost=json.loads((case/'cost.mesh.json').read_text())
    total_mesh_elapsed=time.monotonic()-outer_start
    dump(case/'cost.host-mesh-audit.json',{'elapsed_s':total_mesh_elapsed-mesh_cost['elapsed_s'],'cpu_s':time.process_time()-outer_cpu,'scope':'Host dictionary/mesh acceptance and exact support hash; excludes Docker stage already accounted'})
    if p.get('post_snap_agglomeration') and total_mesh_elapsed>600:
        raise ValueError('Agglomeration exceeded the mesh budget including support audit')

def sample(case):
    commands=''
    if not (case/'log.reconstructPar').exists():commands+='reconstructPar -latestTime > log.reconstructPar 2>&1 && '
    runner.run_stage(case,'sample',commands+'simpleFoam -postProcess -dict system/roiSampling -latestTime > log.roiSampling 2>&1',120)

def raw_vectors(expected, rows, *, allow_other_points=False):
    """Join cloud rows by coordinate, allowing only raw writer round-off.

    OpenFOAM raw sets can print six significant digits, including rotated
    positions. Index the rounding interval per axis in one-metre buckets.
    Missing rows stay missing. Ambiguous or repeated rows cannot supply an ID.
    """
    buckets={}
    for key,point in expected.items():
        intervals=[]
        for x in point:
            half_unit=0.5*10**(math.floor(math.log10(abs(x)))-5) if x else 0.
            eps=half_unit+1e-9
            intervals.append((x-eps,x+eps))
        bins=[range(math.floor(a),math.floor(b)+1) for a,b in intervals]
        for index in itertools.product(*bins):
            buckets.setdefault(index,[]).append((key,intervals))
    vectors={}
    seen=set()
    for row in rows:
        if len(row)!=6 or not all(math.isfinite(v) for v in row[:3]):
            raise ValueError('Malformed raw cloud sample')
        index=tuple(math.floor(v) for v in row[:3])
        matches=[key for key,bounds in buckets.get(index,())
                 if all(lo<=v<=hi for v,(lo,hi) in zip(row[:3],bounds))]
        if not matches and allow_other_points:
            continue
        if len(matches)!=1 or matches[0] in seen:
            raise ValueError(f'Raw cloud coordinate has no unique query ID: {row[:3]}')
        key=matches[0]
        seen.add(key)
        if all(math.isfinite(v) and abs(v)<1e10 for v in row[3:]):
            vectors[key]=row[3:]
    return vectors


def samples(case):
    p=json.loads((case/'parameters.json').read_text())
    queries=json.loads((case/'queries.json').read_text())
    folder=city.numeric_times(case/'postProcessing/roi')[-1];data={}
    for region in dict.fromkeys(q['region'] for q in queries):
        combined=folder/f'{region}_k_U.xy'
        rows=[r[:3]+r[4:7] for r in city.raw_rows(combined)] if combined.exists() else city.raw_rows(folder/f'{region}_U.xy')
        requested={q['id']:q for q in queries if q['region']==region}
        expected={key:q['solver_xyz_m'] for key,q in requested.items()}
        for key,u in raw_vectors(expected,rows).items():
            data[(region,key)]={'xyz':requested[key]['world_xyz_m'],
                                'U':city.transform(u,p['wind_from_deg'],inverse=True)}
    return data,queries,folder.name


def force_audit(case,p,data):
    rows=[]
    for line in (case/'postProcessing/rotor/0/actuationDiskSource.dat').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):rows.append([float(v) for v in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',line)])
    sample_time=float(city.numeric_times(case/'postProcessing/roi')[-1].name)
    row=max((r for r in rows if r[0]<=sample_time),key=lambda r:r[0]);monitor=row[1:4];disk=row[6:9];thrust=row[11]
    normal=p['diskDir']
    disk_axial=sum(v*n for v,n in zip(disk,normal))
    monitor_axial=sum(v*n for v,n in zip(monitor,normal))
    disk_squared=sum(v*v for v in disk)
    expected=(.5*p['area_m2']*p['Ct']*sum(v*v for v in monitor)*disk_axial**2/disk_squared
              if disk_squared>0 else 0.)
    relative=abs(thrust-expected)/expected if expected>0 else None
    valid=(row[0]==sample_time and ('rotor_monitor',0) in data and disk_axial>0 and monitor_axial>0
           and thrust>0 and relative is not None and relative<1e-8)
    return {'passed':valid,'force_iteration':row[0],'sample_iteration':sample_time,
            'monitor_is_valid_fluid_sample':('rotor_monitor',0) in data,
            'monitor_vector_m_s':city.transform(monitor,p['wind_from_deg'],inverse=True),
            'disk_mean_vector_m_s':city.transform(disk,p['wind_from_deg'],inverse=True),
            'vector_frame':'world east/north/up','kinematic_thrust':thrust,
            'physical_thrust_N':p['density_kg_m3']*thrust,'formula_relative_error':relative,
            'sink_convention':'Native sink opposes flow along the yawed disk normal. The thrust audit uses velocity projected onto that normal in solver coordinates.'}


def wake_settling(case,p):
    folders=city.numeric_times(case/'postProcessing/monitor')
    if len(folders)<2:return {'passed':False,'reason':'Less than two terminal samples'}
    a,b=folders[-2:]
    expected=dict(enumerate(p['additional_terminal_settling_points_m'][1:]))
    def points(folder):
        return raw_vectors(expected,city.raw_rows(folder/'probes_U.xy'),allow_other_points=True)
    aa,bb=points(a),points(b);keys=aa.keys()&bb.keys()
    error=math.sqrt(sum(math.dist(aa[k],bb[k])**2 for k in keys)/len(keys))/p['Uref_m_s'] if keys else None
    return {'iterations':[a.name,b.name],'valid_wake_probes':len(keys),'vector_change_RMSE_over_Uref':error,'passed':len(keys)==4 and error is not None and error<=p['probe_change_RMSE_over_Uref_threshold']}

def elapsed_recorded(cost, *, allow_zero=False):
    elapsed=cost.get('elapsed_s')
    return type(elapsed) in (int,float) and math.isfinite(elapsed) and (
        elapsed>=0 if allow_zero else elapsed>0)

def generation_timing(case, diagnostics):
    audits={path.stem.removeprefix('cost.'):json.loads(path.read_text())
            for path in case.glob('cost.host-*.json')}
    elapsed=diagnostics['complete_map_generation_elapsed_s']
    diagnostics.update(
        recorded_pipeline_elapsed_s=elapsed,
        supplemental_host_verification_cost=audits,
        complete_map_generation_elapsed_s=elapsed+sum(
            cost['elapsed_s'] for cost in audits.values() if elapsed_recorded(cost,allow_zero=True)),
        generation_timing_complete=(
            all(elapsed_recorded(diagnostics['cost'].get(stage,{}))
                for stage in ('prepare','mesh','solve','sample'))
            and all(elapsed_recorded(audits.get(f'host-{stage}-audit',{}),allow_zero=True)
                    for stage in ('mesh','solve'))
            and all(elapsed_recorded(cost,allow_zero=True) for cost in audits.values())),
        generation_timing_note=(
            'Recorded stage wall time plus measured automated host verification. '
            'Incomplete totals omit unrecorded intervals. The CPU totals cover the four runner '
            'stages; host audit CPU is reported separately. Comparative scientific evaluation '
            'and manual implementation work excluded.'))

def cost_environment_matches(a,b):
    resources=json.loads(city.PROTOCOL.read_text())['resources']
    stages=[result['cost'].get(stage,{}) for result in (a,b) for stage in ('mesh','solve','sample')]
    host=a.get('host')
    platform=stages[0].get('host_platform')
    expected={'host':host,'host_platform':platform,'image':runner.IMAGE,
              'cpu_limit':resources['cpu_quota'],'memory_limit_bytes':resources['memory_limit_GiB']*1024**3}
    return bool(host and platform and host==b.get('host')
                and all(stage.get(key)==value for stage in stages for key,value in expected.items()))

def strict_reference_configuration(case, parameters):
    protocol_bytes=city.PROTOCOL.read_bytes()
    convergence=json.loads(protocol_bytes)['convergence']
    matches=(parameters.get('role')=='reference'
             and parameters.get('protocol_sha256')==hashlib.sha256(protocol_bytes).hexdigest()
             and parameters.get('residual_thresholds')==convergence['reference_residual_control']
             and parameters.get('flux_threshold')==convergence['relative_flux_imbalance_max']
             and parameters.get('probe_change_RMSE_over_Uref_threshold')==
                 convergence['terminal_probe_vector_change_RMSE_over_Uref_max'])
    if matches:
        try:
            validate_solver(case,parameters)
        except (OSError,TypeError,ValueError):
            return False
    return matches

def in_building_wake(world_xyz,degrees):
    _,downwind_coordinate,z=city.transform(world_xyz,degrees)
    return downwind_coordinate<=-100 and z<=40


def summarize_vector_errors(errors, reference_valid, uref):
    """Compare common valid vectors while keeping the reference-valid denominator."""
    errors = sorted(errors)
    count = len(errors)
    rmse = p95 = None
    if count:
        rmse = math.sqrt(sum(e*e for e in errors) / count)
        position = .95 * (count - 1)
        lo = int(position)
        p95 = errors[lo] + (errors[math.ceil(position)] - errors[lo]) * (position - lo)
    return {
        'reference_valid': reference_valid, 'common_valid': count,
        'coverage': count / reference_valid if reference_valid else 0,
        'vector_RMSE_m_s': rmse, 'P95_m_s': p95,
        'vector_RMSE_over_Uref': rmse / uref if rmse is not None else None,
        'P95_over_Uref': p95 / uref if p95 is not None else None,
    }


def region_errors(data, ref, keys, uref):
    errors = (math.dist(data[k]['U'], ref[k]['U']) for k in keys if k in data)
    return summarize_vector_errors(errors, len(keys), uref)


def descriptive_city_errors(data, ref, queries, uref):
    requested_by_height={}
    for query in queries:
        if query['region']=='city':
            x,y,z=query['world_xyz_m']
            requested_by_height.setdefault(z,{})[x,y]=('city',query['id'])

    heights={}
    narrow_gap=set()
    for height, requested in sorted(requested_by_height.items()):
        valid={key for key in requested.values() if key in ref}
        invalid={xy for xy,key in requested.items() if key not in ref}
        heights[f'{height:g}']={'requested':len(requested),**region_errors(data,ref,valid,uref)}
        for (x,y),key in requested.items():
            if key in valid and (
                (any((x-d,y) in invalid for d in (5,10))
                 and any((x+d,y) in invalid for d in (5,10)))
                or (any((x,y-d) in invalid for d in (5,10))
                    and any((x,y+d) in invalid for d in (5,10)))
            ):
                narrow_gap.add(key)
    return {
        'scope':'Descriptive diagnostics outside the four selection gates.',
        'city_by_height_agl_m':heights,
        'narrow_gap_proxy_definition':(
            'Reference-valid city points bracketed along x or y at the same requested height '
            'by reference-invalid requested city points within 10 m on each side. '
            'Reference sampling validity is a geometric proxy.'),
        'narrow_gap_proxy':region_errors(data,ref,narrow_gap,uref),
    }

def evaluate(case,reference=None):
    p=json.loads((case/'parameters.json').read_text());d=city.diagnostics(case);data,queries,iteration=samples(case)
    support=support_audit(case,p);force=force_audit(case,p,data)
    generation_timing(case,d)
    d.update(scope=p['scope'],joint=True,role=p['role'],direction=p['direction'],wind_from_deg=p['wind_from_deg'],sample_iteration=iteration,force_support=support,force_audit=force,cheap150k_target_met=d['cells']<=150000,cell_limit_adjustment=p['cell_limit_adjustment'],transfer_limitation=p['transfer_limitation'])
    d['sample_matches_final_iteration']=float(iteration)==float(d['solver_final_iteration'])
    d['rotor_wake_settling']=wake_settling(case,p)
    d['checked_coupled_result']=d['frozen_numerical_checks_pass'] and d['sample_matches_final_iteration'] and d['rotor_wake_settling']['passed'] and force['passed'] and support['all_selected_cells_at_fixed_level']
    d['sampling_counts']={region:{'requested':sum(q['region']==region for q in queries),'valid':sum(k[0]==region for k in data)} for region in dict.fromkeys(q['region'] for q in queries)}
    if reference:
        reference,rp=matched_reference(p,reference);rd=city.diagnostics(reference);ref,ref_queries,ref_iteration=samples(reference);rs=support_audit(reference,rp);rf=force_audit(reference,rp,ref)
        for key in ('protocol_sha256','query_sha256','scaled_stl_sha256','Cp','Ct','closure','Uref_m_s','z0_m'):
            if p[key]!=rp[key]:raise ValueError('Unmatched joint '+key)
        for path in ('0/U','0/p','0/k','0/epsilon','0/nut','constant/fvOptions',
                     'constant/transportProperties','constant/turbulenceProperties','system/fvSchemes'):
            if (case/path).read_bytes()!=(reference/path).read_bytes():
                raise ValueError('Unmatched coupled case input '+path)
        if queries!=ref_queries:
            raise ValueError('Unmatched joint query records')
        support_match=all(support[k]==rs[k] for k in ('configuration_sha256','selected_cells','spatial_cell_vertices_sha256'))
        strict_reference=strict_reference_configuration(reference,rp)
        reference_aligned=float(ref_iteration)==float(rd['solver_final_iteration'])
        checked=(strict_reference and rd['frozen_numerical_checks_pass'] and reference_aligned
                 and wake_settling(reference,rp)['passed'] and rf['passed'] and rs['all_selected_cells_at_fixed_level'])
        generation_timing(reference,rd)
        timing_complete=d['generation_timing_complete'] and rd['generation_timing_complete']
        environment_matches=cost_environment_matches(d,rd)
        elapsed_ratio=(rd['complete_map_generation_elapsed_s']/d['complete_map_generation_elapsed_s']
                       if elapsed_recorded({'elapsed_s':rd['complete_map_generation_elapsed_s']})
                       and elapsed_recorded({'elapsed_s':d['complete_map_generation_elapsed_s']}) else None)
        comparison={
            'reference':reference.name,'reference_checked':checked,
            'reference_strict_configuration':strict_reference,
            'reference_sample_iteration':ref_iteration,
            'reference_solver_final_iteration':rd['solver_final_iteration'],
            'reference_sample_matches_final_iteration':reference_aligned,
            'reference_complete_map_generation_elapsed_s':rd['complete_map_generation_elapsed_s'],
            'reference_generation_timing_complete':rd['generation_timing_complete'],
            'reference_supplemental_host_verification_cost':rd['supplemental_host_verification_cost'],
            'both_generation_timings_complete':timing_complete,
            'cost_environment_matches':environment_matches,
            'support_exact_match':support_match,'regions':{},
            'elapsed_ratio_descriptive':elapsed_ratio,
            'speedup':elapsed_ratio if (p['role'] in ('baseline','candidate') and checked and d['checked_coupled_result'] and support_match
                                        and timing_complete and environment_matches) else None,
        }
        for region in ('city','building_wake','rotor_wake','flight','domain'):
            base='city' if region=='building_wake' else region
            keys={k for k in ref if k[0]==base and (region!='building_wake' or in_building_wake(ref[k]['xyz'],p['wind_from_deg']))}
            metrics=region_errors(data,ref,keys,rp['Uref_m_s'])
            rmse=metrics['vector_RMSE_over_Uref']
            metrics['targets_pass']=rmse is not None and rmse<=.05 and metrics['P95_over_Uref']<=.15 and metrics['coverage']>=.95
            comparison['regions'][region]=metrics
        comparison['descriptive_diagnostics']=descriptive_city_errors(data,ref,ref_queries,rp['Uref_m_s'])
        comparison['all_declared_ROI_targets_pass']=all(comparison['regions'][k]['targets_pass'] for k in ('city','building_wake','rotor_wake','flight'))
        comparison['overall_target_pass']=bool(checked and d['checked_coupled_result'] and support_match and comparison['all_declared_ROI_targets_pass'] and comparison['speedup'] is not None and comparison['speedup']>=2)
        d['comparison']=comparison
    out=ROOT/'results'/case.name;dump(out/'metrics.json',d);dump(out/'parameters.json',p)
    (out/'samples.json').write_text(json.dumps([{'region':key[0],'id':key[1],**value} for key,value in sorted(data.items())],separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('cost',)},indent=2))
    return d

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--common-mesh',action='store_true',help='Start from the reference wake refinement and rotor-support lattice')
    parser.add_argument('stage',choices=['build','mesh','solve','sample','evaluate'])
    parser.add_argument('--case',required=True)
    parser.add_argument('--direction',choices=JOINT_DIRECTIONS,help='Meteorological wind-from direction for build (default N)')
    parser.add_argument('--recipe',choices=RECIPES,help='Mesh and solver preset. Defaults to production for N and NE, strict for NbE, NNE and NW (no checked reference).')
    parser.add_argument('--role',choices=['reference','baseline','candidate'])
    parser.add_argument('--reference',help='Reference case for the mesh/solve support audit or field evaluation')
    parser.add_argument('--tight-transport',action='store_true',
                        help='Build a common-mesh reference with inner transport absolute tolerance 1e-10')
    parser.add_argument('--agglomerate',action='store_true',
                        help='Build a baseline with interior sibling agglomeration after common-mesh snapping')
    parser.add_argument('--standard-simple',action='store_true',
                        help='Reproduce the strict NNE standard-SIMPLE trial that exceeded the 600 s solve limit')
    parser.add_argument('--simplec-relaxation',type=float,choices=SIMPLEC_RELAXATION_FACTORS,
                        help='SIMPLEC with unrelaxed pressure and U/".*" equation relaxation: 0.9 as in the '
                             'v2306 pitzDaily tutorial, 0.7 the pre-registered fallback. Strict NNE or NE/NW production only')
    args=parser.parse_args()
    if args.stage!='build' and (args.direction or args.recipe or args.simplec_relaxation is not None):
        parser.error('--direction, --recipe and --simplec-relaxation apply to build')
    if args.stage=='build' and not (args.recipe or args.role):
        args.recipe=DEFAULT_RECIPES[args.direction or 'N']
    if args.recipe:
        if args.role or args.common_mesh or args.tight_transport or args.agglomerate:
            parser.error('--recipe already selects role, mesh and solver flags')
        args.role,args.common_mesh,args.tight_transport,args.agglomerate=RECIPES[args.recipe]
    if args.tight_transport and (args.stage!='build' or args.role!='reference' or not args.common_mesh):
        parser.error('--tight-transport requires build --role reference --common-mesh')
    if args.agglomerate and (args.stage!='build' or args.role!='baseline' or not args.common_mesh):
        parser.error('--agglomerate requires build --role baseline --common-mesh')
    if args.standard_simple and not (args.stage=='build' and args.tight_transport and args.direction=='NNE'):
        parser.error('--standard-simple requires build --direction NNE with the strict recipe')
    case=ROOT/'runs'/args.case
    reference=ROOT/'runs'/args.reference if args.reference else None
    if args.stage=='build':build(args.case,args.role,args.common_mesh,args.tight_transport,args.agglomerate,args.direction or 'N',args.standard_simple,args.simplec_relaxation)
    elif args.stage=='mesh':mesh(case,reference)
    elif args.stage=='solve':
        validate_solver(case,json.loads((case/'parameters.json').read_text()))
        gate_start=time.monotonic();gate_cpu=time.process_time()
        p=json.loads((case/'parameters.json').read_text())
        if p['role']!='reference':
            reference,rp=matched_reference(p,reference)
            a=support_audit(case,p);b=support_audit(reference,rp)
            if not all(a[key]==b[key] for key in ('configuration_sha256','selected_cells','spatial_cell_vertices_sha256')):raise ValueError('Exact force support gate failed; solve prohibited')
        dump(case/'cost.host-solve-audit.json',{'elapsed_s':time.monotonic()-gate_start,'cpu_s':time.process_time()-gate_cpu,'scope':'Mandatory exact support comparison to reference before solve; original reference requires no cross-recipe gate'})
        runner.accept_mesh(case);runner.check_budget(case,json.loads((case/'parameters.json').read_text()));runner.run_stage(case,'solve','decomposePar -force > log.decomposePar 2>&1 && mpirun --allow-run-as-root -np 3 simpleFoam -parallel > log.simpleFoam 2>&1 && reconstructPar -latestTime > log.reconstructPar 2>&1',600)
    elif args.stage=='sample':sample(case)
    else:evaluate(case,reference)
if __name__=='__main__':main()
