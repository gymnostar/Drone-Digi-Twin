#!/usr/bin/env python3
"""Bounded coupled numerical transfer; reuses city builder and turbine force audit."""
import argparse
import hashlib
import importlib.util
import json
import math
import re
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
    relaxation=text.split('relaxationFactors',1)[1]
    expected={'p':.3,'U':.6,'k':.6,'epsilon':.6}
    actual={key:float(value) for key,value in re.findall(r'(p|U|k|epsilon) ([0-9.eE+-]+);',relaxation)}
    if actual!=expected:raise ValueError('Retained relaxation factors changed')
    if ('consistent yes;' in text)!=(p['role']!='baseline'):raise ValueError('Unmatched SIMPLE convention')

def build(name,role,common_mesh=False,tight_transport=False,coarsen_city_floor=False):
    started=time.monotonic();cpu=time.process_time()
    if Path(name).name!=name:raise ValueError('One fresh case name required')
    if tight_transport and (role!='reference' or not common_mesh):
        raise ValueError('Tight transport requires --role reference --common-mesh')
    if coarsen_city_floor and (role!='baseline' or not common_mesh):
        raise ValueError('Coarsening the city floor requires --role baseline --common-mesh')
    case=ROOT/'runs'/name
    args=SimpleNamespace(case=name,closure='realizableKE',direction='N',resolution='reference',empty=False,background_m=20.,mesh_repair=False,solver_repair=role!='baseline',consistent_simple=role!='baseline')
    city.build_poc(args, case_root=ROOT/'runs')
    p=json.loads((case/'parameters.json').read_text())
    p.update(joint=True,role=role,scope='Coupled full Niigata geometry and hypothetical full-size SWiFT rotor placement; numerical transfer only, no interaction observations.',max_cells=500000)
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
    fv.write_text(text)
    wake_level=3 if common_mesh else {'reference':3,'baseline':2,'candidate':1}[role]
    support_box=((-20.25,-286.75,6.75),(20.25,-273.25,60.75))
    wake_box=((-40.5,-455.5,0),(40.5,-253,81))
    snappy=case/'system/snappyHexMeshDict';text=snappy.read_text()
    text=text.replace('geometry {\n','geometry {\n sourceSupport {type searchableBox; '+city.box_string(*support_box)+'}\n rotorWake {type searchableBox; '+city.box_string(*wake_box)+'}\n',1)
    text=text.replace('refinementRegions {','refinementRegions { sourceSupport {mode inside; levels ((1E15 4));} rotorWake {mode inside; levels ((1E15 '+str(wake_level)+'));}',1)
    if coarsen_city_floor:
        original='floorRegion {mode inside; levels ((1E15 3));}'
        if text.count(original)!=1:raise ValueError('Expected one original city floor refinement')
        text=text.replace(original,'floorRegion {mode inside; levels ((1E15 2));}')
        p['city_floor_refinement_level']=2
    snappy.write_text(text)
    lock=json.loads((SCENES/'turbine/reference/data-sources.json').read_text())
    for item in lock['files']:
        path=SCENES/'turbine/data/raw'/item['path']
        if path.stat().st_size!=item['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('Turbine input hash mismatch')
    cp,ct=turbine_prepare.rotor_coefficients(SCENES/'turbine/data/raw',8.7)
    area=math.pi*13.5**2
    city.write(case/'system/topoSetDict','actions ({name rotor; type cellSet; action new; source cylinderToCell; p1 (0 -283.375 32.1); p2 (0 -276.625 32.1); radius 13.5;});')
    city.write(case/'constant/fvOptions',f'rotor {{type actuationDiskSource; variant variableScaling; selectionMode cellSet; cellSet rotor; diskArea {area:.12g}; diskDir (0 -1 0); sink true; Cp {cp:.12g}; Ct {ct:.12g}; monitorMethod points; monitorCoeffs {{points ((0 -215 32.1));}} writeToFile true; writeFileStart 0; writeFileEnd 1000000;}}')
    settling_points=[(0.,y,32.1) for y in (-215.,-334.,-361.,-388.,-415.)]
    monitor_text=' '.join('('+' '.join(map(str,point))+')' for point in settling_points)
    control=case/'system/controlDict';control.write_text(control.read_text().replace('points (','points ('+monitor_text+' ',1))
    p['additional_terminal_settling_points_m']=settling_points
    support_config={'base_m':20,'origin_m':[-600,-1200,0],'level':4,'cell_m':1.25,'box':support_box,'cylinder':{'p1':[0,-283.375,32.1],'p2':[0,-276.625,32.1],'radius':13.5}}
    p['common_reference_mesh_recipe']=(common_mesh or role=='reference') and not coarsen_city_floor
    p.update(rotor_D_m=27.,rotor_center_m=[0,-280,32.1],monitor_m=[0,-215,32.1],diskDir=[0,-1,0],Cp=cp,Ct=ct,area_m2=area,density_kg_m3=1.02,source_cell_m=1.25,source_refinement_level=4,wake_level=wake_level,force_support={'configuration_sha256':hashlib.sha256(json.dumps(support_config,sort_keys=True).encode()).hexdigest(),'configuration':support_config},rotor_coefficient_convention='Locked v3 Cp/Ct interpolated at prescribed8.7m/s; native variableScaling uses solved local monitor vector; no joint fitting.',transfer_limitation='Separate30m cheap recipe cannot share the20m source lattice. This adapted20m candidate is not an exact replay or established optimized method.',cell_limit_adjustment={'original_cheap_target':150000,'joint_cap':500000,'reason':'One allowed adjustment: full city geometry accepted20m background plus common1.25m force support. Cheap150k target remains reported unmet if exceeded.'},physical_validation='No joint observations; numerical agreement is not physical or flight validation.')
    queries=json.loads((case/'queries.json').read_text())
    for region,points in (
        ('rotor_wake',[(x/2,y/2,z/2) for z in range(10,131,5) for y in range(-830,-613,5) for x in range(-54,55,5)]),
        ('flight',[(float(x),float(y),float(z)) for z in (20,40,60) for y in range(-450,-219,5) for x in range(-50,51,5)]),
        ('rotor_monitor',[(0.,-215.,32.1)])):
        queries.extend({'region':region,'id':i,'world_xyz_m':q,'solver_xyz_m':q} for i,q in enumerate(points))
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

def mesh(case):
    outer_start=time.monotonic();outer_cpu=time.process_time()
    validate_solver(case,json.loads((case/'parameters.json').read_text()))
    runner.run_stage(case,'mesh','blockMesh > log.blockMesh 2>&1 && snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 && checkMesh -allTopology -allGeometry > log.checkMesh 2>&1 && checkMesh > log.checkMesh.standard 2>&1 && checkMesh -meshQuality > log.checkMesh.quality 2>&1 && topoSet > log.topoSet 2>&1',600)
    runner.accept_mesh(case);p=json.loads((case/'parameters.json').read_text());runner.check_budget(case,p)
    audit=support_audit(case,p)
    if not audit['all_selected_cells_at_fixed_level']:raise ValueError('Source support not uniformly fixed level')
    dump(case/'support-audit.json',audit)
    mesh_cost=json.loads((case/'cost.mesh.json').read_text())
    dump(case/'cost.host-mesh-audit.json',{'elapsed_s':time.monotonic()-outer_start-mesh_cost['elapsed_s'],'cpu_s':time.process_time()-outer_cpu,'scope':'Host dictionary/mesh acceptance and exact support hash; excludes Docker stage already accounted'})

def sample(case):
    commands=''
    if not (case/'log.reconstructPar').exists():commands+='reconstructPar -latestTime > log.reconstructPar 2>&1 && '
    runner.run_stage(case,'sample',commands+'simpleFoam -postProcess -dict system/roiSampling -latestTime > log.roiSampling 2>&1',120)

def samples(case):
    queries=json.loads((case/'queries.json').read_text())
    folder=city.numeric_times(case/'postProcessing/roi')[-1];data={}
    for region in dict.fromkeys(q['region'] for q in queries):
        combined=folder/f'{region}_k_U.xy'
        rows=[r[:3]+r[4:7] for r in city.raw_rows(combined)] if combined.exists() else city.raw_rows(folder/f'{region}_U.xy')
        values={tuple(round(v,7) for v in r[:3]):r[3:6] for r in rows if len(r)>=6 and all(math.isfinite(v) and abs(v)<1e10 for v in r[3:6])}
        for q in (q for q in queries if q['region']==region):
            u=values.get(tuple(round(v,7) for v in q['solver_xyz_m']))
            if u is not None:data[(region,q['id'])]={'xyz':q['world_xyz_m'],'U':u}
    return data,queries,folder.name

def force_audit(case,p,data):
    rows=[]
    for line in (case/'postProcessing/rotor/0/actuationDiskSource.dat').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):rows.append([float(v) for v in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',line)])
    sample_time=float(city.numeric_times(case/'postProcessing/roi')[-1].name)
    row=max((r for r in rows if r[0]<=sample_time),key=lambda r:r[0]);monitor=row[1:4];disk=row[6:9];thrust=row[11]
    expected=.5*p['area_m2']*p['Ct']*sum(v*v for v in monitor)*disk[1]**2/sum(v*v for v in disk)
    relative=abs(thrust-expected)/expected
    valid=('rotor_monitor',0) in data and monitor[1]<0 and disk[1]<0 and thrust>0 and relative<1e-8
    return {'passed':valid,'force_iteration':row[0],'sample_iteration':sample_time,'monitor_is_valid_fluid_sample':('rotor_monitor',0) in data,'monitor_vector_m_s':monitor,'disk_mean_vector_m_s':disk,'kinematic_thrust':thrust,'physical_thrust_N':1.02*thrust,'formula_relative_error':relative,'sink_convention':'Native source sink=true, diskDir(0 -1 0), actual upstream/disk Uy negative; force opposes prescribed flow. Source formula uses Uy squared, not standalone turbine Ux.'}

def wake_settling(case):
    folders=city.numeric_times(case/'postProcessing/monitor')
    if len(folders)<2:return {'passed':False,'reason':'Less than two terminal samples'}
    a,b=folders[-2:]
    def points(folder):
        return {tuple(round(v,7) for v in r[:3]):r[3:6] for r in city.raw_rows(folder/'probes_U.xy') if r[1]<=-280 and all(math.isfinite(v) and abs(v)<1e10 for v in r[3:6])}
    aa,bb=points(a),points(b);keys=aa.keys()&bb.keys()
    error=math.sqrt(sum(math.dist(aa[k],bb[k])**2 for k in keys)/len(keys))/8.7 if keys else None
    return {'iterations':[a.name,b.name],'valid_wake_probes':len(keys),'vector_change_RMSE_over_Uref':error,'passed':len(keys)==4 and error is not None and error<=.002}

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
            'Incomplete totals omit unrecorded intervals. Comparative scientific evaluation '
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

def region_errors(data, ref, keys, uref):
    common=sorted(keys&data.keys())
    errors=sorted(math.dist(data[k]['U'],ref[k]['U']) for k in common)
    rmse=math.sqrt(sum(e*e for e in errors)/len(errors)) if errors else None
    position=.95*(len(errors)-1)
    lo=int(position) if errors else 0
    p95=errors[lo]+(errors[min(lo+1,len(errors)-1)]-errors[lo])*(position-lo) if errors else None
    return {'reference_valid':len(keys),'common_valid':len(common),
            'coverage':len(common)/len(keys) if keys else 0,
            'vector_RMSE_m_s':rmse,'P95_m_s':p95,
            'vector_RMSE_over_Uref':rmse/uref if rmse is not None else None,
            'P95_over_Uref':p95/uref if p95 is not None else None}

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
    d.update(scope=p['scope'],joint=True,role=p['role'],sample_iteration=iteration,force_support=support,force_audit=force,cheap150k_target_met=d['cells']<=150000,cell_limit_adjustment=p['cell_limit_adjustment'],transfer_limitation=p['transfer_limitation'])
    d['sample_matches_final_iteration']=float(iteration)==float(d['solver_final_iteration'])
    d['rotor_wake_settling']=wake_settling(case)
    d['checked_coupled_result']=d['frozen_numerical_checks_pass'] and d['sample_matches_final_iteration'] and d['rotor_wake_settling']['passed'] and force['passed'] and support['all_selected_cells_at_fixed_level']
    d['sampling_counts']={region:{'requested':sum(q['region']==region for q in queries),'valid':sum(k[0]==region for k in data)} for region in dict.fromkeys(q['region'] for q in queries)}
    if reference:
        rp=json.loads((reference/'parameters.json').read_text());rd=city.diagnostics(reference);ref,ref_queries,ref_iteration=samples(reference);rs=support_audit(reference,rp);rf=force_audit(reference,rp,ref)
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
                 and wake_settling(reference)['passed'] and rf['passed'] and rs['all_selected_cells_at_fixed_level'])
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
            keys={k for k in ref if k[0]==base and (region!='building_wake' or ref[k]['xyz'][1]<=-100 and ref[k]['xyz'][2]<=40)}
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
    parser.add_argument('--common-mesh',action='store_true',help='Use identical reference wake mesh to enforce exact support; solver/stopping comparison only')
    parser.add_argument('stage',choices=['build','mesh','solve','sample','evaluate'])
    parser.add_argument('--case',required=True)
    parser.add_argument('--role',choices=['reference','baseline','candidate'])
    parser.add_argument('--reference')
    parser.add_argument('--tight-transport',action='store_true',
                        help='Build a common-mesh reference with inner transport absolute tolerance 1e-10')
    parser.add_argument('--coarsen-city-floor',action='store_true',
                        help='Build a baseline with city floor volume refinement level 2 and the reference wake recipe')
    args=parser.parse_args()
    if args.tight_transport and (args.stage!='build' or args.role!='reference' or not args.common_mesh):
        parser.error('--tight-transport requires build --role reference --common-mesh')
    if args.coarsen_city_floor and (args.stage!='build' or args.role!='baseline' or not args.common_mesh):
        parser.error('--coarsen-city-floor requires build --role baseline --common-mesh')
    case=ROOT/'runs'/args.case
    if args.stage=='build':build(args.case,args.role,args.common_mesh,args.tight_transport,args.coarsen_city_floor)
    elif args.stage=='mesh':mesh(case)
    elif args.stage=='solve':
        validate_solver(case,json.loads((case/'parameters.json').read_text()))
        gate_start=time.monotonic();gate_cpu=time.process_time()
        p=json.loads((case/'parameters.json').read_text())
        if p['role']!='reference':
            reference=ROOT/'runs/poc-N-reference-fixed'
            rp=json.loads((reference/'parameters.json').read_text())
            a=support_audit(case,p);b=support_audit(reference,rp)
            if not all(a[key]==b[key] for key in ('configuration_sha256','selected_cells','spatial_cell_vertices_sha256')):raise ValueError('Exact force support gate failed; solve prohibited')
        dump(case/'cost.host-solve-audit.json',{'elapsed_s':time.monotonic()-gate_start,'cpu_s':time.process_time()-gate_cpu,'scope':'Mandatory exact support comparison to reference before solve; original reference requires no cross-recipe gate'})
        runner.accept_mesh(case);runner.check_budget(case,json.loads((case/'parameters.json').read_text()));runner.run_stage(case,'solve','decomposePar -force > log.decomposePar 2>&1 && mpirun --allow-run-as-root -np 3 simpleFoam -parallel > log.simpleFoam 2>&1 && reconstructPar -latestTime > log.reconstructPar 2>&1',600)
    elif args.stage=='sample':sample(case)
    else:evaluate(case,ROOT/'runs'/args.reference if args.reference else None)
if __name__=='__main__':main()
