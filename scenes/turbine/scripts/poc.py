"""Prepare or assess the bounded same-physics SWiFT comparison; never launch CFD."""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import prepare

ROOT = prepare.ROOT
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def prepare_pair(prefix, uhub, roughness, condition_role, empty):
    protocol = json.loads((ROOT.parents[1]/'scene-selection/poc-protocol.json').read_text())
    closure = protocol['common_toy_closure']['model']
    for role, level in [('reference',3),('cheap',1)]:
        prepare.main(name=prefix+'-'+role, level=level, iterations=2000,
            closure=closure, source_level=4, poc_role=role, uhub=uhub,
            roughness=roughness, condition_role=condition_role)
    if empty:
        prepare.main(name=prefix+'-empty',level=1,iterations=2000,empty=True,
            closure=closure,source_level=4,poc_role='cheap',uhub=uhub,
            roughness=roughness,condition_role=condition_role)


def label_list(path):
    text = re.sub(r'/\*.*?\*/|//[^\n]*','',path.read_text(),flags=re.S)
    match = re.search(r'^\s*(\d+)\s*\(\s*(.*?)\s*\)\s*;?\s*$',text,re.M|re.S)
    if not match:
        raise ValueError(f'Expected an ASCII OpenFOAM label list: {path}')
    values = [int(v) for v in match[2].split()]
    if len(values) != int(match[1]):
        raise ValueError(f'Label-list length mismatch: {path}')
    return values


def foam_list_body(path):
    text = re.sub(r'/\*.*?\*/|//[^\n]*', '', path.read_text(), flags=re.S)
    match = re.search(r'^\s*\d+\s*\(\s*(.*?)\s*\)\s*;?\s*\Z', text, re.M|re.S)
    if not match:
        raise ValueError(f'Expected ASCII mesh list: {path}')
    return match[1]


def spatial_support_hash(case, indices):
    # The mesh is unsnapped Cartesian refinement. Exact cell vertex sets prove
    # the selected spatial support independently of numbering and refinement history.
    mesh = case/'constant/polyMesh'
    points = [tuple(float(v) for v in row.split()) for row in re.findall(r'\(([^()]+)\)', foam_list_body(mesh/'points'))]
    faces = [[int(v) for v in row.split()] for row in re.findall(r'\d+\(([^()]*)\)', foam_list_body(mesh/'faces'))]
    owners, neighbours = label_list(mesh/'owner'), label_list(mesh/'neighbour')
    cells = {i:set() for i in indices}
    for face_index, vertices in enumerate(faces):
        for cell in (owners[face_index], neighbours[face_index] if face_index < len(neighbours) else -1):
            if cell in cells:
                cells[cell].update(vertices)
    geometry = sorted(tuple(sorted(tuple(round(v, 9) or 0.0 for v in points[j]) for j in vertices)) for vertices in cells.values())
    if any(len(vertices) < 8 for vertices in geometry):
        raise ValueError('Incomplete selected-cell geometry')
    return hashlib.sha256(json.dumps(geometry, separators=(',',':')).encode()).hexdigest()


def support_audit(case, params):
    indices = label_list(case/'constant/polyMesh/sets/rotor')
    levels = label_list(case/'constant/polyMesh/cellLevel')
    selected_levels = sorted({levels[i] for i in indices})
    wanted = params['source_refinement_level']
    return dict(configuration_sha256=params['force_support']['configuration_sha256'],
        selected_cells=len(indices), spatial_cell_vertices_sha256=spatial_support_hash(case, indices), selected_refinement_levels=selected_levels,
        selected_volume_m3=len(indices)*params['source_cell_m']**3,
        all_selected_cells_at_fixed_level=selected_levels==[wanted],
        interpretation='Same axis-aligned base lattice/refinement box and uniform selected level establish the intended common support; compare selected counts as well as total force.')


def diagnostics(case, params, np):
    log = (case/'log.simpleFoam').read_text()
    lines = [line for line in (case/'postProcessing/probes/0/U').read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    a = np.array([[float(v) for v in re.findall(NUMBER,line)] for line in lines])
    vectors = a[:,1:].reshape(len(a),-1,3)
    window = vectors[-min(50,len(a)):]
    changes = np.sqrt(np.mean(np.sum((window-vectors[-1])**2,axis=2),axis=1))
    flux=[]
    for patch in ['inletFlux','outletFlux']:
        f=np.atleast_2d(np.loadtxt(case/f'postProcessing/{patch}/0/surfaceFieldValue.dat'))
        flux.append(float(f[-1,1]))
    imbalance = abs(sum(flux))/abs(flux[0])
    last_residuals={}
    terminal_log=re.split(r'^Time = ',log,flags=re.M)[-1]
    for field,value in re.findall(r'Solving for (p|Ux|Uy|Uz|k|epsilon), Initial residual = ('+NUMBER+')',terminal_log):
        last_residuals[field]=max(last_residuals.get(field,0.),float(value))
    protocol=json.loads((ROOT.parents[1]/'scene-selection/poc-protocol.json').read_text())
    normalized=float(changes.max()/params['Uhub_m_s'])
    limits=params['poc_study']['residual_control']
    residual_checks={field: all(last_residuals.get(component, float('inf'))<=bound for component in (['Ux','Uy','Uz'] if field=='U' else [field])) for field,bound in limits.items()}
    passes=('SIMPLE solution converged' in log
        and imbalance<=protocol['convergence']['relative_flux_imbalance_max']
        and normalized<=protocol['convergence']['terminal_probe_vector_change_RMSE_over_Uref_max']
        and all(residual_checks.values())
        and np.isfinite(window).all() and (np.abs(window)<1e10).all())
    return dict(solver_iteration=int(a[-1,0]),residual_control_converged='SIMPLE solution converged' in log,
        final_initial_residuals=last_residuals,frozen_residual_checks=residual_checks,relative_flux_imbalance=imbalance,
        terminal_window_iterations=len(window),terminal_window_max_vector_RMSE_change_over_Uref=normalized,
        frozen_convergence_checks_pass=bool(passes))


def pipeline(case):
    costs={key:json.loads((case/f'cost.{key}.json').read_text()) for key in ['prepare','mesh','solve','study']}
    for key in ['mesh','solve','study']:
        if costs[key]['returncode'] or not costs[key].get('accepted_stage',False):
            raise ValueError(f'Unaccepted POC stage: {case.name}/{key}')
    return dict(stages=costs,complete_map_generation_wall_seconds=sum(x['elapsed_wall_seconds'] for x in costs.values()),
        cpu_hours=sum(costs[key]['cpu_hours'] for key in ['mesh','solve','study']),
        peak_cgroup_memory_bytes=max(costs[key]['peak_cgroup_memory_bytes'] for key in ['mesh','solve','study']),
        scope='Source verification/case preparation + mesh + solve/decompose/reconstruct + fixed-query map export, including container starts. Acquisition and observational assessment are separate.',
        memory_scope='Peak cgroup memory including cache; not process RSS')


def query_field(case, region, np):
    paths=sorted((case/'postProcessing/sampleStudy').glob(f'*/{region}_k_U.xy'),key=lambda p:float(p.parent.name))
    a=np.atleast_2d(np.loadtxt(paths[-1]))
    if a.shape[1]!=7:
        raise ValueError('Expected x,y,z,k,Ux,Uy,Uz query columns')
    order=np.lexsort((a[:,2],a[:,1],a[:,0]));a=a[order]
    return a[:,:3],a[:,4:7]


def compare(reference, cheap, output):
    import numpy as np
    rc,cc=[ROOT/'runs'/name for name in [reference,cheap]]
    rp,cp=[json.loads((case/'parameters.json').read_text()) for case in [rc,cc]]
    assert rp['poc_study']['role']=='reference' and cp['poc_study']['role']=='cheap'
    assert rp['rotor_enabled'] and cp['rotor_enabled']
    assert rp['poc_study']['protocol_sha256']==cp['poc_study']['protocol_sha256']
    keys=['rotor_D_m','hub_m','Uhub_m_s','z0_m','ustar_m_s','k_equilibrium_m2_s2','Cp','Ct','area_m2','disk_thickness_m','density_kg_m3','closure','force_variant','velocity_convection_scheme']
    assert all(rp[key]==cp[key] for key in keys),'Physical inputs differ'
    for rel in ['constant/fvOptions','0/U','0/k','0/epsilon','0/nut']:
        assert (rc/rel).read_bytes()==(cc/rel).read_bytes(),rel
    assert rp['poc_study']['fixed_ground_band']==cp['poc_study']['fixed_ground_band']
    support=[support_audit(case,p) for case,p in [(rc,rp),(cc,cp)]]
    support_ok=(all(s['all_selected_cells_at_fixed_level'] for s in support)
        and support[0]['configuration_sha256']==support[1]['configuration_sha256']
        and support[0]['selected_cells']==support[1]['selected_cells']
        and support[0]['spatial_cell_vertices_sha256']==support[1]['spatial_cell_vertices_sha256'])
    costs=[pipeline(case) for case in [rc,cc]]
    for stage in ['mesh','solve','study']:
        for key in ['host_model','architecture','image','container_cpu_limit','container_memory_limit_bytes']:
            assert costs[0]['stages'][stage][key]==costs[1]['stages'][stage][key],key
    protocol=json.loads((ROOT.parents[1]/'scene-selection/poc-protocol.json').read_text())
    assert hashlib.sha256((ROOT.parents[1]/'scene-selection/poc-protocol.json').read_bytes()).hexdigest()==rp['poc_study']['protocol_sha256'],'Protocol changed after preparation'
    targets=protocol['targets'];regions={}
    for region in ['wakeROI','background']:
        xyz,ref=query_field(rc,region,np);cxyz,method=query_field(cc,region,np)
        assert xyz.shape==cxyz.shape and np.allclose(xyz,cxyz,rtol=0,atol=1e-8),'Query coordinates differ'
        rv=np.isfinite(ref).all(1)&(np.abs(ref)<1e10).all(1);common=rv&np.isfinite(method).all(1)&(np.abs(method)<1e10).all(1)
        assert rv.any() and common.any(),'Empty reference/common comparison population'
        e=np.linalg.norm(method[common]-ref[common],axis=1);coverage=float(common.sum()/rv.sum())
        rmse=float(np.sqrt(np.mean(e**2)));p95=float(np.quantile(e,.95));u=rp['Uhub_m_s']
        regions[region]=dict(total_queries=len(xyz),reference_valid_queries=int(rv.sum()),common_valid_queries=int(common.sum()),
            common_fraction_of_reference_valid=coverage,vector_RMSE_m_s=rmse,vector_RMSE_over_Uref=rmse/u,
            pointwise_vector_error_P95_m_s=p95,pointwise_vector_error_P95_over_Uref=p95/u,
            maximum_vector_error_m_s=float(e.max()),
            targets_pass=coverage>=targets['minimum_common_valid_fraction_of_reference_valid_ROI_points']
                and rmse/u<=targets['common_region_vector_RMSE_over_prescribed_Uref_max']
                and p95/u<=targets['common_region_pointwise_vector_error_P95_over_Uref_max'])
    numeric=[diagnostics(case,p,np) for case,p in [(rc,rp),(cc,cp)]]
    observation=dict(status='not_applicable',reason='This is a separately specified numerical condition.')
    if rp['poc_study']['observation_comparison_eligible']:
        paths=[ROOT/'results'/name/'metrics.json' for name in [reference,cheap]]
        if all(p.exists() for p in paths):
            measured=[json.loads(p.read_text())['primary_ensemble_mean_accuracy'] for p in paths]
            ratio=measured[1]['RMSE_m_s']/measured[0]['RMSE_m_s']
            observation=dict(status='assessed_development_reference',reference=measured[0],cheap=measured[1],RMSE_ratio=ratio,
                preservation_target_pass=ratio<=targets['measured_RMSE_ratio_to_matched_reference_max'])
        else:observation=dict(status='not_assessed',reason='Run sampleGrid and the existing observational evaluator for both matched neutral cases; the numerical result alone does not establish the measured-error gate.')
    speedup=costs[0]['complete_map_generation_wall_seconds']/costs[1]['complete_map_generation_wall_seconds']
    report=dict(scope='Same-physics numerical RANS fidelity/cost; not physical truth, flight acceptance or confirmatory field generalization.',
        reference=reference,cheap=cheap,normalization_Uhub_m_s=rp['Uhub_m_s'],
        condition_role=rp['poc_study']['condition_role'],force_support=dict(reference=support[0],cheap=support[1],compatible=support_ok),
        convergence=dict(reference=numeric[0],cheap=numeric[1]),cost=dict(reference=costs[0],cheap=costs[1],speedup=speedup,
            speedup_target_pass=speedup>=targets['complete_map_generation_speedup_vs_same_physics_reference_min']),
        numerical_regions=regions,observational_agreement=observation,
        numerical_target_pass=support_ok and all(d['frozen_convergence_checks_pass'] for d in numeric)
            and all(v['targets_pass'] for v in regions.values()),
        limitations=['One timing per recipe; repeat only if timing noise affects a close decision.',
            'Source consistency uses identical lattice/box configuration, uniform selected level and selected-cell count; no source parameters are fitted.',
            'Cached map-query latency is not measured here and must not be conflated with repeated CFD generation.'])
    out=ROOT/'results'/output
    if out.exists():raise FileExistsError(f'Preserve existing comparison: {out}')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(out),speedup=speedup,numerical_regions=regions,numerical_target_pass=report['numerical_target_pass'],observational_status=observation['status']),indent=2))


def domain_sensitivity(reference, enlarged, output):
    import numpy as np
    cases=[ROOT/'runs'/name for name in (reference,enlarged)]
    params=[json.loads((c/'parameters.json').read_text()) for c in cases]
    support=[support_audit(c,p) for c,p in zip(cases,params)]
    regions={}
    for region in ('wakeROI','background'):
        xyz,u=query_field(cases[0],region,np);other,v=query_field(cases[1],region,np)
        assert xyz.shape==other.shape and np.allclose(xyz,other,rtol=0,atol=1e-8)
        valid=np.isfinite(u).all(1)&np.isfinite(v).all(1)&(np.abs(u)<1e10).all(1)&(np.abs(v)<1e10).all(1)
        e=np.linalg.norm(u[valid]-v[valid],axis=1)
        regions[region]=dict(total=len(xyz),common_valid=int(valid.sum()),vector_RMSE_over_Uref=float(np.sqrt(np.mean(e**2))/params[0]['Uhub_m_s']),P95_over_Uref=float(np.quantile(e,.95)/params[0]['Uhub_m_s']))
    report=dict(scope='Enlarged-domain sensitivity with unchanged local lattice/forcing/queries; no asymptotic independence claim',reference=reference,enlarged=enlarged,
        support=support,exact_support_equal=support[0]['spatial_cell_vertices_sha256']==support[1]['spatial_cell_vertices_sha256'],
        convergence=[diagnostics(c,p,np) for c,p in zip(cases,params)],regions=regions,cost=[pipeline(c) for c in cases])
    out=ROOT/'results'/output
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),regions=regions,exact_support_equal=report['exact_support_equal'])))


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('prepare-pair');p.add_argument('--prefix',default='poc-neutral');p.add_argument('--uhub',type=float,default=8.7);p.add_argument('--roughness',type=float,default=.01158)
    p.add_argument('--condition-role',default='development',choices=['development','numerical-holdout']);p.add_argument('--with-empty',action='store_true')
    c=sub.add_parser('compare');c.add_argument('--reference',default='poc-neutral-reference');c.add_argument('--cheap',default='poc-neutral-cheap');c.add_argument('--output',default='poc-neutral-comparison.json')
    d=sub.add_parser('domain-sensitivity');d.add_argument('--reference',required=True);d.add_argument('--enlarged',required=True);d.add_argument('--output',required=True)
    args=parser.parse_args()
    if args.action=='prepare-pair':prepare_pair(args.prefix,args.uhub,args.roughness,args.condition_role,args.with_empty)
    elif args.action=='compare':compare(args.reference,args.cheap,args.output)
    else:domain_sensitivity(args.reference,args.enlarged,args.output)


if __name__=='__main__':main()
