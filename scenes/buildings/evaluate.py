#!/usr/bin/env python3
"""Compare only matched AIJ E ac/N samples; retain coverage and limitations."""
import argparse
import json
import math
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parent

def rows(path):
    return [list(map(float,line.split())) for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]

def velocity_rows(folder,name):
    separate=folder/f'{name}_U.xy'
    if separate.exists():return rows(separate)
    combined=folder/f'{name}_k_U.xy'
    return [row[:3]+row[4:7] for row in rows(combined)]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',default='baseline')
    parser.add_argument('--plot',action='store_true',help='Optional matplotlib comparison plot')
    parser.add_argument('--reference',help='Matched numerical-ABL reference case name')
    args=parser.parse_args()
    case=ROOT/'runs'/args.case
    params=json.loads((case/'parameters.json').read_text())
    if params.get('scenario')=='numerical-abl':
        from poc import evaluate_poc
        evaluate_poc(case,ROOT/'runs'/args.reference if args.reference else None);return
    measurements=json.loads((case/'measurements.json').read_text())
    base=case/'postProcessing/sample'
    times=sorted([p for p in base.iterdir() if p.is_dir()],key=lambda p:float(p.name))
    sample=times[-1]
    if params.get('empty_inflow_verification'):
        import csv
        with (ROOT/'data/raw/AF_caseE.csv').open(encoding='utf-8-sig') as stream:af=list(csv.DictReader(stream))
        approach=velocity_rows(sample,'approach')
        combined=rows(sample/'approach_k_U.xy')
        profile=[]
        for row,raw,ref in zip(approach,combined,af):
            speed=math.sqrt(sum(v*v for v in row[3:6]));target=float(ref['U(m/s)']);target_k=float(ref['k(m2/s2)'])
            profile.append({'z_model_m':row[2],'z_full_equivalent_m':250*row[2],
                            'prescribed_U_m_s':target,'downstream_U_m_s':speed,'relative_U_drift':(speed-target)/target,
                            'prescribed_k_m2_s2':target_k,'downstream_k_m2_s2':raw[3],'relative_k_drift':(raw[3]-target_k)/target_k})
        result={'study':'empty-domain inlet-transport verification; no comparison to city observations',
                'sample_iteration':sample.name,'approach_location_model_xy_m':[0,1.2],
                'distance_from_inlet_model_m':0.8,'profiles':profile,
                'cost':{stage:json.loads((case/f'cost.{stage}.json').read_text()) for stage in ('mesh','solve')},
                'limitations':['Empirical inlet U/k with assumed epsilon and roughness need not be equilibrium profiles.','A refined-profile transport result does not validate city-flow physics.']}
        out=ROOT/'results';out.mkdir(exist_ok=True)
        (out/f'{args.case}-metrics.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2));return
    values=velocity_rows(sample,'measurements')
    comparisons=[]
    for obs in measurements:
        xyz=obs['xyz_model_m']
        found=min(values,key=lambda r:sum((r[i]-xyz[i])**2 for i in range(3)))
        if math.dist(found[:3],xyz)>1e-7:
            comparisons.append({**obs,'U_mean_model_m_s':None,
                                'predicted_mean_vector_magnitude_ratio':None,'error_ratio':None,
                                'missing_reason':'OpenFOAM cloud sampler omitted this point; no valid fluid sample at its exact location'})
            continue
        U=found[3:6]
        valid=all(math.isfinite(v) and abs(v)<1e10 for v in U)
        speed=math.sqrt(sum(v*v for v in U)) if valid else None
        prediction=speed/params['Uref_m_s'] if valid else None
        comparisons.append({**obs,'U_mean_model_m_s':U if valid else None,
                            'predicted_mean_vector_magnitude_ratio':prediction,
                            'error_ratio':prediction-obs['reference_ratio'] if valid else None})
    valid=[r for r in comparisons if r['error_ratio'] is not None]
    errors=[r['error_ratio'] for r in valid]
    n=len(valid)
    metric={'bias_ratio':sum(errors)/n,'MAE_ratio':sum(abs(e) for e in errors)/n,
            'RMSE_ratio':math.sqrt(sum(e*e for e in errors)/n),
            'valid_samples':n,'total_matched_observations':len(measurements)}
    metric.update({key.replace('_ratio','_m_s'):value*params['Uref_m_s'] for key,value in list(metric.items()) if key.endswith('_ratio')})
    approach=velocity_rows(sample,'approach')
    import csv
    with (ROOT/'data/raw/AF_caseE.csv').open(encoding='utf-8-sig') as stream:af=list(csv.DictReader(stream))
    inlet_diagnostic=[]
    for row,ref in zip(approach,af):
        inlet_diagnostic.append({'z_model_m':row[2],'downstream_mean_speed_m_s':math.sqrt(sum(v*v for v in row[3:6])),
                                 'prescribed_inlet_speed_m_s':float(ref['U(m/s)'])})
    history={}
    for t in times[-3:]:
        history[t.name]=velocity_rows(t,'measurements')
    settling=[]
    for (a,ra),(b,rb) in zip(list(history.items()),list(history.items())[1:]):
        deltas=[math.dist(x[3:6],y[3:6]) for x,y in zip(ra,rb) if all(abs(v)<1e10 for v in x[3:6]+y[3:6])]
        settling.append({'iterations':[a,b],'vector_change_RMSE_m_s':math.sqrt(sum(d*d for d in deltas)/len(deltas)),
                         'vector_change_RMSE_over_Uref':math.sqrt(sum(d*d for d in deltas)/len(deltas))/params['Uref_m_s']})
    log=(case/'log.simpleFoam').read_text()
    latest_residual={}
    for field,value in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+)',log):latest_residual[field]=float(value)
    meshlog=(case/'log.checkMesh').read_text()
    cells=re.search(r'^\s*cells:\s+(\d+)',meshlog,re.M)
    iterations=re.findall(r'^Time = ([\d.]+)',log,re.M)
    continuity=re.findall(r'time step continuity errors : sum local = ([\deE+.-]+), global = ([\deE+.-]+), cumulative = ([\deE+.-]+)',log)
    costs={stage:json.loads((case/f'cost.{stage}.json').read_text()) for stage in ('mesh','solve')}
    flux={}
    for name in ('flux','outflux'):
        data=rows(case/'postProcessing'/name/'0/surfaceFieldValue.dat')
        flux[name]=data[-1][-1]
    result={'dataset':'AIJ E v15429018','configuration':'after-construction','wind_direction':'N',
            'fit_performed':False,'sample_iteration':sample.name,'solver_final_iteration':iterations[-1],
            'metrics':metric,'cells':int(cells[1]) if cells else None,
            'last_initial_residuals':latest_residual,'last_continuity_errors':list(map(float,continuity[-1])) if continuity else None,
            'boundary_flux_m3_s':flux,'normalized_flux_imbalance':abs(flux['flux']+flux['outflux'])/abs(flux['flux']),
            'solver_reported_convergence':'SIMPLE solution converged' in log,
            'sample_settling':settling,'approach_profile_diagnostic':inlet_diagnostic,'cost':costs,
            'total_meshing_plus_solve_elapsed_s':sum(c['elapsed_s'] for c in costs.values()),
            'total_cpu_s':sum(c['cpu_s'] for c in costs.values()),
            'limitations':[params['measurement_operator_limitation'],'Single wind direction, no held-out generalization result.',
                           'No measured roughness/epsilon/facility confinement or quantitative experimental uncertainty provided.',
                           'No mesh/domain independence or drone-height/vector validation established.',
                           'AllGeometry check retains concave-cell diagnostic; standard/custom numerical mesh checks pass.']}
    out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/f'{args.case}-metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/f'{args.case}-samples.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,(left,right)=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        obs=[r['reference_ratio'] for r in valid]
        pred=[r['predicted_mean_vector_magnitude_ratio'] for r in valid]
        left.scatter(obs,pred,s=22,color='#1764ac')
        bound=max(obs+pred)*1.05
        left.plot([0,bound],[0,bound],color='black',linewidth=1)
        left.set(xlabel='Measured thermistor mean scalar speed / Uref',
                 ylabel='Magnitude of RANS mean velocity / Uref',xlim=(0,bound),ylim=(0,bound),
                 title=f'N, after construction: {n}/80 matched fluid samples')
        xs=[r['xyz_model_m'][0]*250 for r in valid];ys=[r['xyz_model_m'][1]*250 for r in valid]
        scale=max(abs(e) for e in errors)
        dots=right.scatter(xs,ys,c=errors,cmap='coolwarm',vmin=-scale,vmax=scale,s=42)
        missed=[r for r in comparisons if r['error_ratio'] is None]
        if missed:right.scatter([250*r['xyz_model_m'][0] for r in missed],[250*r['xyz_model_m'][1] for r in missed],marker='x',c='black',label='No valid sample')
        right.set(xlabel='East, full-scale equivalent m',ylabel='North, full-scale equivalent m',title='Predicted minus observed speed ratio',aspect='equal')
        if missed:right.legend(fontsize=8)
        fig.colorbar(dots,ax=right,label='Speed ratio error')
        fig.suptitle('Niigata initial pilot: scalar observation / mean-vector statistic mismatch remains',fontsize=11)
        fig.savefig(out/f'{args.case}-validation.png',dpi=180)
        plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in ('approach_profile_diagnostic','cost')},indent=2))

if __name__=='__main__':main()
