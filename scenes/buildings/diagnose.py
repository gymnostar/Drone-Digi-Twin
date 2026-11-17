#!/usr/bin/env python3
"""Distinguish point interpolation, transported cell fields and inlet assumptions."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parent
IMAGE='sha256:7df8c49ff2f966deb0d8afba54a5deb19befc93ea4346e974f0f13869572c289'

def parse_set(folder,name):
    path=next(folder.glob(name+'_*.xy'))
    names=path.stem[len(name)+1:].split('_')
    result=[]
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        values=list(map(float,line.split()));row={'xyz':values[:3]};index=3
        for field in names:
            size=3 if field=='U' else 1
            row[field]=values[index:index+size] if size==3 else values[index]
            index+=size
        result.append(row)
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',default='inflow-check');args=parser.parse_args()
    case=ROOT/'runs'/args.case
    heights=[.003,.004,.005,.006,.008,.01,.012,.016,.02,.03,.05,.0636,.1,.2,.3,.4,.6,.8,1.0]
    points='\n'.join(f'(0 1.2 {z})' for z in heights)
    measurements=json.loads((case/'measurements.json').read_text())
    measurement_points='\n'.join('('+' '.join(map(str,m['xyz_model_m']))+')' for m in measurements)
    blocks=[]
    for scheme in ('cell','cellPoint'):
        blocks.append(f'''sample_{scheme} {{ type sets; libs (sampling); writeControl timeStep; writeInterval 1;
        interpolationScheme {scheme}; setFormat raw; fields (U k epsilon nut);
        sets {{vertical {{type cloud; axis xyz; points ({points});}}
               measurements {{type cloud; axis xyz; points ({measurement_points});}}}} }}''')
    control='FoamFile { version 2.0; format ascii; class dictionary; object diagnostics;}\nfunctions {\n'+'\n'.join(blocks)+'\n}\n'
    (case/'system/diagnostics').write_text(control)
    command=['/usr/local/bin/docker','run','--rm','--cpus','3','--memory','5g','--network','none',
             '-v',f'{case}:/case','--entrypoint','/bin/bash',IMAGE,'-lc',
             'source /usr/lib/openfoam/openfoam2306/etc/bashrc; cd /case; simpleFoam -postProcess -dict system/diagnostics -latestTime > log.diagnostics 2>&1 && simpleFoam -postProcess -func yPlus -latestTime > log.yPlus 2>&1']
    start=time.monotonic();subprocess.run(command,check=True)
    elapsed=time.monotonic()-start
    data={};measurement_data={}
    for scheme in ('cell','cellPoint'):
        folders=sorted((case/'postProcessing'/f'sample_{scheme}').iterdir(),key=lambda p:float(p.name))
        data[scheme]=parse_set(folders[-1],'vertical')
        measurement_data[scheme]=parse_set(folders[-1],'measurements')
    with (ROOT/'data/raw/AF_caseE.csv').open(encoding='utf-8-sig') as stream:
        profile=[(float(r['z(m)'])/250,float(r['U(m/s)']),float(r['k(m2/s2)'])) for r in csv.DictReader(stream)]
    profile=[(0.,0.,profile[0][2])]+profile
    def interpolate(z,index):
        for a,b in zip(profile,profile[1:]):
            if a[0]<=z<=b[0]:return a[index]+(b[index]-a[index])*(z-a[0])/(b[0]-a[0])
        return profile[-1][index]
    result=[]
    for cell,point in zip(data['cell'],data['cellPoint']):
        z=cell['xyz'][2];u=interpolate(z,1);k=interpolate(z,2)
        row={'z_model_m':z,'z_full_equivalent_m':250*z,'prescribed_U_m_s':u,'prescribed_k_m2_s2':k}
        for scheme,values in [('cell',cell),('cellPoint',point)]:
            speed=math.sqrt(sum(x*x for x in values['U']))
            row[scheme]={'U_mean_vector_magnitude_m_s':speed,'k_m2_s2':values['k'],
                         'epsilon_m2_s3':values['epsilon'],'nut_m2_s':values['nut'],
                         'relative_U_drift':(speed-u)/u,'relative_k_drift':(values['k']-k)/k}
        result.append(row)
    output={'case':args.case,'scope':'postprocessing existing converged/iteration-capped field, no new CFD solve',
            'profile_xy_model_m':[0,1.2],'distance_from_inlet_model_m':0.8,
            'elapsed_s':elapsed,'profiles':result,'yPlus_log':(case/'log.yPlus').read_text().split('yPlus yPlus write:')[-1].strip()}
    params=json.loads((case/'parameters.json').read_text())
    if not params['empty_inflow_verification']:
        error_population={scheme:[] for scheme in data};missing=[]
        for observation in measurements:
            predictions={}
            for scheme,values in measurement_data.items():
                row=min(values,key=lambda r:math.dist(r['xyz'],observation['xyz_model_m']))
                if math.dist(row['xyz'],observation['xyz_model_m'])<1e-7:
                    predictions[scheme]=math.sqrt(sum(v*v for v in row['U']))/params['Uref_m_s']
            if len(predictions)!=2:missing.append(observation['point_id']);continue
            for scheme,value in predictions.items():error_population[scheme].append(value-observation['reference_ratio'])
        output['city_sampler_sensitivity']={'scope':'Same physical field and common finite observation mask; not an accuracy optimization or fitted correction.',
                                            'valid_samples':len(next(iter(error_population.values()))),'total_observations':80,'missing_ids':missing,
                                            'metrics':{scheme:{'bias_ratio':sum(e)/len(e),'MAE_ratio':sum(abs(x) for x in e)/len(e),'RMSE_ratio':math.sqrt(sum(x*x for x in e)/len(e))} for scheme,e in error_population.items()}}
    field_times=sorted([p for p in case.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)],key=lambda p:float(p.name))
    yplus=(field_times[-1]/'yPlus').read_text()
    distributions={}
    for patch in ('ground','buildings'):
        match=re.search(r'\b'+patch+r'\s*\{(.*?)\n\s*\}',yplus,re.S)
        if not match:continue
        values=re.search(r'value\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)',match[1],re.S)
        if not values:continue
        vs=list(map(float,values[1].split()))
        distributions[patch]={'boundary_faces':len(vs),'min':min(vs),'max':max(vs),
                              'faces_yplus_below_30':sum(v<30 for v in vs),
                              'fraction_faces_yplus_below_30':sum(v<30 for v in vs)/len(vs),
                              'denominator':'patch faces, not area weighting; y+<30 is a high-Re wall-function diagnostic'}
    output['yPlus_distributions']=distributions
    (ROOT/'results'/f'{args.case}-diagnostics.json').write_text(json.dumps(output,indent=2)+'\n')
    for row in result:
        if row['z_model_m'] in (.005,.008,.01,.02,.0636):print(json.dumps(row))
    print(output['yPlus_log'])

if __name__=='__main__':main()
