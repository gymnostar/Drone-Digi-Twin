#!/usr/bin/env python3
"""Run bounded OpenFOAM stages; preserve logs and measured cgroup costs."""
import argparse
import platform
import socket
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parent
IMAGE='sha256:7df8c49ff2f966deb0d8afba54a5deb19befc93ea4346e974f0f13869572c289'
DOCKER='/usr/local/bin/docker'

def run_stage(case,stage,commands,limit):
    name=f'niigata-{stage}-{int(time.time())}'
    log=case/f'log.{stage}'
    if log.exists():raise SystemExit(f'Refusing to overwrite {log}')
    script='''source /usr/lib/openfoam/openfoam2306/etc/bashrc
set -o pipefail
cd /case
'''+commands+'''
rc=$?
echo CGROUP_CPU_STAT_BEGIN
cat /sys/fs/cgroup/cpu.stat
echo CGROUP_MEMORY_PEAK_BYTES
cat /sys/fs/cgroup/memory.peak
exit "$rc"
'''
    command=[DOCKER,'run','--rm','--name',name,'--cpus','3','--memory','5g','--network','none',
             '-v',f'{case}:/case','--entrypoint','/bin/bash',IMAGE,'-lc',script]
    start=time.monotonic();timed_out=False
    with log.open('w') as stream:
        proc=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        try:rc=proc.wait(timeout=limit)
        except subprocess.TimeoutExpired:
            timed_out=True
            subprocess.run([DOCKER,'kill',name],capture_output=True)
            rc=proc.wait(timeout=30)
    text=log.read_text()
    import re
    cpu=re.findall(r'^usage_usec (\d+)$',text,re.M)
    peak=re.findall(r'CGROUP_MEMORY_PEAK_BYTES\n(\d+)',text)
    result={'stage':stage,'elapsed_s':time.monotonic()-start,'exit_code':rc,'timed_out':timed_out,
            'cpu_s':int(cpu[-1])/1e6 if cpu else None,
            'peak_cgroup_memory_bytes_including_cache':int(peak[-1]) if peak else None,
            'cpu_limit':3,'memory_limit_bytes':5*1024**3,'image':IMAGE,
            'host':socket.gethostname(),'host_platform':platform.platform()}
    (case/f'cost.{stage}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    if rc or timed_out:raise SystemExit(f'{stage} failed; inspect {log}')
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['mesh','solve','sample'])
    parser.add_argument('--case',default='baseline');args=parser.parse_args()
    case=ROOT/'runs'/args.case
    params=json.loads((case/'parameters.json').read_text())
    if args.stage=='mesh':
        commands='blockMesh > log.blockMesh 2>&1 && '
        commands+='snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 && '
        commands+='checkMesh -allTopology -allGeometry > log.checkMesh 2>&1 && checkMesh > log.checkMesh.standard 2>&1 && checkMesh -meshQuality > log.checkMesh.quality 2>&1'
        run_stage(case,'mesh',commands,600)
        accept_mesh(case)
        check_budget(case,params)
    elif args.stage=='solve':
        accept_mesh(case)
        check_budget(case,params)
        run_stage(case,'solve','decomposePar -force > log.decomposePar 2>&1 && mpirun --allow-run-as-root -np 3 simpleFoam -parallel > log.simpleFoam 2>&1 && reconstructPar -latestTime > log.reconstructPar 2>&1',params.get('solve_timeout_s',1200))
        (case/'niigata.foam').touch()
    else:
        if params.get('scenario')!='numerical-abl':raise SystemExit('Explicit sample stage is for numerical-ABL cases only.')
        run_stage(case,'sample','simpleFoam -postProcess -dict system/roiSampling -latestTime > log.roiSampling 2>&1',120)

def check_budget(case,params):
    if 'max_cells' not in params:return
    import re
    text=(case/'log.checkMesh.standard').read_text()
    found=re.search(r'^\s*cells:\s+(\d+)',text,re.M)
    if not found or int(found[1])>params['max_cells']:
        raise SystemExit('Frozen cell limit exceeded or mesh cell count missing.')

def accept_mesh(case):
    check=(case/'log.checkMesh').read_text()
    if 'Mesh OK.' in check:return
    for suffix in ('standard','quality'):
        path=case/f'log.checkMesh.{suffix}'
        if not path.exists() or 'Mesh OK.' not in path.read_text():
            raise SystemExit('Standard/custom numerical mesh quality not accepted.')
    if 'Failed 1 mesh checks.' not in check or 'Concave cells (using face planes) found' not in check:
        raise SystemExit('Additional allGeometry failure requires inspection.')
    (case/'mesh-acceptance.json').write_text(json.dumps({
        'status':'provisional: standard and custom numerical quality pass; allGeometry concavity remains diagnostic',
        'allGeometry_passed':False,'parent_approved':True,
        'reason':'Nonconvex polyhedral cells alone do not invalidate finite-volume geometry when numerical mesh criteria pass; solution/flux/measurement checks remain required.'},indent=2)+'\n')

if __name__=='__main__':main()
