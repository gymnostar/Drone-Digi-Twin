"""Execute one bounded OpenFOAM stage in the pinned local image; retain costs."""
import argparse
import hashlib
import json
import platform
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGE = 'sha256:7df8c49ff2f966deb0d8afba54a5deb19befc93ea4346e974f0f13869572c289'
DOCKER = '/usr/local/bin/docker'
STAGES = {
    'mesh': ('blockMesh > log.blockMesh 2>&1 && snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 && topoSet > log.topoSet 2>&1 && checkMesh > log.checkMesh 2>&1 && grep -Fq "Mesh OK." log.checkMesh', 600),
    'solve': ('decomposePar > log.decomposePar 2>&1 && mpirun --allow-run-as-root -np 3 simpleFoam -parallel > log.simpleFoam 2>&1 && reconstructPar -latestTime > log.reconstructPar 2>&1', 1200),
    'sample': ('postProcess -func sampleGrid -latestTime > log.sampleGrid 2>&1', 600),
    'study': ('postProcess -func sampleStudy -latestTime > log.sampleStudy 2>&1', 600),
}


def run(name, stage):
    case = ROOT / 'runs' / name
    command, seconds = STAGES[stage]
    params = json.loads((case/'parameters.json').read_text())
    poc = params.get('poc_study')
    if poc is not None and stage in {'mesh','solve'}:
        seconds = poc['max_'+stage+'_seconds']
    if stage == 'study' and not (case/'system/sampleStudy').exists():
        raise ValueError('Study sampling requires a prepared POC common-query dictionary.')
    metrics_path = case / f'cost.{stage}.json'
    if metrics_path.exists():
        raise FileExistsError(f'Stage cost exists: {metrics_path}; preserve retained runs.')
    if stage == 'solve':
        if poc is not None:
            mesh_cost = json.loads((case/'cost.mesh.json').read_text())
            if not mesh_cost.get('accepted_stage',False):
                raise RuntimeError('Solver blocked: POC mesh did not pass its quality and role-size gates.')
        if 'Mesh OK.' not in (case/'log.checkMesh').read_text():
            raise RuntimeError('Solver blocked: mesh did not explicitly pass checkMesh.')
    if stage == 'mesh':
        hashes={str(p.relative_to(case)):hashlib.sha256(p.read_bytes()).hexdigest()
                for folder in ['0','constant','system'] for p in (case/folder).rglob('*') if p.is_file()}
        hashes['parameters.json']=hashlib.sha256((case/'parameters.json').read_bytes()).hexdigest()
        (case/'configuration-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    script = f'''source /usr/lib/openfoam/openfoam2306/etc/bashrc
cd /case
timeout --signal=TERM --kill-after=15s {seconds} /bin/bash -c '{command}'
code=$?
awk '$1=="usage_usec" {{print "CGROUP_CPU_USEC=" $2}}' /sys/fs/cgroup/cpu.stat
printf 'CGROUP_MEMORY_PEAK_BYTES='
cat /sys/fs/cgroup/memory.peak
exit "$code"
'''
    start = time.perf_counter()
    proc = subprocess.run([DOCKER, 'run', '--rm', '--network', 'none', '--cpus', '3',
                           '--memory', '5g', '--entrypoint', '/bin/bash',
                           '-v', f'{case}:/case', IMAGE, '-lc', script],
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter()-start
    (case / f'log.container.{stage}').write_text(proc.stdout)
    result = dict(stage=stage, elapsed_wall_seconds=elapsed, returncode=proc.returncode,
                  container_cpu_limit=3, container_memory_limit_bytes=5*1024**3,
                  image=IMAGE, host=platform.platform(), architecture=platform.machine(),
                  solver='OpenCFD v2306 build fbf00d6b-20230626 linuxARM64GccDPInt32Opt',
                  timing_scope='Host elapsed includes container startup, sourcing and complete stage command.',
                  memory_scope='cgroup v2 memory.peak includes cache; not process RSS',
                  timeout_seconds=seconds)
    result['host_model']=subprocess.check_output(['/usr/sbin/sysctl','-n','hw.model'],text=True).strip()
    result['host_physical_memory_bytes']=int(subprocess.check_output(['/usr/sbin/sysctl','-n','hw.memsize'],text=True))
    result['source_lock_sha256']=hashlib.sha256((ROOT/'reference/data-sources.json').read_bytes()).hexdigest()
    for label, key in [('CPU_USEC','cpu_usage_usec'), ('MEMORY_PEAK_BYTES','peak_cgroup_memory_bytes')]:
        match = re.search(rf'CGROUP_{label}=(\d+)', proc.stdout)
        result[key] = int(match[1]) if match else None
    if result['cpu_usage_usec'] is not None:
        result['cpu_hours'] = result['cpu_usage_usec']/3.6e9
    if poc is not None:
        result['poc_role'] = poc['role']
        result['protocol_sha256'] = poc['protocol_sha256']
        result['gate_failures'] = []
        if stage == 'mesh' and proc.returncode == 0:
            text = (case/'log.checkMesh').read_text()
            match = re.search(r'cells:\s+(\d+)',text)
            result['cells'] = int(match[1]) if match else None
            result['cell_limit'] = poc['max_cells']
            if result['cells'] is None or result['cells'] > poc['max_cells']:
                result['gate_failures'].append('Actual cell count exceeds frozen role limit or is unavailable.')
            if 'Mesh OK.' not in text:
                result['gate_failures'].append('checkMesh did not explicitly report Mesh OK.')
        result['accepted_stage'] = proc.returncode == 0 and not result['gate_failures']
    metrics_path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    if proc.returncode:
        print(proc.stdout[-5000:])
        raise SystemExit(proc.returncode)
    if poc is not None and result['gate_failures']:
        raise SystemExit('; '.join(result['gate_failures']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=STAGES)
    p.add_argument('--name', default='baseline-v2')
    args = p.parse_args()
    run(args.name, args.stage)
