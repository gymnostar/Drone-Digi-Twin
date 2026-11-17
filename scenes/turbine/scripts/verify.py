"""Reproduce undisturbed-inflow and fixed-observation-population grid checks."""
import argparse
import json
import re
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
NUMBER=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def main(names,empty,output='verification.json'):
    records=[]
    arrays=[]
    for name in names:
        r=json.loads((ROOT/f'results/{name}/metrics.json').read_text())
        a=np.load(ROOT/f'results/{name}/comparison.npz')
        records.append(dict(run=name,cell_count=r['numerical_diagnostics']['cells'],
                            physical_accuracy=r['primary_ensemble_mean_accuracy'],
                            wall_seconds=r['stage_costs']['total_elapsed_wall_seconds'],
                            cpu_hours=r['stage_costs']['total_cpu_hours']))
        arrays.append(a)
    for a in arrays[1:]:
        for key in ['y_m','z_m','distances_D','common_mask']:
            assert np.array_equal(a[key],arrays[0][key]),key
        assert np.allclose(a['observed_period_mean'],arrays[0]['observed_period_mean'],equal_nan=True)
    finest=arrays[-1]['CFD_matched_period_mean'].mean(0)
    common=arrays[0]['common_mask']
    for record,a in zip(records,arrays):
        model=a['CFD_matched_period_mean'].mean(0)
        mse_by_plane=[float(np.mean((model[i][mask]-finest[i][mask])**2)) for i,mask in enumerate(common)]
        record['numerical_difference_vs_finest_mean_bins_RMSE_m_s']=float(np.sqrt(np.mean(mse_by_plane)))
    report=dict(same_observation_reference_coordinates_bins_and_masks=True,
                grid_check_scope='Three fixed-physics meshes, matched observed bins; finest is a numerical comparison, not physical truth.',
                grid_sequence=records,
                limitation='Nonuniform refinement ratios, changing cylinder-cell support approximation and no asymptotic-order proof; no automatic GCI claim.')
    report['force_support_discretization_warning'] = (
        'The middle grid places cell centers on the cylinder end boundaries, so binary center selection represents only about 71% of nominal cylinder volume. '
        'Thrust normalization still passes, but source-support aliasing confounds a clean three-grid convergence estimate. '
        'Do not call this study mesh independent; coarse/fine differences remain numerical sensitivity evidence only.')
    scheme_path=ROOT/'results/scheme-linearUpwind/metrics.json'
    if scheme_path.exists():
        scheme=json.loads(scheme_path.read_text())
        report['fixed_mesh_velocity_scheme_sensitivity']=dict(
            run='scheme-linearUpwind',cell_count=scheme['numerical_diagnostics']['cells'],
            velocity_scheme='bounded Gauss linearUpwind grad(U)',
            physical_accuracy=scheme['primary_ensemble_mean_accuracy'],
            wall_seconds=scheme['stage_costs']['total_elapsed_wall_seconds'],
            interpretation='Separate same-grid, same-physics scheme comparison; this did not improve measured ensemble mean accuracy.')
    c=ROOT/'runs'/empty
    p=json.loads((c/'parameters.json').read_text())
    a=np.loadtxt(next((c/'postProcessing/sampleGrid').glob('*/*.xy')))
    reference=p['ustar_m_s']/.41*np.log((a[:,2]+p['z0_m'])/p['z0_m'])
    mask=(np.abs(a[:,1])<27)&(a[:,2]>=13.5)&(a[:,2]<=67.5)
    e=a[mask,4]-reference[mask]
    probes=np.array([[float(x) for x in re.findall(NUMBER,line)]
        for line in (c/'postProcessing/probes/0/U').read_text().splitlines()
        if line.strip() and not line.startswith('#')])
    empty_result=dict(run=empty,scope='Undisturbed ABL Ux versus prescribed inlet log profile; |y|<27m,z13.5..67.5m',
        n_sampled_points=int(mask.sum()),bias_m_s=float(e.mean()),RMSE_m_s=float(np.sqrt(np.mean(e*e))),
        maximum_absolute_error_m_s=float(np.max(np.abs(e))),
        final_probe_Ux=probes[-1,1::3].tolist(),probe_x_m=[-65,0,54,81,108,135],
        iterations=int(probes[-1,0]),reference_Uhub_m_s=p['Uhub_m_s'],
        tested_closure=p['closure'],
        core_k_mean_m2_s2=float(a[mask,3].mean()),prescribed_k_m2_s2=p['k_equilibrium_m2_s2'],
        core_k_RMSE_m2_s2=float(np.sqrt(np.mean((a[mask,3]-p['k_equilibrium_m2_s2'])**2))),
        core_k_mean_relative_drift=float((a[mask,3].mean()-p['k_equilibrium_m2_s2'])/p['k_equilibrium_m2_s2']),
        limitations=['cellPoint/grid sampling alters values in coarse first ground-adjacent cells; z<13.5m is excluded only from this explicitly defined transport diagnostic, not from physical wake errors.',
                     'Top slip and chosen closure/default coefficients do not exactly preserve analytical equilibrium everywhere.'])
    (ROOT/f'results/{empty}').mkdir(exist_ok=True)
    (ROOT/f'results/{empty}/metrics.json').write_text(json.dumps(empty_result,indent=2)+'\n')
    report['empty_ABL']=empty_result
    (ROOT/'results'/output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',default=['baseline-v2','verification-mid','verification-fine'])
    p.add_argument('--empty',default='empty-v2')
    p.add_argument('--output',default='verification.json')
    args=p.parse_args()
    main(args.runs,args.empty,args.output)
