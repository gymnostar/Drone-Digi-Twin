"""Compare sampled CFD with matched measured projected-LOS mean bins, not gusts."""
import argparse
import csv
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def rows(path):
    return [line for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]


def metrics(error):
    a = np.asarray(error)
    return dict(n=int(a.size), bias_m_s=float(a.mean()),
                MAE_m_s=float(np.abs(a).mean()), RMSE_m_s=float(np.sqrt(np.mean(a*a))))


def probe_array(path):
    return np.array([[float(x) for x in re.findall(NUMBER, row)] for row in rows(path)])


def main(name, plot):
    start = time.perf_counter()
    case = ROOT / 'runs' / name
    params = json.loads((case / 'parameters.json').read_text())
    if params.get('poc_study') and not params['poc_study']['observation_comparison_eligible']:
        raise ValueError('This numerical condition does not match the measured SWiFT ensemble; use the numerical common-query comparison.')
    outputs = sorted((case / 'postProcessing/sampleGrid').glob('*/grid_k_U.xy'), key=lambda p: float(p.parent.name))
    source = outputs[-1]
    sample = np.loadtxt(source)
    axes = [np.unique(sample[:, i]) for i in range(3)]
    shape = tuple(map(len, axes))
    assert np.prod(shape) == len(sample)
    assert np.isfinite(sample).all()
    grid = sample.reshape(*shape, 7)
    expected_xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    assert np.allclose(grid[..., :3], expected_xyz, atol=1e-8), 'Unexpected cloud row order'
    ux = RegularGridInterpolator(axes, grid[..., 4], bounds_error=False, fill_value=np.nan)
    # Grid: y=-40.5..40.5 and z=0..81, 3.375m bins; all data retain actual x.
    step, ny, nz = 3.375, 24, 24
    ys = -40.5 + (np.arange(ny)+.5)*step
    zs = (np.arange(nz)+.5)*step
    distance_values = [2, 3, 4, 5]
    periods, means_obs, means_cfd, populations, period_bin_counts = [], [], [], [], []
    for path in sorted((ROOT / 'data/raw/lidar').glob('*.mat')):
        mat = loadmat(path, squeeze_me=True, struct_as_record=False)
        scan, avg = mat['scan'], mat['scan_avg']
        obs, cfd, count_by_d, pop_by_d = [], [], [], []
        for distance in distance_values:
            sums_o, sums_p, counts = [np.zeros(ny*nz) for _ in range(3)]
            n_slots = n_finite = n_spatial = 0
            minimum_projection = float('inf')
            selected = np.flatnonzero((avg.focus_dist_set_D == distance)
                        & (avg.bad_scan_flag == 0) & (avg.low_data_available_flag == 0)
                        & (avg.WD_smooth_flag == 0))
            for i in selected:
                xyz = np.column_stack([getattr(scan.streamwiseCS, a)[i] for a in 'xyz'])
                vlos, projection = np.asarray(scan.vlos[i]), np.asarray(scan.streamwiseCS.vlos_projection[i])
                n_slots += len(vlos)
                finite = (np.isfinite(xyz).all(1) & np.isfinite(vlos) & np.isfinite(projection) & (projection > 0))
                n_finite += int(finite.sum())
                xyz, observed = xyz[finite], (vlos[finite]/projection[finite])
                predicted = ux(xyz)
                ij = np.floor(np.column_stack(((xyz[:,1]+40.5)/step, xyz[:,2]/step))).astype(int)
                valid = (np.isfinite(predicted) & (ij[:,0] >= 0) & (ij[:,0] < ny)
                         & (ij[:,1] >= 0) & (ij[:,1] < nz))
                n_spatial += int(valid.sum())
                if valid.any():
                    minimum_projection = min(minimum_projection,float(projection[finite][valid].min()))
                index = ij[valid, 0]*nz + ij[valid, 1]
                sums_o += np.bincount(index, weights=observed[valid], minlength=ny*nz)
                sums_p += np.bincount(index, weights=predicted[valid], minlength=ny*nz)
                counts += np.bincount(index, minlength=ny*nz)
            enough = counts >= 10
            o, p = [np.full(ny*nz, np.nan) for _ in range(2)]
            o[enough], p[enough] = sums_o[enough]/counts[enough], sums_p[enough]/counts[enough]
            obs.append(o.reshape(ny,nz)); cfd.append(p.reshape(ny,nz)); count_by_d.append(counts.reshape(ny,nz))
            pop_by_d.append(dict(distance_D=distance, selected_scans=len(selected), sample_slots=n_slots,
                                 finite_xyz_vlos_positive_projection=n_finite, inside_comparison_window_and_CFD=n_spatial,
                                 minimum_accepted_projection=minimum_projection,
                                 bins_with_at_least_10_samples=int(enough.sum())))
        periods.append(path.name)
        means_obs.append(obs); means_cfd.append(cfd); populations.append(pop_by_d); period_bin_counts.append(count_by_d)
    observed, predicted = np.array(means_obs), np.array(means_cfd)
    common = np.isfinite(observed).all(0) & np.isfinite(predicted).all(0)
    report = dict(run=name, tested_closure=params['closure'], CFD_sample_time=float(source.parent.name),
        current_population_role='DEVELOPMENT reference used for method selection; prior date-based groups are exploratory, not a confirmatory test of the selected method.',
        validation_scope='Unfitted six-period near-neutral ensemble; binned mean projected-LOS/streamwise proxy, not full-vector validation or a held-out correction test.',
        measurement_operator=dict(reference='released vlos / streamwiseCS.vlos_projection',
            prediction='CFD Ux interpolated at each original measured x,y,z; sampled CFD grid spacing 3.375m',
            limitations=['Projection assumes negligible transverse/vertical LOS contribution; no virtual lidar beam-volume convolution.',
                        '3D interpolation of an exported grid adds sampling error beyond the native mesh.',
                        'Samples share bins but retain actual curved focal-surface x coordinates.',
                        'Binned means are sample-weighted within each period; periods and distances are equally weighted in summary.'],
            minimum_samples_per_period_bin=10, interpolation_extrapolation=False,
            comparison_window_m=dict(y=[-40.5,40.5], z=[0,81]), bin_width_m=step,
            common_mask='Bins finite in all six periods at each downstream distance'),
        periods=periods, populations=populations, by_period_distance=[], by_distance=[])
    report['hub_profile_definition'] = (
        'At the 32.0625m bin center, positive deficit is max(1-U/8.7,0); '
        'center and second-moment width use only the common observed crosswind bins. '
        'These are coverage-limited descriptors, not complete-wake fitted widths.')
    mses, biases, maes = [], [], []
    for pi, period in enumerate(periods):
        for di, distance in enumerate(distance_values):
            err = (predicted[pi,di]-observed[pi,di])[common[di]]
            if not len(err):
                raise ValueError(f'Empty common validation mask at {distance}D')
            vals = metrics(err)
            report['populations'][pi][di]['sample_count_in_common_mask'] = int(np.array(period_bin_counts)[pi,di][common[di]].sum())
            report['by_period_distance'].append(dict(period=period, distance_D=distance, **vals))
            mses.append(vals['RMSE_m_s']**2); biases.append(vals['bias_m_s']); maes.append(vals['MAE_m_s'])
    report['secondary_period_variability_summary'] = dict(n_groups=len(mses),
        bias_m_s=float(np.mean(biases)), MAE_m_s=float(np.mean(maes)), RMSE_m_s=float(np.sqrt(np.mean(mses))),
        RMSE_divided_by_prescribed_hub_speed=float(np.sqrt(np.mean(mses))/params['Uhub_m_s']),
        normalization_m_s=params['Uhub_m_s'], normalization_kind='published ensemble hub speed',
        interpretation='Same prescribed ensemble CFD compared with individual periods; includes between-period variability, not matched-condition model accuracy.')
    mean_o, mean_p = observed.mean(0), predicted.mean(0)
    hub_index = int(np.argmin(np.abs(zs-params['hub_m'])))
    for di, distance in enumerate(distance_values):
        vals = metrics((mean_p[di]-mean_o[di])[common[di]])
        profiles = {}
        for label, data in [('observed',mean_o),('CFD',mean_p)]:
            profile_mask = common[di,:,hub_index]
            yy, speed = ys[profile_mask], data[di,:,hub_index][profile_mask]
            deficits = np.maximum(1-speed/params['Uhub_m_s'], 0)
            center = float(np.sum(yy*deficits)/np.sum(deficits)) if deficits.sum() else None
            width = float(np.sqrt(np.sum(deficits*(yy-center)**2)/deficits.sum())) if deficits.sum() else None
            profiles[label] = dict(height_m=float(zs[hub_index]), valid_bins=int(profile_mask.sum()),
                                  center_m=center, deficit_second_moment_width_m=width,
                                  maximum_deficit_fraction=float(deficits.max()) if len(deficits) else None)
        report['by_distance'].append(dict(distance_D=distance, common_bins=int(common[di].sum()),
            total_window_bins=ny*nz, equal_period_mean_field_errors=vals, hub_profile=profiles))
    plane_metrics = [p['equal_period_mean_field_errors'] for p in report['by_distance']]
    primary_rmse = float(np.sqrt(np.mean([p['RMSE_m_s']**2 for p in plane_metrics])))
    report['primary_ensemble_mean_accuracy'] = dict(
        n_periods=len(periods), n_distance_groups=len(distance_values),
        common_bins_by_distance=[int(v.sum()) for v in common],
        bias_m_s=float(np.mean([p['bias_m_s'] for p in plane_metrics])),
        MAE_m_s=float(np.mean([p['MAE_m_s'] for p in plane_metrics])), RMSE_m_s=primary_rmse,
        RMSE_divided_by_prescribed_hub_speed=primary_rmse/params['Uhub_m_s'],
        weighting='Form equal-period observed and actual-location CFD mean for each bin first; equal-bin means within each distance; equal four-distance means of MSE, then square root.',
        normalization_m_s=params['Uhub_m_s'], normalization_kind='published ensemble hub speed')
    probes = probe_array(case / 'postProcessing/probes/0/U')
    last = probes[-min(50,len(probes)):]
    report['numerical_diagnostics'] = dict(iterations=int(probes[-1,0]),
        convergence_status='SIMPLE residualControl reported convergence' if 'SIMPLE solution converged' in (case/'log.simpleFoam').read_text() else 'iteration cap or incomplete convergence',
        last_50_iteration_Ux_probe_range_m_s=np.ptp(last[:,1::3],axis=0).tolist(),
        final_U_probe_vectors=probes[-1,1:].reshape(-1,3).tolist())
    fluxes = [np.loadtxt(case/f'postProcessing/{key}/0/surfaceFieldValue.dat')[-1,1] for key in ['inletFlux','outletFlux']]
    report['numerical_diagnostics'].update(inlet_flux_m3_s=float(fluxes[0]), outlet_flux_m3_s=float(fluxes[1]),
        absolute_flux_imbalance_over_inlet=float(abs(sum(fluxes))/abs(fluxes[0])))
    force = probe_array(case/'postProcessing/rotor/0/actuationDiskSource.dat')[-1]
    # t, Uref(3), Cp, Ct, Udisk(3), CpStar,CtStar,T,P; incompressible source has rho=1.
    uref, udisk, thrust_kinematic = force[1:4], force[6:9], force[11]
    expected = .5*params['area_m2']*params['Ct']*np.dot(uref,uref)*udisk[0]**2/np.dot(udisk,udisk)
    check = (case/'log.checkMesh').read_text()
    n_cells = int(re.search(r'cells:\s+(\d+)', check)[1])
    selected_cells = int(re.search(r'cellSet rotor now size (\d+)',(case/'log.topoSet').read_text())[1])
    selected_volume = selected_cells*params.get('source_cell_m',params['smallest_wake_cell_m'])**3
    report['numerical_diagnostics'].update(cells=n_cells, selected_rotor_cells=selected_cells,
        selected_rotor_volume_m3=selected_volume,
        geometric_disk_cylinder_volume_m3=params['area_m2']*params['disk_thickness_m'],
        selected_volume_note='Selected cells lie in the uniform innermost refinement; native V_i/totalV weights sum to one.',
        mesh_OK='Mesh OK.' in check, force_monitor_U_m_s=uref.tolist(), disk_volume_mean_U_m_s=udisk.tolist(),
        source_thrust_kinematic_m4_s2=float(thrust_kinematic),
        physical_thrust_N=float(thrust_kinematic*params['density_kg_m3']),
        source_formula_relative_error=float(abs(thrust_kinematic-expected)/expected),
        force_direction_check='upstream monitoring flow positive x; disk and downstream probes have lower Ux, confirming momentum sink',
        force_density_note='simpleFoam native source uses unit density; multiply logged thrust by chosen density for N.')
    report['model_limitations'] = params['model_limitations']
    report['stage_costs'] = {stage:json.loads((case/f'cost.{stage}.json').read_text()) for stage in ['mesh','solve','sample']}
    report['stage_costs']['total_elapsed_wall_seconds'] = sum(report['stage_costs'][stage]['elapsed_wall_seconds'] for stage in ['mesh','solve','sample'])
    report['stage_costs']['total_cpu_hours'] = sum(report['stage_costs'][stage]['cpu_hours'] for stage in ['mesh','solve','sample'])
    report['stage_costs']['total_timing_scope'] = 'Meshing + solve/decompose/reconstruct + native sampling, each including container startup; excludes source acquisition, preparation and this observational assessment.'
    report['stage_costs']['evaluation_wall_seconds'] = time.perf_counter()-start
    result = ROOT/'results'/name
    result.mkdir(parents=True, exist_ok=True)
    if (result/'metrics.json').exists() and not (result/'metrics.v1.period-target.json').exists():
        old = json.loads((result/'metrics.json').read_text())
        if 'equal_period_distance_summary' in old:
            (result/'metrics.v1.period-target.json').write_text(json.dumps(old,indent=2)+'\n')
    (result/'metrics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(result/'comparison.npz', y_m=ys,z_m=zs,distances_D=distance_values,
                        observed_period_mean=observed,CFD_matched_period_mean=predicted,common_mask=common,
                        per_period_bin_sample_count=np.array(period_bin_counts))
    with (result/'matched_mean_bins.csv').open('w') as stream:
        writer=csv.writer(stream)
        writer.writerow(['period','x_D','y_m','z_m','measured_projected_proxy_m_s','CFD_Ux_m_s'])
        for pi,period in enumerate(periods):
            for di,distance in enumerate(distance_values):
                for yi,zi in np.argwhere(common[di]):
                    writer.writerow([period,distance,ys[yi],zs[zi],observed[pi,di,yi,zi],predicted[pi,di,yi,zi]])
    if plot:
        make_plot(result,ys,zs,distance_values,mean_o,mean_p,common,params['Uhub_m_s'])
    report['stage_costs']['complete_assessment_wall_seconds'] = time.perf_counter()-start
    (result/'metrics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(summary=report['primary_ensemble_mean_accuracy'],numerical=report['numerical_diagnostics'],
                          cost=report['stage_costs']['total_elapsed_wall_seconds'],result=str(result)),indent=2))


def make_plot(result,ys,zs,ds,observed,predicted,common,uref):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,4,figsize=(12,8),sharex=True,sharey=True,layout='constrained')
    for i,d in enumerate(ds):
        for j,field in enumerate([observed,predicted,predicted-observed]):
            data=np.where(common[i],field[i],np.nan).T
            im=axes[j,i].pcolormesh(ys,zs,data,shading='nearest',cmap='coolwarm' if j==2 else 'viridis',
                                    vmin=-2 if j==2 else 3,vmax=2 if j==2 else 11)
            axes[j,i].set_aspect('equal')
            if j==0:axes[j,i].set_title(f'{d}D')
            if i==0:axes[j,i].set_ylabel(['Measured proxy\nz [m]','CFD Ux\nz [m]','CFD − measured\nz [m]'][j])
            if j==2:axes[j,i].set_xlabel('y [m]')
            if i==3:fig.colorbar(im,ax=axes[j,:],label='m/s',shrink=.75)
    fig.suptitle(f'SWiFT neutral: {result.name} | matched bins, equal six-period means\n'
                 'Projected-LOS proxy; observed coverage only; not full-vector or gust validation',fontsize=12)
    fig.savefig(result/'validation.png',dpi=140)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name',default='baseline-v2')
    p.add_argument('--plot',action='store_true')
    args=p.parse_args()
    main(args.name,args.plot)
