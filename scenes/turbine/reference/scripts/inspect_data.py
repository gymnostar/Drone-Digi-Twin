"""Audit the released neutral data; no flow simulation or model fitting.

Requires NumPy/SciPy (available in the QGIS Python bundle). All populations and
exclusions are recorded in data-quality.json. Finite lidar samples are availability,
not an independent certification of the publisher's physical quality control.
"""
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.io import loadmat
from scipy.interpolate import griddata

from collect import ROOT, collect

DATA = ROOT / 'data/raw'


def seconds(text):
    h, m, s = text.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)


def read_csv(path):
    with path.open() as stream:
        rows = list(csv.reader(stream))
    t = np.array([seconds(r[0]) for r in rows[1:]])
    values = np.array([[float(v) for v in row[1:]] for row in rows[1:]])
    assert np.all(np.diff(t) > 0), path
    return rows[0][1:], t, values


def csv_audit(path):
    columns, t, values = read_csv(path)
    dt = float(np.median(np.diff(t)))
    return dict(file=str(path.relative_to(DATA)), rows=len(t),
                first_relative_s=float(t[0]), last_relative_s=float(t[-1]),
                median_dt_s=dt,
                expected_slots_between_first_and_last=int(round((t[-1] - t[0]) / dt)) + 1,
                gaps_over_1p5_dt=int(np.sum(np.diff(t) > 1.5 * dt)),
                finite_per_column=dict(zip(columns, np.isfinite(values).sum(0).tolist())))


def concat(values):
    return np.concatenate([np.atleast_1d(v) for v in values])


def stats(values):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    return dict(mean=float(values.mean()), min=float(values.min()), max=float(values.max()))


def wake_preview(report):
    """Availability/shape inspection, not a CFD comparison or paper reproduction.

    Follows the published loader's LOS/projection convention, using explicit
    SciPy linear interpolation instead of the notebook's natural-neighbour option.
    Missing samples/outside-hull points stay NaN. Average scans within each period,
    then average available period means equally at each cell.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    y, z = np.meshgrid(np.arange(-40., 40.1, 2), np.arange(0., 80.1, 2))
    means, coverages = {d: [] for d in [2, 3, 4, 5]}, {d: [] for d in [2, 3, 4, 5]}
    for path in sorted((DATA / 'lidar').glob('*.mat')):
        m = loadmat(path, squeeze_me=True, struct_as_record=False)
        scan, avg = m['scan'], m['scan_avg']
        for distance in means:
            fields = []
            indices = np.flatnonzero((avg.focus_dist_set_D == distance)
                                    & (avg.bad_scan_flag == 0)
                                    & (avg.low_data_available_flag == 0)
                                    & (avg.WD_smooth_flag == 0))
            for i in indices:
                vel = scan.vlos[i]
                projection = scan.streamwiseCS.vlos_projection[i]
                points = np.column_stack((scan.streamwiseCS.y[i], scan.streamwiseCS.z[i]))
                mask = np.isfinite(points).all(1) & np.isfinite(vel) & np.isfinite(projection) & (projection > 0)
                assert mask.sum() >= 3
                fields.append(griddata(points[mask], vel[mask] / projection[mask], (y, z), method='linear'))
            a = np.array(fields)
            count = np.isfinite(a).sum(0)
            mean = np.divide(np.nansum(a, axis=0), count, out=np.full(y.shape, np.nan), where=count > 0)
            means[distance].append(mean)
            coverages[distance].append(count)
    fig, panels = plt.subplots(2, 2, figsize=(8.4, 7.3), sharex=True, sharey=True, layout='constrained')
    axes = panels.ravel()
    fields, period_counts = [], []
    for ax, distance in zip(axes, means):
        a = np.array(means[distance])
        count = np.isfinite(a).sum(0)
        mean = np.divide(np.nansum(a, axis=0), count, out=np.full(y.shape, np.nan), where=count > 0)
        # Display only cells represented in all six periods; retain all coverage in NPZ.
        shown = np.where(count == 6, mean, np.nan)
        fields.append(mean)
        period_counts.append(count)
        im = ax.pcolormesh(y / 27, (z - 32.1) / 27, shown, vmin=3, vmax=11, shading='nearest', cmap='viridis')
        ax.set_title(f'x/D = {distance}')
        ax.set_xlabel('Crosswind y / D')
        ax.set_aspect('equal')
        ax.set_xlim(-1.48, 1.48)
        ax.set_ylim(-1.19, 1.78)
        ax.axhline(0, color='white', alpha=.45, linewidth=.6)
        ax.axvline(0, color='white', alpha=.45, linewidth=.6)
    axes[0].set_ylabel('(Height − hub height) / D')
    axes[2].set_ylabel('(Height − hub height) / D')
    fig.colorbar(im, ax=axes, label='Projected LOS estimate [m/s]', shrink=.8)
    fig.suptitle('SWiFT measurements: equal mean of six periods, only cells covered in all six\n'
                 'Linear interpolation; released streamwise frame; no CFD or wake-centre tracking', fontsize=12)
    output = ROOT / 'reference/wake-data-preview.png'
    fig.savefig(output, dpi=140)
    plt.close(fig)
    np.savez_compressed(DATA / 'wake-preview.npz', y_m=y, z_m=z, distances_D=[2, 3, 4, 5],
                        projected_los_m_s=fields, contributing_periods=period_counts,
                        per_period_scan_count=np.array(list(coverages.values())))
    report['wake_preview'] = dict(method='linear interpolation, vlos/projection, per-period means equally weighted',
                                  displayed_mask='finite in all six period means',
                                  extrapolation=False, CFD=False, wake_tracking=False,
                                  artifact='wake-data-preview.png')
    print('PASS: measured wake fields inspected on a common grid; coverage retained, no extrapolation.')


def main(plot=False):
    collect(offline=True)
    report = dict(status='suitable_for_neutral_mean_wake_benchmark_preparation',
                  cfd_run=False,
                  source_lock_sha256=hashlib.sha256((ROOT / 'reference/data-sources.json').read_bytes()).hexdigest(),
                  versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  csv=[], inflow=[], lidar=[], turbine_model_candidates={})
    paths = sorted((DATA / 'upstream/met_data').glob('*neutral*.csv'))
    for path in paths:
        report['csv'].append(csv_audit(path))
        _, t, values = read_csv(path)
        assert np.isfinite(values).all(), path
        horizontal = np.hypot(values[:, 0], values[:, 1])
        # Fit trends against actual timestamps, so missing samples do not compress time.
        design = np.column_stack((np.ones(len(t)), t - t.mean()))
        fluct = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
        horizontal_fluct = horizontal - design @ np.linalg.lstsq(design, horizontal, rcond=None)[0]
        report['inflow'].append(dict(
            index=path.stem.split('_')[-3], rows=len(t),
            mean_u_m_s=float(values[:, 0].mean()),
            mean_horizontal_speed_m_s=float(horizontal.mean()),
            detrended_sigma_uvw_m_s=fluct.std(0).tolist(),
            streamwise_TI_percent=float(100 * fluct[:, 0].std() / values[:, 0].mean()),
            horizontal_TI_percent=float(100 * horizontal_fluct.std() / horizontal.mean()),
            detrended_k_from_three_components_m2_s2=float(.5 * np.var(fluct, axis=0).sum())))
    for path in sorted((DATA / 'upstream/wt_data/measurements').glob('*neutral*.csv')):
        if 'summary' not in path.name:
            report['csv'].append(csv_audit(path))

    # These files retain original UTC times and geolocated samples. CSV times
    # start at zero: chronological CSV-to-UTC association is only an inference.
    for path in sorted((DATA / 'lidar').glob('*.mat')):
        m = loadmat(path, squeeze_me=True, struct_as_record=False)
        scan, avg = m['scan'], m['scan_avg']
        assert len(scan.vlos) == len(avg.time)
        assert np.all(np.diff(avg.time) > 0)
        valid_scan = ((avg.bad_scan_flag == 0) & (avg.low_data_available_flag == 0)
                      & (avg.WD_smooth_flag == 0))
        by_focus = []
        for focus in np.unique(avg.focus_dist_set_D):
            indices = np.flatnonzero(avg.focus_dist_set_D == focus)
            samples = finite = selected = 0
            for i in indices:
                vel = np.atleast_1d(scan.vlos[i])
                xyz = np.column_stack([getattr(scan.streamwiseCS, a)[i] for a in 'xyz'])
                projection = np.atleast_1d(scan.streamwiseCS.vlos_projection[i])
                assert len(vel) == len(xyz) == len(projection)
                mask = np.isfinite(vel) & np.isfinite(xyz).all(1) & np.isfinite(projection)
                samples += len(vel)
                finite += int(mask.sum())
                if valid_scan[i]:
                    selected += int(mask.sum())
            by_focus.append(dict(focus_D=float(focus), scans=len(indices),
                                 sample_slots=samples, finite_velocity_xyz_projection=finite,
                                 finite_after_released_scan_flags=selected))
        times = concat(scan.time)
        report['lidar'].append(dict(
            file=path.name, scans=len(avg.time),
            first_sample_utc=datetime.fromtimestamp(float(np.nanmin(times)), timezone.utc).isoformat(),
            last_sample_utc=datetime.fromtimestamp(float(np.nanmax(times)), timezone.utc).isoformat(),
            released_flag_rejections=int((~valid_scan).sum()),
            direction_from_deg=stats(avg.meta1.SonicWD_32m),
            yaw_heading_deg=stats(avg.wtga1.YawHeading),
            smoothed_yaw_offset_deg=stats(avg.wtga1.yaw_offset_smooth),
            by_focus=by_focus))
    for version in ['v1', 'v2', 'v3']:
        path = DATA / f'upstream/wt_data/{version}/Simulated_data/RotorPerformances_OF2.csv'
        a = np.loadtxt(path, skiprows=1)
        # Aerodynamic nondimensional coefficients, not electrical Cp nor total tower thrust.
        report['turbine_model_candidates'][version] = dict(
            source=str(path.relative_to(DATA)), Uref_m_s=8.7,
            interpolated_AeroCt=float(np.interp(8.7, a[:, 0], a[:, 5])),
            interpolated_AeroCp=float(np.interp(8.7, a[:, 0], a[:, 4])),
            pitch_deg=float(np.interp(8.7, a[:, 0], a[:, 2])),
            source_kind='OpenFAST simulated performance, not measured thrust')
    report['ensemble_recomputed'] = {
        key: float(np.mean([period[key] for period in report['inflow']]))
        for key in ['mean_u_m_s', 'mean_horizontal_speed_m_s', 'streamwise_TI_percent',
                    'horizontal_TI_percent', 'detrended_k_from_three_components_m2_s2']}
    report['limitations'] = [
        'Six selected periods, not a wind-climate archive or a generalization dataset.',
        'CSV indices follow release order; exact UTC pairing not encoded in CSV headers.',
        'Released vlos is a scalar LOS observation; its projection factor is dimensionless.',
        'Finite sample counts measure availability; no new SNR threshold or gap filling applied.',
        '1D field of view is small; publisher processing excludes <=1.5D and 2.5D.',
        'Steady RANS can compare fixed-frame mean wakes, not wake meandering or gust dynamics.',
        'Documentation still names v1; 2020 paper section 2.3 supports v3 (-0.75 deg pitch). No OpenFOAM forcing configured.',
        'Higher-level raw mast profiles, exact mast survey and instrumental uncertainty remain incomplete.',
    ]
    if plot:
        wake_preview(report)
    out = ROOT / 'reference/data-quality.json'
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print('PASS: CSV uniqueness/cadence/finiteness; all six lidar structures, UTC order, flags and sample populations.')
    print(json.dumps(report['ensemble_recomputed'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wake-preview', action='store_true')
    main(parser.parse_args().wake_preview)
