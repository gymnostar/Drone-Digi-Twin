from contextlib import redirect_stdout
from io import StringIO
import json
import math
import re
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import study


class CoupledMetricsTest(unittest.TestCase):
    def test_cold_map_cost_includes_both_host_audits(self):
        def stages():
            return {
                'cost': {stage: {'elapsed_s': 1.0}
                         for stage in ('prepare', 'mesh', 'solve', 'sample')},
                'complete_map_generation_elapsed_s': 4.0,
            }

        with TemporaryDirectory() as folder:
            case = Path(folder)
            (case / 'cost.host-mesh-audit.json').write_text(json.dumps({'elapsed_s': 2.0}))
            incomplete = stages()
            study.generation_timing(case, incomplete)
            self.assertFalse(incomplete['generation_timing_complete'])
            self.assertEqual(incomplete['complete_map_generation_elapsed_s'], 6.0)

            (case / 'cost.host-solve-audit.json').write_text(json.dumps({'elapsed_s': 3.0}))
            complete = stages()
            study.generation_timing(case, complete)
            self.assertTrue(complete['generation_timing_complete'])
            self.assertEqual(complete['complete_map_generation_elapsed_s'], 9.0)

    def test_gap_coverage_uses_reference_and_requested_grid(self):
        # At 5 m the center is bracketed. A second missing sample off the
        # requested grid would incorrectly mark the valid edge point as a gap.
        positions = [(-5, 0, 5), (0, 0, 5), (5, 0, 5), (10, 0, 5), (0, 0, 10)]
        queries = [{'region': 'city', 'id': i, 'world_xyz_m': xyz}
                   for i, xyz in enumerate(positions)]
        ref = {('city', i): {'U': [1.0, 0.0, 0.0]} for i in (1, 3, 4)}
        candidate = {key: value for key, value in ref.items() if key != ('city', 1)}

        result = study.descriptive_city_errors(candidate, ref, queries, 8.7)
        gap = result['narrow_gap_proxy']
        self.assertEqual((gap['reference_valid'], gap['common_valid']), (1, 0))
        self.assertIsNone(gap['vector_RMSE_over_Uref'])
        self.assertEqual(result['city_by_height_agl_m']['5']['coverage'], 0.5)
        self.assertEqual(result['city_by_height_agl_m']['10']['coverage'], 1.0)


class DirectionalStudyTest(unittest.TestCase):
    def test_build_presets_and_world_geometry(self):
        def read_point(text, entry):
            match = re.search(r'\b' + entry + r' \(([^()]+)\);', text)
            self.assertIsNotNone(match)
            return tuple(map(float, match[1].split()))

        def approx_point(actual, expected):
            self.assertLess(math.dist(actual, expected), 1e-10)

        with TemporaryDirectory() as tmp:
            scenes = Path(tmp)
            root = scenes / 'joint'
            root.mkdir()
            (root / 'agglomeration').symlink_to(study.ROOT / 'agglomeration')
            city = scenes / 'buildings'
            (city / 'data/raw').mkdir(parents=True)
            stl = city / 'data/raw/BD_caseE.stl'
            stl.write_text('solid fixture\nendsolid fixture\n')
            (city / 'inputs.json').write_text(json.dumps({'files': [{
                'file': stl.name, 'sha256': study.hashlib.sha256(stl.read_bytes()).hexdigest()
            }]}))
            study.dump(scenes / 'turbine/reference/data-sources.json', {'files': []})

            for direction, recipe, recovery in [('N', 'strict', False), ('NE', None, False),
                                                ('NbE', None, False), ('NNE', None, False), ('NNE', 'strict', True),
                                                ('NW', None, False)]:
                with self.subTest(direction=direction, recovery=recovery):
                    name = direction + ('-simple' if recovery else '')
                    args = ['study.py', 'build', '--case', name,
                            '--direction', direction]
                    if recovery:
                        args.append('--standard-simple')
                    if recipe:
                        args.extend(['--recipe', recipe])
                    with (patch.object(study, 'ROOT', root),
                          patch.object(study, 'SCENES', scenes),
                          patch.object(study.city, 'ROOT', city),
                          patch.object(study.turbine_prepare, 'rotor_coefficients',
                                       return_value=(.48146, .76927)),
                          patch.object(sys, 'argv', args), redirect_stdout(StringIO())):
                        study.main()

                    case = root / 'runs' / name
                    p = json.loads((case / 'parameters.json').read_text())
                    queries = json.loads((case / 'queries.json').read_text())
                    first = {q['region']: q for q in queries if q['id'] == 0}
                    theta = study.city.DIRECTIONS[direction]
                    normal = study.city.transform((0, -1, 0), theta, inverse=True)
                    hub = study.HUB_WORLD
                    approx_point(p['rotor_center_m'], hub)
                    approx_point(p['disk_direction_world'], normal)
                    approx_point(p['monitor_m'], [h - 65*n for h,n in zip(hub,normal)])
                    approx_point(first['rotor_monitor']['world_xyz_m'], p['monitor_m'])
                    approx_point(first['flight']['world_xyz_m'], (-50,-450,20))
                    approx_point(first['city']['world_xyz_m'], (-200,-200,5))
                    for q in first.values():
                        approx_point(q['solver_xyz_m'], study.city.transform(q['world_xyz_m'],theta))
                    h_solver = p['rotor_center_solver_m']
                    approx_point(first['rotor_wake']['solver_xyz_m'],
                                 (h_solver[0]-27,h_solver[1]-135,5))
                    topo = (case / 'system/topoSetDict').read_text()
                    p1,p2 = (read_point(topo,name) for name in ('p1','p2'))
                    approx_point([(a+b)/2 for a,b in zip(p1,p2)], h_solver)
                    approx_point([b-a for a,b in zip(p1,p2)], (0,6.75,0))
                    geom = (case / 'system/snappyHexMeshDict').read_text()
                    floor = geom.split('floorRegion {', 1)[1].split('}', 1)[0]
                    if direction == 'N':
                        archived = study.ROOT / 'results/poc-N-reference-tight-transport/parameters.json'
                        old = json.loads(archived.read_text())
                        self.assertEqual(p['query_sha256'], old['query_sha256'])
                        self.assertEqual(p['force_support'], old['force_support'])
                        self.assertIn('type searchableBox;', floor)
                    else:
                        self.assertIn('type searchableRotatedBox;', floor)
                        approx_point(read_point(floor, 'origin'), study.city.transform((-205,-205,0),theta))
                        approx_point(read_point(floor, 'span'), (410,410,10))
                        approx_point(read_point(floor, 'e1'), study.city.transform((1,0,0),theta))
                    self.assertEqual(p.get('standard_simple_recovery', False), recovery)
                    if recovery:
                        self.assertFalse(p['repair_recipe']['consistent_simple'])
                        original = root / 'runs/NNE/system/fvSolution'
                        self.assertEqual((case / 'system/fvSolution').read_text(),
                                         original.read_text().replace('consistent yes;', 'consistent no;'))
                        with self.assertRaisesRegex(ValueError, 'strict NNE'):
                            study.validate_solver(case, dict(p, direction='N'))
                        legacy = dict(p)
                        del legacy['standard_simple_recovery']
                        with self.assertRaisesRegex(ValueError, 'Unmatched SIMPLE'):
                            study.validate_solver(case, legacy)
                    if direction == 'NNE':
                        self.assertEqual(len(p['directional_mesh_repairs_world_m']), 2)
                    if direction == 'NbE':
                        self.assertEqual(theta, 11.25)
                        self.assertEqual(p['directional_mesh_repairs_world_m'], [[[10,17.5,0],[25,30,7.5]]])
                        corner = geom.split('buildingCorner0 {', 1)[1].split('}', 1)[0]
                        self.assertIn('type searchableRotatedBox;', corner)
                        approx_point(read_point(corner, 'origin'), study.city.transform((10,17.5,0),theta))
                        approx_point(read_point(corner, 'span'), (15,12.5,7.5))
                        approx_point(read_point(corner, 'e1'), study.city.transform((1,0,0),theta))
                        self.assertIn('buildingCorner0 {mode inside; levels ((1E15 4));}', geom)
                    if direction in ('N', 'NbE', 'NNE', 'NW'):
                        self.assertEqual(p['role'], 'reference')
                        self.assertTrue(p['tight_transport'])
                    else:
                        self.assertEqual(p['role'], 'baseline')
                        self.assertTrue(p['post_snap_agglomeration'])
                        box = (case / 'system/agglomerationBox').read_text()
                        approx_point(read_point(box, 'min'), (h_solver[0]-26,h_solver[1]-12,1.5))
                        approx_point(read_point(box, 'max'), (h_solver[0]+26,h_solver[1]+12,66))

    def test_conventional_simplec_relaxation_only_for_two_recipes(self):
        with TemporaryDirectory() as tmp:
            scenes = Path(tmp)
            root = scenes / 'joint'
            root.mkdir()
            (root / 'agglomeration').symlink_to(study.ROOT / 'agglomeration')
            city = scenes / 'buildings'
            (city / 'data/raw').mkdir(parents=True)
            stl = city / 'data/raw/BD_caseE.stl'
            stl.write_text('solid fixture\nendsolid fixture\n')
            (city / 'inputs.json').write_text(json.dumps({'files': [{
                'file': stl.name, 'sha256': study.hashlib.sha256(stl.read_bytes()).hexdigest()
            }]}))
            study.dump(scenes / 'turbine/reference/data-sources.json', {'files': []})

            def build(name, direction, *extra):
                args = ['study.py', 'build', '--case', name, '--direction', direction, *extra]
                with (patch.object(study, 'ROOT', root),
                      patch.object(study, 'SCENES', scenes),
                      patch.object(study.city, 'ROOT', city),
                      patch.object(study.turbine_prepare, 'rotor_coefficients',
                                   return_value=(.48146, .76927)),
                      patch.object(sys, 'argv', args), redirect_stdout(StringIO())):
                    study.main()
                case = root / 'runs' / name
                return case, json.loads((case / 'parameters.json').read_text())

            for direction in ('NE', 'NW', 'NNE'):
                with self.subTest(direction=direction):
                    # NW defaults to strict until it has a checked reference.
                    recipe = ('--recipe', 'production') if direction == 'NW' else ()
                    default, dp = build(direction, direction, *recipe)
                    self.assertNotIn('simplec_conventional_relaxation', dp)
                    for factor in ('0.9', '0.7'):
                        case, p = build(f'{direction}-{factor}', direction, *recipe, '--simplec-relaxation', factor)
                        self.assertEqual(p['simplec_conventional_relaxation'], float(factor))
                        self.assertEqual(p['residual_thresholds'], dp['residual_thresholds'])
                        text = (case / 'system/fvSolution').read_text()
                        expected = (default / 'system/fvSolution').read_text().replace(
                            study.RETAINED_RELAXATION,
                            f'relaxationFactors {{equations {{U {factor}; ".*" {factor};}}}}'
                        ).replace('consistent no;', 'consistent yes;')
                        self.assertEqual(text, expected)
                        self.assertIn('consistent yes;', text)
                        self.assertNotIn('fields', text.split('relaxationFactors', 1)[1])
                        if direction == 'NNE':
                            self.assertTrue(study.strict_reference_configuration(case, p))
                        # Any other heading, recipe or relaxation edit stays rejected.
                        with self.assertRaisesRegex(ValueError, 'strict NNE reference or the NE or NW production'):
                            study.validate_solver(case, dict(p, direction='N'))
                        with self.assertRaisesRegex(ValueError, 'Invalid simplec'):
                            study.validate_solver(case, dict(p, simplec_conventional_relaxation=0.8))
                        with self.assertRaisesRegex(ValueError, 'Invalid simplec'):
                            study.validate_solver(case, dict(p, simplec_conventional_relaxation=True))
                        legacy = dict(p)
                        del legacy['simplec_conventional_relaxation']
                        with self.assertRaisesRegex(ValueError, 'Retained relaxation'):
                            study.validate_solver(case, legacy)
                        other = '0.7' if factor == '0.9' else '0.9'
                        with self.assertRaisesRegex(ValueError, 'Conventional SIMPLEC relaxation factors changed'):
                            study.validate_solver(case, dict(p, simplec_conventional_relaxation=float(other)))
            with self.assertRaisesRegex(ValueError, 'strict NNE reference or the NE or NW production'):
                build('N-simplec', 'N', '--simplec-relaxation', '0.9')
            with self.assertRaisesRegex(ValueError, 'strict NNE reference or the NE or NW production'):
                build('NE-strict-simplec', 'NE', '--recipe', 'strict', '--simplec-relaxation', '0.9')
            with self.assertRaisesRegex(ValueError, 'strict NNE reference or the NE or NW production'):
                build('NNE-simple-simplec', 'NNE', '--standard-simple', '--simplec-relaxation', '0.9')
            for direction, recipe in (('NW', 'strict'), ('NbE', 'strict'), ('NbE', 'production'),
                                      ('NNE', 'production'), ('N', 'strict')):
                with self.assertRaisesRegex(ValueError, 'strict NNE reference or the NE or NW production'):
                    build(f'{direction}-{recipe}-simplec', direction, '--recipe', recipe, '--simplec-relaxation', '0.9')

    def test_raw_samples_rotate_velocity_and_preserve_missing_ids(self):
        with TemporaryDirectory() as tmp:
            case = Path(tmp)
            degrees = 45
            study.dump(case / 'parameters.json', {'wind_from_deg': degrees})
            queries = []
            for i, world in enumerate([(160.,-1000.,40.), (0.,-990.,40.), (5.,-995.,40.)]):
                queries.append({'region': 'domain', 'id': i, 'world_xyz_m': world,
                                'solver_xyz_m': study.city.transform(world,degrees)})
            study.dump(case / 'queries.json', queries)
            folder = case / 'postProcessing/roi/123'
            folder.mkdir(parents=True)
            raw = folder / 'domain_k_U.xy'
            # Six-digit raw coordinates, out of query order, with one solid point absent.
            lines = [' '.join(f'{v:.6g}' for v in queries[i]['solver_xyz_m']) + ' 0 0 -8 1'
                     for i in (2,0)]
            raw.write_text('\n'.join(lines))
            data, _, t = study.samples(case)
            self.assertEqual(t, '123')
            self.assertEqual(set(data), {('domain',0), ('domain',2)})
            for row in data.values():
                self.assertLess(math.dist(row['U'],(-math.sqrt(32),-math.sqrt(32),1)),1e-10)
            raw.write_text('\n'.join(lines + [lines[0]]))
            with self.assertRaisesRegex(ValueError, 'unique query ID'):
                study.samples(case)
            raw.write_text('9999 0 0 0 0 -8 1')
            with self.assertRaisesRegex(ValueError, 'unique query ID'):
                study.samples(case)

    def test_probe_and_force_audits_use_yawed_rotor_frame(self):
        with TemporaryDirectory() as tmp:
            case = Path(tmp)
            degrees = 45
            settling = [study.rotor_point((0,y,32.1),degrees) for y in (-215,-334,-361,-388,-415)]
            p = {'additional_terminal_settling_points_m': settling, 'Uref_m_s': 8.7,
                 'probe_change_RMSE_over_Uref_threshold': .002,'diskDir': [0,-1,0],
                 'area_m2': 2.,'Ct': .75,'wind_from_deg': degrees,'density_kg_m3': 1.02}
            for time in ('25','50'):
                folder = case / 'postProcessing/monitor' / time
                folder.mkdir(parents=True)
                points = settling + [(200,-350,10)]
                (folder / 'probes_U.xy').write_text(''.join(
                    ' '.join(f'{v:.6g}' for v in point)+' 0 -8 0\n' for point in points))
            wake = study.wake_settling(case,p)
            self.assertTrue(wake['passed'])
            self.assertEqual(wake['valid_wake_probes'],4)
            force = case / 'postProcessing/rotor/0/actuationDiskSource.dat'
            force.parent.mkdir(parents=True)
            (case / 'postProcessing/roi/50').mkdir(parents=True)
            force.write_text('50 0 -8 0 0 0 0 -4 0 0 0 48\n')
            report = study.force_audit(case,p,{('rotor_monitor',0): {}})
            self.assertTrue(report['passed'])
            self.assertLess(math.dist(report['monitor_vector_m_s'],(-math.sqrt(32),-math.sqrt(32),0)),1e-10)
            force.write_text('49 0 -8 0 0 0 0 -4 0 0 0 48\n')
            self.assertFalse(study.force_audit(case,p,{('rotor_monitor',0): {}})['passed'])

    def test_reference_selection_and_building_wake(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / 'runs/poc-NE-reference-tight-transport'
            p = {'direction': 'NE', 'wind_from_deg': 45., 'post_snap_agglomeration': True}
            study.dump(reference / 'parameters.json', dict(p,role='reference'))
            with patch.object(study, 'ROOT', root):
                self.assertEqual(study.matched_reference(p,None)[0], reference)
            with self.assertRaisesRegex(ValueError, 'same wind direction'):
                study.matched_reference(dict(p,direction='N',wind_from_deg=0.),reference)
        self.assertTrue(study.in_building_wake((-160,0,40),45))
        self.assertFalse(study.in_building_wake((160,-160,40),45))


if __name__ == '__main__':
    unittest.main()
