import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

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


if __name__ == '__main__':
    unittest.main()
