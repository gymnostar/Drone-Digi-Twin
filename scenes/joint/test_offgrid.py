import math
import unittest

import offgrid_check as off
from lookup import WindField


def synthetic_field(velocity, masked=()):
    """A small WindField on the table's corner axes with a given velocity function."""
    field = WindField.__new__(WindField)
    field.metadata = {"Uref_m_s": 8.7, "wind_from_deg": 0}
    axes = ((-200, 0, 200), (-450, -220, 200), (5, 20, 80))
    field.points = {(x, y, z): None if (x, y, z) in masked else velocity(x, y, z)
                    for x in axes[0] for y in axes[1] for z in axes[2]}
    return field


def linear(x, y, z):
    return (0.01 * x + 1, -0.02 * y, 0.001 * z)


class PositionsTest(unittest.TestCase):
    def test_fixed_seed_gives_the_same_positions_inside_their_boxes(self):
        table, flight = off.draw_positions(500, 7)
        self.assertEqual((table, flight), off.draw_positions(500, 7))
        self.assertNotEqual(table, off.draw_positions(500, 8)[0])
        self.assertEqual((len(table), len(flight)), (500, 500))
        self.assertTrue(all(off.in_box(p, off.TABLE_BOX) for p in table))
        self.assertTrue(all(off.in_box(p, off.FLIGHT_BOX) for p in flight))
        self.assertTrue(all(round(v, 3) == v for p in table for v in p))

    def test_regions(self):
        self.assertTrue(off.in_rotor_wake((0, -350, 30), 0))
        self.assertFalse(off.in_rotor_wake((0, -250, 30), 0))
        self.assertFalse(off.in_rotor_wake((0, -350, 70), 0))
        table = [(0, -350, 30), (0, 0, 19.9), (0, 0, 80), (0, -300, 5)]
        flight = [(10, -300, 50)]
        groups = off.regions(table, flight, 0)
        self.assertEqual(groups["flight_volume"], [(0, -350, 30), (10, -300, 50)])
        self.assertEqual(groups["city_below_20m"], [(0, 0, 19.9)])
        self.assertEqual(groups["height_60_to_80m"], [(0, 0, 80)])
        self.assertEqual(groups["height_5_to_10m"], [(0, -300, 5)])
        self.assertEqual(groups["height_20_to_32.1m"], [(0, -350, 30)])
        self.assertEqual(groups["rotor_wake"], [(0, -350, 30)])


class CompareTest(unittest.TestCase):
    def test_error_statistics_and_mask_counts(self):
        positions = [(i, 0, 0) for i in range(6)]
        cfd = {(i, 0, 0): (0., 0., 0.) for i in range(4)}
        cfd.update({(4, 0, 0): (1., 0., 0.), (5, 0, 0): None})
        table = {(0, 0, 0): (3., 4., 0.), (1, 0, 0): (1., 0., 0.), (2, 0, 0): (0., 2., 0.),
                 (3, 0, 0): (0., 0., 0.), (4, 0, 0): None, (5, 0, 0): (1., 1., 1.)}
        report = off.compare(cfd, table, positions, 10.)
        self.assertEqual({k: report[k] for k in ("positions", "cfd_valid", "both_valid",
                                                  "table_masked_cfd_air", "table_air_cfd_masked")},
                         {"positions": 6, "cfd_valid": 5, "both_valid": 4,
                          "table_masked_cfd_air": 1, "table_air_cfd_masked": 1})
        # Errors 5, 1, 2, 0.
        self.assertAlmostEqual(report["vector_RMSE_m_s"], math.sqrt(30 / 4))
        self.assertAlmostEqual(report["P95_m_s"], 2 + 3 * 0.85)
        self.assertEqual((report["max_m_s"], report["max_at_world_m"]), (5., [0, 0, 0]))
        self.assertAlmostEqual(report["max_over_Uref"], 0.5)

    def test_no_common_positions_reports_none(self):
        report = off.compare({(0, 0, 0): None}, {(0, 0, 0): None}, [(0, 0, 0)], 8.7)
        self.assertEqual((report["both_valid"], report["vector_RMSE_m_s"]), (0, None))

    def test_linear_field_interpolates_exactly_and_masks_are_counted(self):
        table_set, flight_set = off.draw_positions(300, 3)
        positions = table_set + flight_set
        field = synthetic_field(linear)
        cfd = {p: linear(*p) for p in positions}
        report, worst = off.regional_report(field, cfd, table_set, flight_set)
        self.assertEqual(report["all_table_box"]["both_valid"], 300)
        self.assertLess(report["all_table_box"]["max_m_s"], 1e-9)
        self.assertIn(worst, [name for name, *_ in off.height_bands()])

        masked = synthetic_field(linear, masked={(0, -220, 20)})
        cfd[table_set[0]] = None
        report, _ = off.regional_report(masked, cfd, table_set, flight_set)
        whole = report["all_table_box"]
        expected_lost = sum(cfd[p] is not None and masked.interpolate(*p) is None for p in table_set)
        self.assertGreater(expected_lost, 0)
        self.assertEqual(whole["table_masked_cfd_air"], expected_lost)
        self.assertEqual(whole["table_air_cfd_masked"], int(masked.interpolate(*table_set[0]) is not None))
        self.assertEqual(whole["cfd_valid"], 299)


if __name__ == "__main__":
    unittest.main()
