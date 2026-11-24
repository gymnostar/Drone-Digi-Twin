import hashlib
import json
import math
import subprocess
import sys
import unittest
from types import SimpleNamespace

from evaluate_lookup import score
from lookup import ROOT, WindField, WindSector, read_json


class LookupTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.north = WindField.for_direction("N")

    def test_published_fields_match_their_checked_solver_results(self):
        for label, record in read_json(ROOT / "lookup/index.json")["fields"].items():
            with self.subTest(direction=label):
                field = WindField.for_direction(label.lower())
                metrics = ROOT / "results" / record["case"] / "metrics.json"
                self.assertEqual(field.metadata["case"], record["case"])
                self.assertEqual(field.metadata["metrics_sha256"], hashlib.sha256(metrics.read_bytes()).hexdigest())
                self.assertTrue(read_json(metrics)["checked_coupled_result"])

    def test_archived_north_wake_and_mask(self):
        self.assertEqual(self.north.at(0, -350, 20), (0.01480031823, -4.839241765, 0.05637215366))
        self.assertIsNone(self.north.at(-195, 10, 5))
        with self.assertRaises(KeyError):
            self.north.at(0, -349, 20)

    def test_same_short_transect_in_both_released_fields(self):
        # Sample in order across a building mask, including a between-grid point.
        waypoints = [(-195, 5, 5), (-195, 10, 5), (-195, 11, 5), (-195, 25, 5)]
        endpoint_velocities = {
            "N": [(-0.08491, -1.37402, -0.17237), (-1.02957, -3.67691, -0.29512)],
            "NE": [(-0.38890, -0.13720, -0.10408), (-3.01119, -0.85734, -0.28689)],
        }
        for label, (first, last) in endpoint_velocities.items():
            expected = [("valid_air", first), ("masked", None),
                        ("unsampled", None), ("valid_air", last)]
            for xyz, (status, velocity) in zip(waypoints, expected):
                with self.subTest(direction=label, xyz=xyz):
                    command = [sys.executable, str(ROOT / "lookup.py"), "query", label, *map(str, xyz)]
                    result = json.loads(subprocess.check_output(command))
                    self.assertEqual(result["wind_from_deg"], {"N": 0, "NE": 45}[label])
                    self.assertEqual(result["xyz_m"], list(xyz))
                    self.assertEqual(result["status"], status)
                    if velocity is None:
                        self.assertIsNone(result["U_enu_m_s"])
                    else:
                        self.assertLess(math.dist(result["U_enu_m_s"], velocity), 1e-5)

    def test_any_heading_uses_released_or_labelled_nearest_field(self):
        cases = {0: ("N", "released"), 360: ("N", "released"), 45: ("NE", "released"),
                 10: ("N", "nearest_released_unvalidated"), 22.5: ("N", "nearest_released_unvalidated"),
                 30: ("NE", "nearest_released_unvalidated")}
        for degrees, (label, basis) in cases.items():
            with self.subTest(degrees=degrees):
                field, actual = WindField.for_heading(degrees)
                self.assertEqual((field.metadata["direction"], actual), (label, basis))
        for degrees in (-10, 90, 180, 350, math.nan):
            with self.subTest(degrees=degrees), self.assertRaises(ValueError):
                WindField.for_heading(degrees)
        command = [sys.executable, str(ROOT / "lookup.py"), "heading", "30", "0", "-350", "20"]
        result = json.loads(subprocess.check_output(command))
        self.assertEqual((result["field"], result["field_wind_from_deg"], result["heading_basis"]),
                         ("NE", 45, "nearest_released_unvalidated"))
        self.assertEqual(result["U_enu_m_s"], list(WindField.for_direction("NE").at(0, -350, 20)))
        refused = subprocess.run([sys.executable, str(ROOT / "lookup.py"), "heading", "180", "0", "-350", "20"],
                                 capture_output=True, text=True)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("No released or fallback field", refused.stderr)

    def test_sector_reproduces_endpoint_and_requires_both_masks_inside(self):
        points = {(0, 0, 5): (4., -4., 0.), (5, 0, 5): (1., 0., 0.)}
        def field(angle, values):
            value = WindField.__new__(WindField)
            value.metadata = dict(self.north.metadata, wind_from_deg=angle)
            value.points = values
            return value
        end = field(45., {(0, 0, 5): (-4., -4., 2.), (5, 0, 5): None})
        model = WindSector(field(0., points), end)
        self.assertEqual(model.at_angle(22.5, wind_aligned=False)(0, 0, 5), (0., -4., 1.))
        base = (0., -8., 2.)
        end.points[0, 0, 5] = (-math.sqrt(32.), -math.sqrt(32.), 2.)
        points[0, 0, 5] = base
        radians = math.radians(22.5)
        expected = (-8*math.sin(radians), -8*math.cos(radians), 2.)
        self.assertLess(math.dist(model.at(22.5, 0, 0, 5), expected), 1e-12)
        self.assertEqual(model.at(0., 5, 0, 5), (1., 0., 0.))
        self.assertIsNone(model.at(22.5, 5, 0, 5))
        self.assertIsNone(model.at(45., 5, 0, 5))
        for angle in (-1., 46., math.nan):
            with self.assertRaises(ValueError):
                model.at(angle, 0, 0, 5)
        end.metadata["Uref_m_s"] = 10.
        with self.assertRaisesRegex(ValueError, "Uref_m_s"):
            WindSector(field(0., points), end)

    def test_regional_accuracy_and_coverage_use_separate_denominators(self):
        points = {(0, 0, 5): (1., 0., 0.), (5, 0, 5): (2., 0., 0.), (10, 0, 5): None}
        truth = SimpleNamespace(points=points, metadata={"Uref_m_s": 8.7})
        predictions = {(0, 0, 5): (1., 0., 0.), (5, 0, 5): None, (10, 0, 5): (1., 0., 0.)}
        targets = read_json(ROOT.parent / "poc-protocol.json")["targets"]
        report = score(truth, predictions, set(points), targets)
        self.assertEqual(report["vector_RMSE_m_s"], 0.)
        self.assertEqual(report["coverage"], .5)
        self.assertEqual(report["predicted_air_at_reference_masked"], 1)
        self.assertFalse(report["targets_pass"])
