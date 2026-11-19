import hashlib
import math
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
