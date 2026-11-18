#!/usr/bin/env python3
"""Score a reserved CFD direction on the common lookup grid and time repeated queries."""
import argparse
import json
import math
import platform
import random
import statistics
import time
from pathlib import Path

from lookup import HEIGHTS, ROOT, WindField, WindSector, read_json


def regions(field):
    angle = math.radians(field.metadata["wind_from_deg"])
    c, s = math.cos(angle), math.sin(angle)
    city, wake, rotor, flight = set(), set(), set(), set()
    for point in field.points:
        x, y, z = point
        solver_y = s * x + c * y
        if abs(x) <= 200 and abs(y) <= 200 and z in (5, 10, 20, 40, 80):
            city.add(point)
            if solver_y <= -100 and z <= 40:
                wake.add(point)
        # Wind-aligned offsets from the fixed hub, in the dense CFD wake bounds.
        across = c * x - s * (y + 280)
        downwind = -(s * x + c * (y + 280))
        if abs(across) <= 27 and 27 <= downwind <= 135 and 5 <= z <= 65:
            rotor.add(point)
        if abs(x) <= 50 and -450 <= y <= -220 and z in (20, 40, 60):
            flight.add(point)
    missing = {p for p in city if field.points[p] is None}
    gaps = {
        (x, y, z) for x, y, z in city - missing
        if ((any((x - d, y, z) in missing for d in (5, 10)) and
             any((x + d, y, z) in missing for d in (5, 10))) or
            (any((x, y - d, z) in missing for d in (5, 10)) and
             any((x, y + d, z) in missing for d in (5, 10))))
    }
    city_all_heights = {p for p in field.points if abs(p[0]) <= 200 and abs(p[1]) <= 200}
    return {"city": city, "building_wake": wake, "rotor_wake": rotor, "flight": flight,
            "narrow_gap_proxy": gaps, "full_table": set(field.points),
            **{f"city_{z}m": {p for p in city_all_heights if p[2] == z} for z in HEIGHTS}}


def score(truth, predictions, positions, targets):
    valid = [p for p in positions if truth.points[p] is not None]
    common = [p for p in valid if predictions[p] is not None]
    errors = sorted(math.dist(predictions[p], truth.points[p]) for p in common)
    speed = truth.metadata["Uref_m_s"]
    rmse = math.sqrt(sum(e * e for e in errors) / len(errors)) if errors else None
    p95 = None
    if errors:
        i = .95 * (len(errors) - 1)
        p95 = errors[int(i)] + (i % 1) * (errors[math.ceil(i)] - errors[int(i)])
    coverage = len(common) / len(valid) if valid else 0.
    false_air = sum(predictions[p] is not None for p in positions if truth.points[p] is None)
    passed = (rmse is not None and rmse / speed <= targets["common_region_vector_RMSE_over_prescribed_Uref_max"]
              and p95 / speed <= targets["common_region_pointwise_vector_error_P95_over_Uref_max"]
              and coverage >= targets["minimum_common_valid_fraction_of_reference_valid_ROI_points"])
    return {"requested": len(positions), "reference_valid": len(valid), "common_valid": len(common),
            "predicted_air_at_reference_masked": false_air, "coverage": coverage,
            "vector_RMSE_m_s": rmse, "P95_m_s": p95,
            "vector_RMSE_over_Uref": rmse / speed if rmse is not None else None,
            "P95_over_Uref": p95 / speed if p95 is not None else None, "targets_pass": passed}


def regional_report(truth, values, targets):
    return {region: score(truth, values, points, targets) for region, points in regions(truth).items()}


def gates_pass(report):
    return all(report[r]["targets_pass"] for r in ("city", "building_wake", "rotor_wake", "flight"))


def query_timing(query, positions):
    seed = 20261118
    workload = random.Random(seed).choices(sorted(positions), k=10000)
    valid = sum(query(*p) is not None for p in workload)
    batches = []
    for _ in range(5):
        start = time.perf_counter()
        count = sum(query(*p) is not None for p in workload)
        elapsed = time.perf_counter() - start
        if count != valid:
            raise ValueError("Inconsistent query results across repeats")
        batches.append(elapsed * 1e6 / len(workload))
    return {"calls_per_batch": len(workload), "valid_per_batch": valid, "seed": seed,
            "batch_mean_us_per_call": batches, "median_batch_mean_us_per_call": statistics.median(batches),
            "scope": "Warm in-process single-position calls, including Python loop and valid check. "
                     "One fixed inflow with heading weights prepared before the calls. Uniform fixed-grid requests with replacement, including mask hits. One warmup plus five batches."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field", type=Path)
    parser.add_argument("--compare-reference", type=Path, help="Check the actual table against a same-direction strict field")
    parser.add_argument("--northeast", type=Path)
    parser.add_argument("--holdout", type=Path)
    args = parser.parse_args()
    if bool(args.northeast) != bool(args.holdout):
        parser.error("Supply both training endpoint and reserved CFD field for directional assessment")
    campaign = read_json(ROOT / "lookup/campaign.json")
    targets = read_json((ROOT / "lookup" / campaign["targets_source"]).resolve())["targets"]
    report = {"platform": platform.platform(), "python": platform.python_version(), "load_s": {}}

    def load_field(path):
        start = time.perf_counter()
        field = WindField(path)
        report["load_s"][path.name] = time.perf_counter() - start
        if any(field.metadata[k] != v for k, v in campaign["inflow"].items()):
            raise ValueError("Inflow differs from the directional campaign")
        return field

    north = load_field(args.field)
    if args.compare_reference:
        if args.northeast or args.holdout:
            parser.error("Choose a same-direction or a held-out comparison")
        ref = load_field(args.compare_reference)
        north.check_grid_and_physics(ref)
        if ref.metadata["role"] != "reference" or north.metadata["wind_from_deg"] != ref.metadata["wind_from_deg"]:
            parser.error("Comparison requires a same-direction strict reference")
        report["cases"] = [north.metadata["case"], ref.metadata["case"]]
        report["regions"] = regional_report(ref, north.points, targets)
        report["all_four_regions_pass"] = gates_pass(report["regions"])
        print(json.dumps(report, indent=2, allow_nan=False))
        return
    if north.metadata["wind_from_deg"] != campaign["training_directions_deg"][0]:
        raise ValueError("Unexpected first training direction")
    report["direct_query"] = query_timing(north.at, north.points)
    if args.holdout:
        ne, truth = load_field(args.northeast), load_field(args.holdout)
        angle = truth.metadata["wind_from_deg"]
        if angle != campaign["holdout_direction_deg"] or ne.metadata["wind_from_deg"] != campaign["training_directions_deg"][1]:
            raise ValueError("Unexpected training/holdout split")
        if truth.metadata["role"] != "reference":
            raise ValueError("A strict holdout reference is required")
        sector = WindSector(north, ne)
        north.check_grid_and_physics(truth)
        predict = sector.at_angle(angle)
        linear = sector.at_angle(angle, wind_aligned=False)
        report["interpolated_query"] = query_timing(predict, truth.points)
        report["holdout_wind_from_deg"] = angle
        report["cases"] = [f.metadata["case"] for f in (north, ne, truth)]
        for model, query in (("wind_aligned", predict), ("world_linear", linear), ("north_only_baseline", north.at)):
            values = {p: query(*p) for p in truth.points}
            report[model] = regional_report(truth, values, targets)
        report["all_four_regions_pass"] = gates_pass(report["wind_aligned"])
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
