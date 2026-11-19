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

from study import in_building_wake, summarize_vector_errors
from lookup import HEIGHTS, ROOT, WindField, WindSector, read_json


def regions(field):
    angle = math.radians(field.metadata["wind_from_deg"])
    c, s = math.cos(angle), math.sin(angle)
    city, wake, rotor, flight = set(), set(), set(), set()
    for point in field.points:
        x, y, z = point
        if abs(x) <= 200 and abs(y) <= 200 and z in (5, 10, 20, 40, 80):
            city.add(point)
            if in_building_wake(point, field.metadata["wind_from_deg"]):
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
    errors = (math.dist(predictions[p], truth.points[p]) for p in valid if predictions[p] is not None)
    metrics = summarize_vector_errors(errors, len(valid), truth.metadata["Uref_m_s"])
    rmse = metrics["vector_RMSE_over_Uref"]
    return {
        **metrics, "requested": len(positions),
        "predicted_air_at_reference_masked": sum(predictions[p] is not None for p in positions if truth.points[p] is None),
        "targets_pass": (rmse is not None and rmse <= targets["common_region_vector_RMSE_over_prescribed_Uref_max"]
                         and metrics["P95_over_Uref"] <= targets["common_region_pointwise_vector_error_P95_over_Uref_max"]
                         and metrics["coverage"] >= targets["minimum_common_valid_fraction_of_reference_valid_ROI_points"]),
    }


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
    if args.holdout and not args.northeast:
        parser.error("Supply the northeast training endpoint to score a held-out CFD field")
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
    if args.northeast:
        ne = load_field(args.northeast)
        if ne.metadata["wind_from_deg"] != campaign["training_directions_deg"][1]:
            raise ValueError("Unexpected second training direction")
        sector = WindSector(north, ne)
        truth = load_field(args.holdout) if args.holdout else None
        angle = truth.metadata["wind_from_deg"] if truth else campaign["holdout_direction_deg"]
        predict = sector.at_angle(angle)
        report["interpolated_query"] = query_timing(predict, north.points)
        report["timed_estimate_wind_from_deg"] = angle
    if args.holdout:
        holdouts = {campaign["holdout_direction_deg"]: "primary",
                    campaign["secondary_holdout"]["wind_from_deg"]: "secondary"}
        if angle not in holdouts or truth.metadata["role"] != "reference":
            raise ValueError("A strict reference from a reserved direction is required")
        north.check_grid_and_physics(truth)
        linear = sector.at_angle(angle, wind_aligned=False)
        report["holdout_wind_from_deg"] = angle
        report["holdout_role"] = holdouts[angle]
        report["cases"] = [f.metadata["case"] for f in (north, ne, truth)]
        for model, query in (("wind_aligned", predict), ("world_linear", linear), ("north_only_baseline", north.at)):
            values = {p: query(*p) for p in truth.points}
            report[model] = regional_report(truth, values, targets)
        report["all_four_regions_pass"] = gates_pass(report["wind_aligned"])
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
