#!/usr/bin/env python3
"""Compare the lookup table's trilinear interpolation with CFD sampled between grid nodes.

Draws fixed-seed random world positions, samples the archived CFD case there with
OpenFOAM cellPoint interpolation (the same scheme the table export uses), evaluates
WindField.interpolate at the same positions and writes regional error statistics.
This is descriptive evidence; the frozen grid-node targets are reported as context only.
"""
import argparse
import hashlib
import json
import math
import platform
import random
import time
from pathlib import Path

from lookup import HEIGHTS, ROOT, WindField, read_json, sample_points

# (low, high) per axis: east, north, height above ground, in metres.
TABLE_BOX = ((-200., 200.), (-450., 200.), (5., 80.))
FLIGHT_BOX = ((-50., 50.), (-450., -220.), (20., 60.))
DEFAULT_SEED = 20261124
ROUNDING_DECIMALS = 3  # millimetre positions keep the OpenFOAM point list readable


def random_points(box, count, rng):
    return [tuple(round(rng.uniform(low, high), ROUNDING_DECIMALS) for low, high in box)
            for _ in range(count)]


def draw_positions(count, seed=DEFAULT_SEED):
    """Return (table_set, flight_set): uniform positions in the table box and the flight volume."""
    rng = random.Random(seed)
    return random_points(TABLE_BOX, count, rng), random_points(FLIGHT_BOX, count, rng)


def in_box(point, box):
    return all(low <= value <= high for value, (low, high) in zip(point, box))


def in_rotor_wake(point, wind_from_deg):
    """Same wind-aligned wake bounds as evaluate_lookup.regions, for continuous heights."""
    x, y, z = point
    angle = math.radians(wind_from_deg)
    c, s = math.cos(angle), math.sin(angle)
    across = c * x - s * (y + 280)
    downwind = -(s * x + c * (y + 280))
    return abs(across) <= 27 and 27 <= downwind <= 135 and 5 <= z <= 65


def in_city(point):
    x, y, z = point
    return abs(x) <= 200 and abs(y) <= 200 and z < 20


def height_bands():
    """Bands between consecutive table levels; the top band includes 80 m."""
    last = len(HEIGHTS) - 2
    return [(f"height_{low:g}_to_{high:g}m", low, high, i == last)
            for i, (low, high) in enumerate(zip(HEIGHTS, HEIGHTS[1:]))]


def regions(table_set, flight_set, wind_from_deg):
    """Region name -> positions. Each region is a uniform sample of its volume."""
    result = {
        "all_table_box": list(table_set),
        # Both sets are uniform inside the flight volume, so their union is too.
        "flight_volume": [p for p in table_set if in_box(p, FLIGHT_BOX)] + list(flight_set),
        "rotor_wake": [p for p in table_set if in_rotor_wake(p, wind_from_deg)],
        "city_below_20m": [p for p in table_set if in_city(p)],
    }
    for name, low, high, top in height_bands():
        result[name] = [p for p in table_set if low <= p[2] < high or (top and p[2] == high)]
    return result


def percentile(sorted_values, fraction):
    """Linear interpolation between order statistics, as in study.summarize_vector_errors."""
    position = fraction * (len(sorted_values) - 1)
    low = int(position)
    high = math.ceil(position)
    return sorted_values[low] + (sorted_values[high] - sorted_values[low]) * (position - low)


def compare(cfd, table, positions, uref):
    """Mask counts and vector errors of table values against CFD values at the same positions."""
    errors = []
    cfd_valid = table_masked = table_air = 0
    worst = None
    for p in positions:
        truth, value = cfd[p], table[p]
        if truth is not None:
            cfd_valid += 1
            if value is None:
                table_masked += 1
                continue
            error = math.dist(truth, value)
            errors.append(error)
            if worst is None or error > worst[0]:
                worst = (error, p)
        elif value is not None:
            table_air += 1
    errors.sort()
    report = {"positions": len(positions), "cfd_valid": cfd_valid, "both_valid": len(errors),
              "table_masked_cfd_air": table_masked, "table_air_cfd_masked": table_air,
              "vector_RMSE_m_s": None, "P95_m_s": None, "max_m_s": None, "max_at_world_m": None,
              "vector_RMSE_over_Uref": None, "P95_over_Uref": None, "max_over_Uref": None}
    if errors:
        rmse = math.sqrt(sum(e * e for e in errors) / len(errors))
        p95 = percentile(errors, .95)
        report.update(vector_RMSE_m_s=rmse, P95_m_s=p95, max_m_s=worst[0], max_at_world_m=list(worst[1]),
                      vector_RMSE_over_Uref=rmse / uref, P95_over_Uref=p95 / uref, max_over_Uref=worst[0] / uref)
    return report


def regional_report(field, cfd, table_set, flight_set):
    table = {p: field.interpolate(*p) for p in dict.fromkeys(table_set + flight_set)}
    uref = field.metadata["Uref_m_s"]
    report = {name: compare(cfd, table, positions, uref)
              for name, positions in regions(table_set, flight_set, field.metadata["wind_from_deg"]).items()}
    bands = [(report[name]["vector_RMSE_m_s"], name) for name, *_ in height_bands()
             if report[name]["vector_RMSE_m_s"] is not None]
    return report, (max(bands)[1] if bands else None)


def check_case(case, field):
    parameters = read_json(case / "parameters.json")
    results = ROOT / "results" / case.name
    metrics_file = results / "metrics.json"
    metrics = read_json(metrics_file)
    meta = field.metadata
    if meta["case"] != case.name or metrics["case"] != case.name:
        raise ValueError("Field, case and evaluated result name different cases")
    if parameters != read_json(results / "parameters.json"):
        raise ValueError("Case parameters differ from the evaluated result")
    if meta["iteration"] != metrics["sample_iteration"]:
        raise ValueError("Field iteration differs from the evaluated sample iteration")
    if meta["metrics_sha256"] != hashlib.sha256(metrics_file.read_bytes()).hexdigest():
        raise ValueError("Field was exported from different metrics")
    for key in ("wind_from_deg", "Uref_m_s", "query_sha256", "protocol_sha256"):
        if meta[key] != parameters[key]:
            raise ValueError(f"Field and case differ in {key}")
    return parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True, help="Archived CFD case directory")
    parser.add_argument("--field", type=Path, required=True, help="Released compressed CSV for that case")
    parser.add_argument("--count", type=int, default=20000, help="Positions per set (table box and flight volume)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")
    start = time.perf_counter()
    field = WindField(args.field)
    parameters = check_case(args.case, field)
    iteration = field.metadata["iteration"]
    table_set, flight_set = draw_positions(args.count, args.seed)
    sample_start = time.perf_counter()
    cfd = sample_points(args.case, parameters, iteration, table_set + flight_set)
    sample_elapsed = time.perf_counter() - sample_start
    report, worst_band = regional_report(field, cfd, table_set, flight_set)
    targets = read_json(ROOT.parent / "poc-protocol.json")["targets"]
    result = {
        "description": "Trilinear lookup interpolation (WindField.interpolate) against OpenFOAM cellPoint "
                       "samples of the same CFD solution at random positions between grid nodes. "
                       "Errors measure interpolation only, not CFD accuracy against a reference or real wind.",
        "case": args.case.name, "iteration": iteration, "direction": field.metadata["direction"],
        "wind_from_deg": field.metadata["wind_from_deg"], "Uref_m_s": field.metadata["Uref_m_s"],
        "field_file": args.field.name,
        "field_sha256": hashlib.sha256(args.field.read_bytes()).hexdigest(),
        "seed": args.seed, "positions_per_set": args.count, "position_rounding_m": 10 ** -ROUNDING_DECIMALS,
        "sets": {"table_box_world_m": TABLE_BOX, "flight_volume_world_m": FLIGHT_BOX},
        "region_definitions": {
            "all_table_box": "Table-box set only.",
            "flight_volume": "|east| <= 50, north -450..-220, height 20..60 m; table-box positions inside it plus the flight set.",
            "rotor_wake": "Table-box positions within evaluate_lookup's wind-aligned rotor wake bounds "
                          "(|across| <= 27 m, 27..135 m downwind of the hub at north -280 m, height 5..65 m).",
            "city_below_20m": "Table-box positions with |east|, |north| <= 200 m and height below 20 m.",
            "height_bands": "Table-box positions with low <= height < high between consecutive table levels; the top band includes 80 m.",
        },
        "masks": {"table_masked_cfd_air": "CFD has a finite sample but the table returns None (a contributing node is masked): coverage lost near walls.",
                  "table_air_cfd_masked": "The table returns a velocity where CFD has no finite sample: the table reports air that CFD does not."},
        "regions": report,
        "worst_height_band_by_RMSE": worst_band,
        "grid_node_targets_for_context": {
            "vector_RMSE_over_Uref_max": targets["common_region_vector_RMSE_over_prescribed_Uref_max"],
            "P95_over_Uref_max": targets["common_region_pointwise_vector_error_P95_over_Uref_max"],
            "note": "Frozen targets for grid-node comparisons against a strict reference. Shown for scale only; "
                    "they were not set for off-grid interpolation and are not evaluated as gates here.",
        },
        "platform": platform.platform(), "python": platform.python_version(),
        "elapsed_s": {"cfd_sampling": sample_elapsed, "total": time.perf_counter() - start},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({name: {k: r[k] for k in ("both_valid", "vector_RMSE_over_Uref", "P95_over_Uref")}
                      for name, r in report.items()}, indent=2))


if __name__ == "__main__":
    main()
