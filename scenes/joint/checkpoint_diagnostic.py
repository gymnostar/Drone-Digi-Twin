#!/usr/bin/env python3
"""Compare unreleased solver checkpoints at the reserved heading after a strict solve timeout."""
import argparse
import filecmp
import json
from pathlib import Path
from types import SimpleNamespace

from evaluate_lookup import regional_report
from study import validate_solver
from lookup import ROOT, WindField, WindSector, read_json, sample_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs=2, type=Path)
    parser.add_argument("--common-iteration", type=int, required=True)
    args = parser.parse_args()
    if args.common_iteration <= 0 or args.cases[0].resolve() == args.cases[1].resolve():
        parser.error("Supply two different cases and a positive shared checkpoint iteration")
    campaign = read_json(ROOT / "lookup/campaign.json")
    angle = campaign["holdout_direction_deg"]
    targets = read_json(ROOT.parent / "poc-protocol.json")["targets"]
    lower, upper = WindField.for_direction("N"), WindField.for_direction("NE")
    predict = WindSector(lower, upper).at_angle(angle)
    prediction = {p: predict(*p) for p in lower.points}
    report = {
        "strict_holdout_qualified": False,
        "scope": "Diagnostic differences against saved iterates from two timed-out strict solves. "
                 "These fields do not establish a converged reference or interpolation accuracy.",
        "wind_from_deg": angle, "common_iteration": args.common_iteration, "iterations_by_case": {}, "cases": [],
        "change_between_checkpoints": {}, "fit_minus_common_checkpoint": {},
    }

    def differences(truth, values):
        scores = regional_report(truth, values, targets)
        for region in scores.values():
            del region["targets_pass"]
        return scores

    # The two checkpoints isolate the coupling algorithm only on the same
    # geometry, stencils, initial conditions and force support.
    for rel in ("constant/polyMesh/points", "constant/polyMesh/faces", "constant/polyMesh/owner",
                "constant/polyMesh/neighbour", "constant/polyMesh/boundary", "constant/fvOptions",
                "system/fvSchemes", "0/U", "0/p", "0/k", "0/epsilon", "0/nut"):
        if not filecmp.cmp(args.cases[0] / rel, args.cases[1] / rel, shallow=False):
            parser.error(f"Checkpoint cases differ in {rel}")
    final_fields = []
    support_hashes = set()
    for case in args.cases:
        p = read_json(case / "parameters.json")
        validate_solver(case, p)
        if any(p[k] != lower.metadata[k] for k in ("protocol_sha256","source_hashes","rotor_D_m","Cp","Ct","closure")):
            parser.error("Checkpoint physics differ from the training fields")
        stage = read_json(case / "cost.solve.json")
        if (not stage["timed_out"] or p["role"] != "reference" or not p["tight_transport"]
                or p["wind_from_deg"] != angle or any(p[k] != v for k, v in campaign["inflow"].items())):
            parser.error(f"{case.name} must be a timed-out strict solve of the reserved inflow")
        support_hashes.add(read_json(case / "support-audit.json")["spatial_cell_vertices_sha256"])
        times = sorted(int(f.name) for f in (case / "processor0").iterdir()
                       if f.name.isdigit() and int(f.name) > 0 and (f / "U").is_file())[-2:]
        if len(times) != 2 or args.common_iteration not in times:
            parser.error(f"{case.name} needs two saved iterates including {args.common_iteration}")
        report["iterations_by_case"][case.name] = times
        earlier, later = (sample_grid(case, p, str(i), reconstruct=True) for i in times)
        if lower.points.keys() != earlier.keys() or earlier.keys() != later.keys():
            parser.error("Checkpoint sampling grid differs from the released fields")
        truth = SimpleNamespace(points=later, metadata=p)
        report["cases"].append(case.name)
        report["change_between_checkpoints"][case.name] = differences(truth, earlier)
        truth = SimpleNamespace(points=(earlier, later)[times.index(args.common_iteration)], metadata=p)
        report["fit_minus_common_checkpoint"][case.name] = differences(truth, prediction)
        final_fields.append(truth)
    if len(support_hashes) != 1:
        parser.error("Checkpoint cases have different rotor support")
    report["between_solvers_at_common_checkpoint"] = differences(final_fields[0], final_fields[1].points)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
