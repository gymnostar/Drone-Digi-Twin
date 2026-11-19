#!/usr/bin/env python3
"""Export and query sampled mean wind at exact world coordinates (metres, east/north/up)."""
import argparse
import csv
import gzip
import hashlib
import io
from importlib import import_module
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HEIGHTS = (5, 10, 20, 32.1, 40, 60, 80)
HEADER = ["x_east_m", "y_north_m", "z_agl_m", "valid_air", "u_east_m_s", "u_north_m_s", "u_up_m_s"]


def read_json(path):
    return json.loads(path.read_bytes())


def planner_grid():
    """The same local grid for every inflow, including the swept rotor wake."""
    return [(x, y, z) for z in HEIGHTS
            for y in range(-450, 201, 5) for x in range(-200, 201, 5)]


def sample_grid(case, parameters, iteration, *, reconstruct=False):
    study = import_module(".study", __package__) if __package__ else import_module("study")

    xyz = planner_grid()
    degrees = parameters["wind_from_deg"]
    expected = {p: study.city.transform(p, degrees) for p in xyz}
    input_dir = case / "processor0" if reconstruct else case
    if not iteration.isdigit() or not (input_dir / iteration / "U").is_file():
        raise ValueError(f"No saved velocity field for iteration {iteration}")
    # Keep archived meshes and fields read-only. The solver writes only into
    # the temporary sampling case, including its postProcessing directory.
    with tempfile.TemporaryDirectory() as directory:
        scratch = Path(directory).resolve()
        system = scratch / "system"
        system.mkdir()
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            shutil.copyfile(case / "system" / name, system / name)
        points = " ".join(study.foam_point(p) for p in expected.values())
        study.city.write(system / "gridSampling",
                         "functions {grid {type sets; libs (sampling); writeControl timeStep; "
                         "writeInterval 1; interpolationScheme cellPoint; setFormat raw; fields (U); "
                         "sets {world {type cloud; axis xyz; points (" + points + ");}}}}")
        if reconstruct:
            shutil.copyfile(case / "system/decomposeParDict", system / "decomposeParDict")
            prepare = 'ln -s /source/processor[0-9]* . && reconstructPar -time "$1"'
        else:
            prepare = 'ln -s "/source/$1" "$1"'
        container = f"niigata-grid-{scratch.name}"
        command = [study.runner.DOCKER, "run", "--rm", "--name", container, "--cpus", "3", "--memory", "5g",
                   "--network", "none", "-v", f"{case.resolve()}:/source:ro",
                   "-v", f"{scratch}:/case", "--workdir", "/case", "--entrypoint", "/bin/bash",
                   parameters["image"], "-lc",
                   "source /usr/lib/openfoam/openfoam2306/etc/bashrc && "
                   "ln -s /source/constant constant && " + prepare + " && "
                   'timeout 110 simpleFoam -postProcess -dict system/gridSampling -time "$1"', "grid", iteration]
        try:
            run = subprocess.run(command, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            subprocess.run([study.runner.DOCKER, "rm", "-f", container], capture_output=True, timeout=20)
            raise
        if run.returncode:
            raise RuntimeError((run.stdout + run.stderr)[-3000:])
        raw = study.city.raw_rows(scratch / "postProcessing" / "grid" / iteration / "world_U.xy")
        sampled = study.raw_vectors(expected, raw)
    if not sampled:
        raise ValueError("No valid grid samples")
    return {p: study.city.transform(sampled[p], degrees, inverse=True) if p in sampled else None
            for p in xyz}


def export(case, results, destination):
    # The standalone reader needs only Python. Import the CFD tools only for export.
    study = import_module(".study", __package__) if __package__ else import_module("study")

    parameters = read_json(case / "parameters.json")
    metrics_file = results / "metrics.json"
    metrics = read_json(metrics_file)
    if parameters != read_json(results / "parameters.json") or metrics["case"] != case.name:
        raise ValueError("Case and evaluated result differ")
    if not metrics["checked_coupled_result"]:
        raise ValueError("Coupled convergence, force or settling checks failed")
    if parameters["role"] == "reference":
        if not parameters.get("tight_transport") or not study.strict_reference_configuration(case, parameters):
            raise ValueError("A checked strict reference is required")
    else:
        comparison = metrics["comparison"]
        if not all(comparison[key] for key in (
            "reference_checked", "support_exact_match", "all_declared_ROI_targets_pass"
        )):
            raise ValueError("Candidate failed the matched-reference field comparison")
    queries_file = case / "queries.json"
    if hashlib.sha256(queries_file.read_bytes()).hexdigest() != parameters["query_sha256"]:
        raise ValueError("Query coordinates differ from the evaluated case")
    iteration = metrics["sample_iteration"]
    if iteration != metrics["solver_final_iteration"]:
        raise ValueError("Sample time and final evaluated iteration differ")
    field = sample_grid(case, parameters, iteration)
    original, queries, sample_time = study.samples(case)
    if sample_time != iteration:
        raise ValueError("Original sample time differs from the final field")
    overlap = 0
    max_difference = 0.
    for q in queries:
        xyz = tuple(q["world_xyz_m"])
        if xyz not in field:
            continue
        old = original.get((q["region"], q["id"]))
        new = field[xyz]
        if (old is None) != (new is None):
            raise ValueError(f"Lookup and original sampling masks disagree at {xyz}")
        if new is not None:
            max_difference = max(max_difference, math.dist(old["U"], new))
            overlap += 1
    if not overlap or max_difference > 1e-6:
        raise ValueError(f"Independent extraction disagrees with the original samples: {max_difference}")
    metadata = {key: parameters[key] for key in (
        "direction", "wind_from_deg", "Uref_m_s", "Uref_height_m", "z0_m", "closure",
        "query_sha256", "protocol_sha256", "source_hashes", "rotor_D_m", "Cp", "Ct",
    )}
    metadata.update(
        checked_overlap_valid_rows=overlap, max_overlap_vector_difference_m_s=max_difference,
        format_version=1, case=case.name, iteration=iteration, role=parameters["role"],
        rotor_center_world_m=parameters["rotor_center_m"],
        disk_direction_world=study.city.transform(parameters["diskDir"], parameters["wind_from_deg"], inverse=True),
        grid_spacing_xy_m=5, heights_agl_m=HEIGHTS,
        coordinate_frame="Local east, north, up. Heights above flat ground in metres.",
        validity="1 = finite CFD sample in air, 0 = requested location without a valid sample",
        unique_positions=len(field), valid_positions=sum(v is not None for v in field.values()),
        metrics_sha256=hashlib.sha256(metrics_file.read_bytes()).hexdigest(),
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as temp:
        path = Path(temp.name)
    try:
        with path.open("wb") as stream, gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as output:
                output.write("# " + json.dumps(metadata, sort_keys=True, allow_nan=False) + "\n")
                writer = csv.writer(output, lineterminator="\n")
                writer.writerow(HEADER)
                for xyz, vector in sorted(field.items()):
                    writer.writerow((*xyz, int(vector is not None), *(vector or ("", "", ""))))
        os.replace(path, destination)
    finally:
        path.unlink(missing_ok=True)
    return metadata


class WindField:
    def __init__(self, path):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
            comment = stream.readline()
            if not comment.startswith("# "):
                raise ValueError("Missing wind field metadata")
            self.metadata = json.loads(comment[2:])
            rows = csv.reader(stream)
            if self.metadata["format_version"] != 1 or next(rows) != HEADER:
                raise ValueError("Unsupported wind field format")
            self.points = {}
            for x, y, z, valid, ue, un, uu in rows:
                xyz = (float(x), float(y), float(z))
                if valid == "1":
                    vector = (float(ue), float(un), float(uu))
                elif valid == "0" and not any((ue, un, uu)):
                    vector = None
                else:
                    raise ValueError(f"Invalid mask or velocity at {xyz}")
                if xyz in self.points or not all(map(math.isfinite, xyz + (vector or ()))):
                    raise ValueError(f"Duplicate or non-finite sample at {xyz}")
                self.points[xyz] = vector
        if (self.metadata["unique_positions"] != len(self.points)
                or self.metadata["valid_positions"] != sum(v is not None for v in self.points.values())):
            raise ValueError("Wind field counts disagree with metadata")

    @classmethod
    def for_direction(cls, direction):
        """Load a released CFD condition by compass label, e.g. N or NE."""
        folder = ROOT / "lookup"
        fields = read_json(folder / "index.json")["fields"]
        label = next((name for name in fields if name.casefold() == direction.casefold()), None)
        if label is None:
            raise ValueError(f"No released CFD field for {direction}")
        record = fields[label]
        path = folder / record["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Published {direction} field checksum mismatch")
        field = cls(path)
        if field.metadata["direction"] != label or field.metadata["wind_from_deg"] != record["wind_from_deg"]:
            raise ValueError("Lookup index and field direction disagree")
        return field

    def check_grid_and_physics(self, other):
        for key in ("Uref_m_s", "Uref_height_m", "z0_m", "closure", "rotor_center_world_m",
                    "protocol_sha256", "source_hashes", "coordinate_frame", "rotor_D_m", "Cp", "Ct"):
            if self.metadata[key] != other.metadata[key]:
                raise ValueError(f"Wind fields have different {key}")
        if self.points.keys() != other.points.keys():
            raise ValueError("Wind fields have different world grids")

    def at(self, x, y, z):
        """Return ENU velocity, None for a masked sample, or raise KeyError off the grid."""
        return self.points[x, y, z]


class WindSector:
    """Interpolate wind-aligned components at each common-air world location."""
    def __init__(self, lower, upper):
        self.lower, self.upper = lower, upper
        self.start = lower.metadata["wind_from_deg"]
        self.end = upper.metadata["wind_from_deg"]
        lower.check_grid_and_physics(upper)
        if not 0 <= self.start < self.end < 360:
            raise ValueError("Expected two ordered wind headings")

    def at(self, wind_from_deg, x, y, z):
        """Query at a single position. Use at_angle for repeated requests at one angle."""
        return self.at_angle(wind_from_deg)(x, y, z)

    def at_angle(self, wind_from_deg, *, wind_aligned=True):
        """Prepare a point query with the heading weights computed once."""
        if not self.start <= wind_from_deg <= self.end:
            raise ValueError(f"Supported sector is {self.start:g} to {self.end:g} degrees")
        if wind_from_deg == self.start:
            return self.lower.at
        if wind_from_deg == self.end:
            return self.upper.at
        w = (wind_from_deg - self.start) / (self.end - self.start)
        # R(-query) * R(endpoint) maps that endpoint vector into the queried
        # world heading while preserving its crosswind and vertical components.
        da = math.radians(self.start - wind_from_deg) if wind_aligned else 0.
        db = math.radians(self.end - wind_from_deg) if wind_aligned else 0.
        ca, sa = (1 - w) * math.cos(da), (1 - w) * math.sin(da)
        cb, sb = w * math.cos(db), w * math.sin(db)

        def query(x, y, z):
            a, b = self.lower.at(x, y, z), self.upper.at(x, y, z)
            if a is None or b is None:
                return None
            return (ca*a[0] - sa*a[1] + cb*b[0] - sb*b[1],
                    sa*a[0] + ca*a[1] + sb*b[0] + cb*b[1],
                    (1 - w)*a[2] + w*b[2])
        return query


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("export")
    build.add_argument("--case", type=Path, required=True)
    build.add_argument("--results", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    query = commands.add_parser("query")
    query.add_argument("field", help="Released compass label or a compressed CSV path")
    query.add_argument("xyz", type=float, nargs=3, metavar="METRES")
    estimate = commands.add_parser("estimate", help="Experimental linear directional estimate")
    estimate.add_argument("lower", type=Path)
    estimate.add_argument("upper", type=Path)
    estimate.add_argument("degrees", type=float)
    estimate.add_argument("xyz", type=float, nargs=3, metavar="METRES")
    args = parser.parse_args()
    if args.command == "export":
        start = time.perf_counter()
        metadata = export(args.case, args.results, args.output)
        print(json.dumps({"field": metadata, "export_elapsed_s": time.perf_counter() - start}, indent=2))
    else:
        if not all(map(math.isfinite, args.xyz)):
            parser.error("Coordinates must be finite")
        interpolated = False
        try:
            if args.command == "estimate":
                field = WindSector(WindField(args.lower), WindField(args.upper))
                degrees = args.degrees
                query_at = field.at_angle(degrees)
                interpolated = field.start < degrees < field.end
            else:
                field = (WindField.for_direction(args.field) if args.field.isalpha() and not Path(args.field).is_file()
                         else WindField(args.field))
                degrees = field.metadata["wind_from_deg"]
                query_at = field.at
        except (ValueError, OSError) as error:
            parser.error(str(error))
        try:
            velocity = query_at(*args.xyz)
        except KeyError:
            velocity, status = None, "unsampled"
        else:
            if interpolated:
                status = "estimated" if velocity is not None else "masked_endpoint"
            else:
                status = "valid_air" if velocity is not None else "masked"
        print(json.dumps({"wind_from_deg": degrees, "xyz_m": args.xyz,
                          "status": status, "U_enu_m_s": velocity}))


if __name__ == "__main__":
    main()
