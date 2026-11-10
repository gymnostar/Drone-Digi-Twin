"""Summarise an EDMO wind CSV: coverage per year, calms, wind rose, speed spread.

The point: pick CFD inlet cases from what the wind actually does here,
instead of guessing a speed and direction.
Usage: python3 summarise_edmo.py ../raw/EDMO_2016-2025.csv
"""
import csv
import statistics
import sys
from collections import Counter

KT = 0.514444  # knots -> m/s
SECTORS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]

rows = list(csv.DictReader(open(sys.argv[1])))
per_year = Counter(r["valid"][:4] for r in rows)
print("reports per year (2 per hour max = 17520):", dict(sorted(per_year.items())))

speeds, rose, calm, variable = [], Counter(), 0, 0
rose_speed = {s: [] for s in SECTORS}
for r in rows:
    if r["sknt"] == "M":
        continue
    u = float(r["sknt"]) * KT
    speeds.append(u)
    if u == 0:
        calm += 1
    elif r["drct"] in ("M", "0.00"):
        variable += 1  # METAR "VRB": light wind, direction not defined
    else:
        s = SECTORS[int(((float(r["drct"]) + 22.5) % 360) // 45)]
        rose[s] += 1
        rose_speed[s].append(u)

n = len(speeds)
print(f"\n{n} reports with speed | calm {calm / n:.1%} | variable direction {variable / n:.1%}")
q = statistics.quantiles(speeds, n=20)
print(f"speed m/s: mean {statistics.mean(speeds):.1f} | median {q[9]:.1f} | 90th pct {q[17]:.1f} | 95th pct {q[18]:.1f} | max {max(speeds):.1f}")
print("\nwind FROM | share of all reports | mean m/s | 90th pct m/s")
for s in SECTORS:
    v = rose_speed[s]
    print(f"  {s:>3}     | {len(v) / n:6.1%} {'#' * round(60 * len(v) / n):<25}| {statistics.mean(v):4.1f}     | {statistics.quantiles(v, n=10)[8]:4.1f}")
