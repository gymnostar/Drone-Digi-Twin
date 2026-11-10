"""Download EDMO (Oberpfaffenhofen airfield) wind reports from the Iowa Environmental Mesonet.

Why this source: EDMO's METAR station sits ~300 m outside the DLR Oberpfaffenhofen scene
(48.081 N, 11.283 E), so the same series can drive the CFD inlet and validate it.
METAR wind is measured at 10 m above ground, as a 10-minute average.

Usage: python3 fetch_edmo.py 2016 2025   -> DLR/wind-data/raw/EDMO_2016-2025.csv
Stdlib only, so it runs without installing anything.
"""
import sys
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"


def fetch(first_year: int, last_year: int) -> Path:
    params = [
        ("station", "EDMO"),
        ("data", "drct"),   # direction the wind blows FROM, degrees true
        ("data", "sknt"),   # mean speed, knots (1 kt = 0.5144 m/s)
        ("data", "gust"),   # gust, knots, only reported when gusty
        ("year1", first_year), ("month1", 1), ("day1", 1),
        ("year2", last_year + 1), ("month2", 1), ("day2", 1),
        ("tz", "Etc/UTC"), ("format", "onlycomma"), ("latlon", "no"),
        ("missing", "M"), ("report_type", 3), ("report_type", 4),  # routine + special METARs
    ]
    out = Path(__file__).resolve().parent.parent / "raw" / f"EDMO_{first_year}-{last_year}.csv"
    with urllib.request.urlopen(f"{URL}?{urllib.parse.urlencode(params)}", timeout=300) as r:
        out.write_bytes(r.read())
    return out


if __name__ == "__main__":
    first, last = (int(a) for a in sys.argv[1:3])
    path = fetch(first, last)
    print(f"wrote {path} ({sum(1 for _ in path.open()) - 1} reports)")
