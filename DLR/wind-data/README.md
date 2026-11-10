# Wind data — DLR Oberpfaffenhofen (EDMO)

The DLR scene's terrain STL header gives origin E670000 N5328000 (UTM32), which is DLR Oberpfaffenhofen.
The airfield's METAR station EDMO (48.081 N, 11.283 E) lies just outside the scene: its archive coordinate is about 50 m south of the scene edge, and the actual mast is about 300 m out (measured by Henry). The terrain between them is flat and open, so it works as the inlet reference. To validate against it, extend the CFD domain so it includes the mast.

- **Source:** Iowa Environmental Mesonet ASOS/METAR archive (`scripts/fetch_edmo.py`).
- **What it is:** 10 m above ground, 10-minute mean speed and direction ("FROM"), every 30 min, in knots, UTC.
- **Rebuild:** `python3 scripts/fetch_edmo.py 2016 2025`, then `python3 scripts/summarise_edmo.py raw/EDMO_2016-2025.csv`.

## Caveats
- 2016–2018 report only about 7.3k per year, roughly daytime opening hours. Full 24 h coverage starts mid-2019.
  For statistics, prefer 2020–2025, otherwise daytime (windier) hours are over-weighted.
- 11.8% of reports are "variable direction" (light wind, VRB). They have a speed but no sector.
- No DWD synoptic station operates at Oberpfaffenhofen. The nearest current DWD hourly wind station is München-Stadt (03379), about 20 km away in a city, so it's a poor substitute.

## Summary 2016–2025 (`summarise_edmo.py`)
Mean 3.2 m/s, median 2.6, 90th percentile 6.2, 95th percentile 7.2, max 19.5 m/s. Prevailing wind is from the SW (23.5%) and W (17.4%), which are also the strongest sectors (90th percentile about 8 m/s).

## Second sensor: DLR weather station (candidate for validation)
- **Position:** 48.087 N, 11.280 E, 578 m asl (from the wdc.dlr.de/data_products/local_weather plots).
  In UTM32 that's E 669787 N 5328485, i.e. scene-local x = -213 m, y = +485 m. It's **inside** the scene, about 700 m NNW of EDMO.
  Three decimals means ±50 m. The mast height and whether it's on a roof are unknown.
- **Data:** the site only publishes PNG plots, no download. The wind plot's last update was 2024-10-28 and shows about 0 m/s, so the anemometer may be out of service.
  Ask DLR (Prof. Rini) for the raw history and the mast details.
