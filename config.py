"""项目参数。"""

import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TLE_DIR = ROOT / "data" / "tle"
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "output"

for d in (TLE_DIR, PROCESSED_DIR, OUTPUT_DIR):
    d.mkdir(parents=True, exist_ok=True)

FORECAST_DAYS = [1, 3, 7]
FORECAST_TOLERANCE_DAYS = 0.5
LOOKBACK_DAYS = 90
BOOTSTRAP_N = 5000
CONFIDENCE_LEVEL = 0.95
SIGMA_K = 3.0

ORBIT_SATELLITES = {
    "LEO": [25544, 20580, 43013, 48274, 33591],
    "MEO": [37753, 28915, 37846, 40105, 43001],
    "GEO": [41866, 41882, 39070, 40367, 28884],
}

SPACETRACK_USER = os.environ.get("SPACETRACK_USER", "")
SPACETRACK_PASSWORD = os.environ.get("SPACETRACK_PASSWORD", "")

LIU2025_ORIGINAL_TLE_RMSE = {
    "LEO": {"mean": 5.20},
    "MEO": {"mean": 6.63},
}


def date_range_utc():
    end = datetime.now(timezone.utc)
    start = end - timedelta(days=LOOKBACK_DAYS)
    return start, end
