"""TLE 解析、半长轴异常剔除、SGP4 差分误差。"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sgp4.api import Satrec

from config import FORECAST_DAYS, FORECAST_TOLERANCE_DAYS, SIGMA_K

MU = 398600.4418  # WGS72 地球引力常数 km^3/s^2


@dataclass
class TLERecord:
    norad_id: int
    line1: str
    line2: str
    epoch: datetime
    satrec: Satrec
    semi_major_km: float
    is_outlier: bool = False

    @property
    def epoch_jd(self) -> tuple[float, float]:
        return self.satrec.jdsatepoch, self.satrec.jdsatepochF


def _epoch_from_satrec(sat: Satrec) -> datetime:
    jd = sat.jdsatepoch + sat.jdsatepochF
    a = int(jd + 0.5)
    frac = jd + 0.5 - a
    if a < 2299161:
        b = 0
    else:
        alpha = int((a - 1867216.25) / 36524.25)
        b = 1 + alpha - int(alpha / 4)
    c = a + b + 1524
    d = int((c - 122.1) / 365.25)
    e = int(365.25 * d)
    g = int((c - e) / 30.6001)
    day = c - e - int(30.6001 * g) + frac
    month = g - 1 if g < 14 else g - 13
    year = d - 4715 if month > 2 else d - 4716
    hours = (day - int(day)) * 24
    h, rem = divmod(hours * 3600, 3600)
    m, s = divmod(rem, 60)
    return datetime(year, month, int(day), int(h), int(m), int(s), tzinfo=timezone.utc)


def semi_major_axis_km(sat: Satrec) -> float:
    """由 TLE 平运动计算半长轴 a (km)。"""
    n_rad_s = sat.no_kozai * 2 * np.pi / 86400.0
    if n_rad_s <= 0:
        return np.nan
    return (MU / n_rad_s ** 2) ** (1 / 3)


def parse_tle_text(text: str) -> list[TLERecord]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    records: list[TLERecord] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            l1, l2 = lines[i], lines[i + 1]
            i += 2
        elif i + 2 < len(lines) and lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            l1, l2 = lines[i + 1], lines[i + 2]
            i += 3
        else:
            i += 1
            continue
        try:
            sat = Satrec.twoline2rv(l1, l2)
            a = semi_major_axis_km(sat)
            if np.isnan(a):
                continue
            records.append(TLERecord(
                int(l1[2:7]), l1, l2, _epoch_from_satrec(sat), sat, a,
            ))
        except (ValueError, IndexError):
            continue
    records.sort(key=lambda r: r.epoch)
    return records


def parse_tle_file(path: Path) -> list[TLERecord]:
    return parse_tle_text(path.read_text(encoding="utf-8", errors="ignore"))


def filter_outlier_tle(records: list[TLERecord]) -> tuple[list[TLERecord], dict]:
    """半长轴二次拟合 + 3σ 野值剔除。"""
    n = len(records)
    stats = {"total": n, "outliers": 0, "used": n}
    if n < 7:
        return records, stats

    a_vals = np.array([r.semi_major_km for r in records])
    residuals = np.full(n, np.nan)

    for i in range(3, n - 3):
        # 窗口内 7 点，跳过中心点 i 做拟合
        idx_fit = [i - 3, i - 2, i - 1, i + 1, i + 2, i + 3]
        t_fit = np.array([k - (i - 3) for k in idx_fit], dtype=float)
        a_fit = a_vals[idx_fit]
        coeffs = np.polyfit(t_fit, a_fit, 2)
        t_target = 3.0  # 第 4 个点在窗口中的位置
        a_pred = np.polyval(coeffs, t_target)
        residuals[i] = a_vals[i] - a_pred

    valid = residuals[~np.isnan(residuals)]
    if len(valid) < 3:
        return records, stats

    sigma = np.std(valid, ddof=1)
    if sigma < 1e-9:
        return records, stats

    outlier_idx = set()
    for i in range(3, n - 3):
        if not np.isnan(residuals[i]) and abs(residuals[i]) > SIGMA_K * sigma:
            outlier_idx.add(i)
            records[i].is_outlier = True

    filtered = [r for j, r in enumerate(records) if j not in outlier_idx]
    stats["outliers"] = len(outlier_idx)
    stats["used"] = len(filtered)
    stats["usage_pct"] = 100.0 * len(filtered) / n if n else 0.0
    return filtered, stats


def calculate_prediction_error(l1_old, l2_old, l1_new, l2_new, jd, fr) -> float:
    sat_old = Satrec.twoline2rv(l1_old, l2_old)
    sat_new = Satrec.twoline2rv(l1_new, l2_new)
    e1, r_pred, _ = sat_old.sgp4(jd, fr)
    e2, r_ref, _ = sat_new.sgp4(jd, fr)
    if e1 != 0 or e2 != 0:
        return np.nan
    return float(np.linalg.norm(np.array(r_pred) - np.array(r_ref)))


def _find_pairs(records: list[TLERecord], days: list[int], tol_days: float):
    pairs = []
    used = set()
    for d in days:
        target = d * 86400.0
        tol = tol_days * 86400.0
        for i, r0 in enumerate(records):
            for j in range(i + 1, len(records)):
                r1 = records[j]
                delta = (r1.epoch - r0.epoch).total_seconds()
                if abs(delta - target) <= tol:
                    key = (r0.epoch.isoformat(), r1.epoch.isoformat(), d)
                    if key not in used:
                        pairs.append((r0, r1, d))
                        used.add(key)
    return pairs


def compute_errors(records: list[TLERecord], orbit: str, norad_id: int) -> pd.DataFrame:
    pairs = _find_pairs(records, FORECAST_DAYS, FORECAST_TOLERANCE_DAYS)
    rows = []
    for r0, r1, d in pairs:
        jd, fr = r1.epoch_jd
        err = calculate_prediction_error(r0.line1, r0.line2, r1.line1, r1.line2, jd, fr)
        if np.isnan(err):
            continue
        rows.append({
            "Orbit": orbit, "NORAD_ID": norad_id, "Day": d,
            "Error_km": err,
        })
    return pd.DataFrame(rows)


def build_dataset(tle_paths: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """构建误差数据集与剔除统计表。"""
    frames, filter_rows = [], []
    for orbit, sats in tle_paths.items():
        for nid, path in sats.items():
            raw = parse_tle_file(path)
            filtered, fst = filter_outlier_tle(raw)
            filter_rows.append({
                "Orbit": orbit, "NORAD_ID": nid,
                "原始TLE数": fst["total"],
                "异常TLE数": fst["outliers"],
                "使用率_%": round(fst.get("usage_pct", 100), 2),
            })
            if len(filtered) < 2:
                continue
            df = compute_errors(filtered, orbit, nid)
            if not df.empty:
                frames.append(df)
    err_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    filt_df = pd.DataFrame(filter_rows)
    return err_df, filt_df
