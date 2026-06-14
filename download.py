"""Space-Track 历史 TLE 下载。"""

import time
from datetime import datetime
from pathlib import Path

import requests

from config import (
    LOOKBACK_DAYS,
    ORBIT_SATELLITES,
    SPACETRACK_PASSWORD,
    SPACETRACK_USER,
    TLE_DIR,
    date_range_utc,
)

LOGIN_URL = "https://www.space-track.org/ajaxauth/login"
QUERY_BASE = "https://www.space-track.org/basicspacedata/query"


class SpaceTrackClient:
    def __init__(self, identity: str, password: str):
        self.session = requests.Session()
        self.identity = identity
        self.password = password
        self._logged_in = False

    def login(self) -> bool:
        resp = self.session.post(
            LOGIN_URL,
            data={"identity": self.identity, "password": self.password},
            timeout=30,
        )
        self._logged_in = resp.status_code == 200
        return self._logged_in

    def query_gp_history(self, norad_id: int, start: datetime, end: datetime) -> str:
        if not self._logged_in and not self.login():
            raise ConnectionError("Space-Track 登录失败")
        url = (
            f"{QUERY_BASE}/class/gp_history/"
            f"NORAD_CAT_ID/{norad_id}/"
            f"EPOCH/{start:%Y-%m-%d}--{end:%Y-%m-%d}/"
            f"orderby/EPOCH%20asc/format/tle"
        )
        resp = self.session.get(url, timeout=120)
        resp.raise_for_status()
        return resp.text.strip()

    def query_gp_recent(self, norad_id: int, days: int = 90) -> str:
        if not self._logged_in and not self.login():
            raise ConnectionError("Space-Track 登录失败")
        url = (
            f"{QUERY_BASE}/class/gp/"
            f"NORAD_CAT_ID/{norad_id}/"
            f"EPOCH/>now-{days}/"
            f"orderby/EPOCH%20asc/format/tle"
        )
        resp = self.session.get(url, timeout=120)
        resp.raise_for_status()
        return resp.text.strip()


def _expected_count() -> int:
    return sum(len(v) for v in ORBIT_SATELLITES.values())


def _count(paths: dict) -> int:
    return sum(len(v) for v in paths.values())


def load_cached() -> dict[str, dict[int, Path]]:
    """从 data/tle/ 加载本地缓存。"""
    paths: dict[str, dict[int, Path]] = {k: {} for k in ORBIT_SATELLITES}
    for orbit, ids in ORBIT_SATELLITES.items():
        for nid in ids:
            matches = sorted(TLE_DIR.glob(f"{nid}_*.tle"))
            if matches and matches[0].stat().st_size > 0:
                paths[orbit][nid] = matches[0]
    return paths


def download_one(
    norad_id: int,
    start: datetime,
    end: datetime,
    client: SpaceTrackClient,
    save_dir: Path | None = None,
) -> Path | None:
    save_dir = save_dir or TLE_DIR
    out = save_dir / f"{norad_id}_{start:%Y%m%d}_{end:%Y%m%d}.tle"
    if out.exists() and out.stat().st_size > 0:
        return out
    try:
        text = client.query_gp_history(norad_id, start, end)
        if not text:
            text = client.query_gp_recent(norad_id, days=LOOKBACK_DAYS)
        if not text:
            return None
        out.write_text(text + "\n", encoding="utf-8")
        return out
    except requests.RequestException:
        return None


def download_all_satellites(start: datetime, end: datetime, delay: float = 1.0) -> dict[str, dict[int, Path]]:
    """循环下载 15 颗卫星 3 个月历史 TLE。"""
    if not SPACETRACK_USER or not SPACETRACK_PASSWORD:
        return {}
    client = SpaceTrackClient(SPACETRACK_USER, SPACETRACK_PASSWORD)
    if not client.login():
        return {}

    print(f"[Space-Track] 下载 {start.date()} ~ {end.date()}，共 {_expected_count()} 颗卫星")
    paths: dict[str, dict[int, Path]] = {k: {} for k in ORBIT_SATELLITES}
    for orbit, ids in ORBIT_SATELLITES.items():
        for nid in ids:
            p = download_one(nid, start, end, client)
            if p:
                paths[orbit][nid] = p
            time.sleep(delay)
    return paths


def get_tle_paths(force_download: bool = False) -> tuple[dict[str, dict[int, Path]], str]:
    """获取 TLE 路径：优先本地缓存，必要时从 Space-Track 下载。"""
    expected = _expected_count()
    start, end = date_range_utc()

    if not force_download:
        cached = load_cached()
        if _count(cached) == expected:
            print(f"[本地缓存] Space-Track TLE：{expected} 颗卫星")
            return cached, f"spacetrack 本地缓存 ({start:%Y%m%d}-{end:%Y%m%d})"

    if SPACETRACK_USER and SPACETRACK_PASSWORD:
        paths = download_all_satellites(start, end)
        if _count(paths) > 0:
            return paths, f"spacetrack 在线下载 ({start:%Y%m%d}-{end:%Y%m%d})"
        cached = load_cached()
        if _count(cached) > 0:
            print("[警告] 下载失败，使用已有缓存")
            return cached, "spacetrack 部分缓存"

    cached = load_cached()
    if _count(cached) == expected:
        print(f"[本地缓存] Space-Track TLE：{expected} 颗卫星")
        return cached, f"spacetrack 本地缓存 ({start:%Y%m%d}-{end:%Y%m%d})"
    if _count(cached) > 0:
        print(f"[警告] 缓存不完整 ({_count(cached)}/{expected})，使用已有数据")
        return cached, f"spacetrack 部分缓存 ({_count(cached)}颗)"

    print("错误：无 TLE 数据。")
    print("  请将 Space-Track 下载的 .tle 放入 data/tle/，或设置凭证后运行 python main.py --download")
    return {}, ""
