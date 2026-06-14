#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TLE 预报精度统计分析。"""

import argparse
import io
import json
import sys
from pathlib import Path

if sys.platform == "win32" and hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from config import FORECAST_DAYS, OUTPUT_DIR, PROCESSED_DIR
from download import get_tle_paths
from plots import plot_boxplot, plot_ci_comparison, plot_filter_stats, plot_growth, plot_liu2025_benchmark
from statistics import (
    ci_to_row,
    compare_gaussian_bootstrap,
    covariance_vs_observed,
    fit_growth,
    orbit_tests,
    summarize,
)
from tle_utils import build_dataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="从 Space-Track 重新下载")
    args = parser.parse_args()
    np.random.seed(42)

    print("=" * 58)
    print(" TLE 预报精度统计")
    print("=" * 58)

    paths, source = get_tle_paths(force_download=args.download)
    if not paths or not any(paths.values()):
        sys.exit(1)
    print(f"数据来源: {source}")

    df, filter_df = build_dataset(paths)
    if df.empty:
        print("错误：无有效误差样本")
        sys.exit(1)

    df.to_csv(PROCESSED_DIR / "errors_all.csv", index=False)
    filter_df.to_csv(PROCESSED_DIR / "tle_filter_stats.csv", index=False)
    print(f"有效样本: {len(df)}")

    summary = summarize(df)
    summary.to_csv(PROCESSED_DIR / "summary.csv", index=False)
    print("\n--- 误差中位数 (km) ---")
    print(summary.pivot(index="Orbit", columns="Day", values="median").round(2))

    print("\n--- 协方差传播 Σ(Δt)=ΦΣ₀Φᵀ+Q ---")
    for day in FORECAST_DAYS:
        err = df[df["Day"] == day]["Error_km"].values
        c = covariance_vs_observed(err, day)
        print(f"  {day}天: 观测方差={c['observed_variance']:.4f}, 理论简化方差={c['theory_variance']:.4f}")

    growth = fit_growth(df)
    print("\n--- 误差增长律 (AIC最优) ---")
    for orbit, info in growth.items():
        print(f"  {orbit}: {info['best']}")

    ci_by_group = {}
    ci_rows = []
    for orbit in ["LEO", "MEO", "GEO"]:
        for day in FORECAST_DAYS:
            sub = df[(df["Orbit"] == orbit) & (df["Day"] == day)]["Error_km"].values
            if len(sub) < 10:
                continue
            ci = compare_gaussian_bootstrap(sub)
            if not ci:
                continue
            ci_by_group[(orbit, day)] = ci
            ci_rows.append(ci_to_row(orbit, day, ci))
            if orbit == "LEO" and day == 7:
                print(f"\n--- Gaussian vs Bootstrap (LEO 7天) ---")
                print(f"  偏度={ci['skewness']:.2f}")
                g, bp, bb = ci["gaussian"], ci["bootstrap_percentile"], ci["bootstrap_bca"]
                print(f"  Gaussian:           [{g['lower']:.2f}, {g['upper']:.2f}] km")
                print(f"  Bootstrap Percentile:[{bp['lower']:.2f}, {bp['upper']:.2f}] km")
                print(f"  Bootstrap BCa:      [{bb['lower']:.2f}, {bb['upper']:.2f}] km")

    ci_df = pd.DataFrame(ci_rows)
    ci_df.to_csv(PROCESSED_DIR / "confidence_intervals.csv", index=False)

    tests = orbit_tests(df)
    tests.to_csv(PROCESSED_DIR / "hypothesis_tests.csv", index=False)
    print("\n--- 轨道类型差异检验 ---")
    for _, r in tests.iterrows():
        p = r["anova_p"] if r["recommended"] == "ANOVA" else r["kruskal_p"]
        print(f"  {int(r['Day'])}天 {r['recommended']} p={p:.2e} → {'显著' if r['significant'] else '不显著'}")

    plot_growth(df, growth, OUTPUT_DIR)
    plot_boxplot(df, OUTPUT_DIR)
    plot_filter_stats(filter_df, OUTPUT_DIR)
    plot_liu2025_benchmark(summary, OUTPUT_DIR)
    for (orbit, day), ci in ci_by_group.items():
        plot_ci_comparison(ci, orbit, day, OUTPUT_DIR)

    report = {
        "data_source": source,
        "n_samples": len(df),
        "summary": summary.to_dict(orient="records"),
        "hypothesis_tests": tests.to_dict(orient="records"),
    }
    (OUTPUT_DIR / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(f"\n完成 → output/ ")


if __name__ == "__main__":
    main()
