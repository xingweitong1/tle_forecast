"""绘图。"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import LIU2025_ORIGINAL_TLE_RMSE, OUTPUT_DIR

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def plot_growth(df, growth, out_dir: Path) -> Path:
    colors = {"LEO": "red", "MEO": "green", "GEO": "blue"}
    fig, ax = plt.subplots(figsize=(9, 5), dpi=120)
    days = sorted(df["Day"].unique())
    for orbit in ["LEO", "MEO", "GEO"]:
        sub = df[df["Orbit"] == orbit]
        if sub.empty:
            continue
        med = sub.groupby("Day")["Error_km"].median().reindex(days)
        ax.scatter(days, med, color=colors[orbit], s=70, label=orbit)
        if orbit in growth and "polynomial" in growth[orbit]["models"]:
            c = growth[orbit]["models"]["polynomial"]["coeffs"]
            x = np.linspace(0.5, 7.5, 50)
            ax.plot(x, np.poly1d(c)(x), color=colors[orbit], linestyle="--")
    ax.set_xlabel("预报时长 (天)")
    ax.set_ylabel("位置误差中位数 (km)")
    ax.set_title("误差增长律拟合（LEO / MEO / GEO）")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    p = out_dir / "error_growth.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def plot_boxplot(df, out_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=120, sharey=True)
    for ax, day in zip(axes, sorted(df["Day"].unique())):
        sub = df[df["Day"] == day]
        ax.boxplot(
            [sub[sub["Orbit"] == o]["Error_km"].values for o in ["LEO", "MEO", "GEO"]],
            tick_labels=["LEO", "MEO", "GEO"],
        )
        ax.set_title(f"{day} 天")
        ax.grid(True, axis="y", alpha=0.4)
    fig.suptitle("三轨道预报误差分布")
    fig.tight_layout()
    p = out_dir / "error_boxplot.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def plot_ci_comparison(ci: dict, orbit: str, day: int, out_dir: Path) -> Path:
    """Gaussian vs Bootstrap 置信区间对比图。"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    labels, est, lo, hi = [], [], [], []
    g = ci["gaussian"]
    labels.append("Gaussian")
    est.append(g["estimate"])
    lo.append(g["lower"])
    hi.append(g["upper"])
    bp = ci["bootstrap_percentile"]
    labels.append("Bootstrap\n(Percentile)")
    est.append(bp["estimate"])
    lo.append(bp["lower"])
    hi.append(bp["upper"])
    bc = ci["bootstrap_bca"]
    labels.append("Bootstrap\n(BCa)")
    est.append(bc["estimate"])
    lo.append(bc["lower"])
    hi.append(bc["upper"])

    x = np.arange(len(labels))
    err = [[e - l for e, l in zip(est, lo)], [u - e for e, u in zip(est, hi)]]
    ax.bar(x, est, color="steelblue", alpha=0.7)
    ax.errorbar(x, est, yerr=err, fmt="none", color="black", capsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("误差均值 (km)")
    ax.set_title(f"{orbit} {day}天：Gaussian vs Bootstrap 95% CI")
    ax.grid(True, axis="y", alpha=0.4)
    fig.tight_layout()
    p = out_dir / f"ci_comparison_{orbit}_day{day}.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def plot_liu2025_benchmark(summary, out_dir: Path) -> Path:
    """与文献原始 TLE RMSE 比较。"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    orbits = ["LEO", "MEO", "GEO"]
    day = 7
    x = np.arange(3)
    ours, liu = [], []
    for o in orbits:
        row = summary[(summary["Orbit"] == o) & (summary["Day"] == day)]
        ours.append(row["median"].values[0] if len(row) else np.nan)
        liu.append(LIU2025_ORIGINAL_TLE_RMSE[o]["mean"] if o in LIU2025_ORIGINAL_TLE_RMSE else np.nan)
    w = 0.35
    ax.bar(x - w / 2, ours, w, label="本研究(7天中位数)")
    ax.bar(x + w / 2, liu, w, label="Liu2025 原始TLE RMSE")
    ax.set_xticks(x)
    ax.set_xticklabels(orbits)
    ax.set_ylabel("位置误差 (km)")
    ax.set_title("与《空间科学学报》Liu2025 预报精度对标")
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", alpha=0.4)
    fig.tight_layout()
    p = out_dir / "liu2025_benchmark.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def plot_filter_stats(filter_df: pd.DataFrame, out_dir: Path) -> Path:
    """TLE 历元清洗后留存比例。"""
    fig, ax = plt.subplots(figsize=(10, 4), dpi=120)
    labels = [f"{int(r.NORAD_ID)}\n({r.Orbit})" for r in filter_df.itertuples()]
    usage = filter_df["使用率_%"].values
    colors = {"LEO": "#e74c3c", "MEO": "#27ae60", "GEO": "#3498db"}
    bar_colors = [colors.get(r.Orbit, "steelblue") for r in filter_df.itertuples()]
    ax.bar(labels, usage, color=bar_colors, alpha=0.85)
    ax.axhline(87, color="gray", linestyle="--", linewidth=1.5, label="Liu2025 文献平均留存≈87%")
    ax.set_xlabel("NORAD ID（轨道类型）")
    ax.set_ylabel("TLE 历元留存比例 (%)")
    ax.set_title("数据清洗后 TLE 留存比例（半长轴 3σ 剔除）")
    ax.set_ylim(80, 102)
    ax.legend()
    ax.grid(True, axis="y", alpha=0.4)
    fig.tight_layout()
    p = out_dir / "tle_filter_usage.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p
