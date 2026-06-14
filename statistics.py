"""统计方法链：协方差传播、Bootstrap、增长律、假设检验。"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

from config import BOOTSTRAP_N, CONFIDENCE_LEVEL, FORECAST_DAYS


# 汇总 
def summarize(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["Orbit", "Day"])["Error_km"]
        .agg(count="count", mean="mean", median="median", std="std")
        .reset_index()
    )


# 协方差传播 Σ(Δt)=ΦΣ₀Φᵀ+Q 
def simplified_covariance(t_days: float, sigma0: float = 0.5, q: float = 0.1) -> float:
    phi = 1.0 + 0.1 * t_days
    return (phi * sigma0 * phi) + q * t_days


def covariance_vs_observed(errors: np.ndarray, t_days: float) -> dict:
    obs_var = float(np.var(errors, ddof=1)) if len(errors) > 1 else np.nan
    theory = simplified_covariance(t_days, sigma0=obs_var * 0.3 if obs_var > 0 else 0.5)
    return {"t_days": t_days, "observed_variance": obs_var, "theory_variance": theory}


# Bootstrap：百分位 + BCa 
def _boot_samples(data, n_boot, stat_func=np.mean):
    n = len(data)
    out = np.empty(n_boot)
    for b in range(n_boot):
        out[b] = stat_func(np.random.choice(data, n, replace=True))
    return out


def bootstrap_percentile_ci(data, n_boot=BOOTSTRAP_N, conf=CONFIDENCE_LEVEL):
    data = np.asarray(data, float)
    n = len(data)
    boot = _boot_samples(data, n_boot)
    alpha = (1 - conf) / 2
    m = np.mean(data)
    return m, np.percentile(boot, alpha * 100), np.percentile(boot, (1 - alpha) * 100)


def bootstrap_bca_ci(data, n_boot=BOOTSTRAP_N, conf=CONFIDENCE_LEVEL):
    data = np.asarray(data, float)
    n = len(data)
    theta = np.mean(data)
    boot = _boot_samples(data, n_boot)
    prop = np.clip(np.mean(boot < theta), 1e-6, 1 - 1e-6)
    z0 = stats.norm.ppf(prop)
    jack = np.array([np.mean(np.delete(data, i)) for i in range(n)])
    jm = np.mean(jack)
    num = np.sum((jm - jack) ** 3)
    den = 6.0 * (np.sum((jm - jack) ** 2) ** 1.5 + 1e-12)
    a = num / den if den else 0.0
    za, z1 = stats.norm.ppf((1 - conf) / 2), stats.norm.ppf((1 + conf) / 2)

    def adj(z):
        nz = z0 + z
        return stats.norm.cdf(z0 + nz / (1 - a * nz + 1e-12))

    lo = np.percentile(boot, np.clip(adj(za) * 100, 0.1, 99.9))
    hi = np.percentile(boot, np.clip(adj(z1) * 100, 0.1, 99.9))
    return theta, lo, hi


def gaussian_ci(data, conf=CONFIDENCE_LEVEL):
    data = np.asarray(data, float)
    n = len(data)
    m = np.mean(data)
    z = stats.norm.ppf(1 - (1 - conf) / 2)
    se = np.std(data, ddof=1) / np.sqrt(n)
    return m, m - z * se, m + z * se


def compare_gaussian_bootstrap(data) -> dict:
    """Gaussian vs Bootstrap（百分位 + BCa）完整对比。"""
    data = np.asarray(data, float)
    data = data[~np.isnan(data)]
    if len(data) < 5:
        return {}
    gm, glo, ghi = gaussian_ci(data)
    pm, plo, phi = bootstrap_percentile_ci(data)
    bm, blo, bhi = bootstrap_bca_ci(data)
    return {
        "n": len(data),
        "sample_mean": float(gm),
        "sample_median": float(np.median(data)),
        "skewness": float(stats.skew(data)),
        "gaussian": {"estimate": gm, "lower": glo, "upper": ghi},
        "bootstrap_percentile": {"estimate": pm, "lower": plo, "upper": phi},
        "bootstrap_bca": {"estimate": bm, "lower": blo, "upper": bhi},
    }


def ci_to_row(orbit: str, day: int, ci: dict) -> dict:
    return {
        "Orbit": orbit, "Day": day,
        "gauss_lo": ci["gaussian"]["lower"], "gauss_hi": ci["gaussian"]["upper"],
        "boot_pct_lo": ci["bootstrap_percentile"]["lower"],
        "boot_pct_hi": ci["bootstrap_percentile"]["upper"],
        "boot_bca_lo": ci["bootstrap_bca"]["lower"],
        "boot_bca_hi": ci["bootstrap_bca"]["upper"],
        "skewness": ci["skewness"],
        "mean": ci["sample_mean"],
        "median": ci["sample_median"],
    }


# 增长律拟合
def _poly2(t, a, b, c):
    return a * t ** 2 + b * t + c


def _power(t, a, b):
    return a * np.power(t, b)


def fit_growth(df: pd.DataFrame) -> dict:
    result = {}
    for orbit in df["Orbit"].unique():
        med = df[df["Orbit"] == orbit].groupby("Day")["Error_km"].median()
        days, errs = med.index.values.astype(float), med.values.astype(float)
        if len(days) < 2:
            continue
        coeffs = np.polyfit(days, errs, min(2, len(days) - 1))
        models = {"polynomial": {"coeffs": coeffs.tolist(), "aic": 0}}
        try:
            popt, _ = curve_fit(_power, days, errs, p0=[1.0, 1.0], maxfev=5000)
            pred = _power(days, *popt)
            aic = len(days) * np.log(np.sum((errs - pred) ** 2) / len(days) + 1e-12) + 4
            models["power_law"] = {"params": popt.tolist(), "aic": aic}
        except (RuntimeError, ValueError):
            pass
        pred_p = np.poly1d(coeffs)(days)
        models["polynomial"]["aic"] = len(days) * np.log(np.sum((errs - pred_p) ** 2) / len(days) + 1e-12) + 6
        best = min(models, key=lambda k: models[k]["aic"])
        result[orbit] = {"best": best, "models": models, "days": days.tolist(), "medians": errs.tolist()}
    return result


# ANOVA + Kruskal-Wallis
def _normal_ok(groups: list[np.ndarray]) -> bool:
    for g in groups:
        if len(g) < 8:
            return False
        s = g if len(g) <= 5000 else np.random.choice(g, 5000, replace=False)
        if stats.shapiro(s).pvalue < 0.05:
            return False
    return True


def orbit_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for day in FORECAST_DAYS:
        groups = [df[(df["Orbit"] == o) & (df["Day"] == day)]["Error_km"].dropna().values for o in ["LEO", "MEO", "GEO"]]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) < 2:
            continue
        normal = _normal_ok(groups)
        f, p_a = stats.f_oneway(*groups) if all(len(g) >= 2 for g in groups) else (np.nan, np.nan)
        h, p_k = stats.kruskal(*groups)
        rows.append({
            "Day": day, "normal_assumption": normal,
            "recommended": "ANOVA" if normal else "Kruskal-Wallis",
            "anova_F": f, "anova_p": p_a,
            "kruskal_H": h, "kruskal_p": p_k,
            "significant": p_k < 0.05,
        })
    return pd.DataFrame(rows)
