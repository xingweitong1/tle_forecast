# TLE 预报精度统计分析

基于 **Space-Track** 历史 TLE 与 **SGP4** 传播器，对 LEO / MEO / GEO 三类轨道卫星进行位置预报精度统计分析。

## 快速开始

克隆仓库后，使用内置 TLE 缓存即可直接运行，**无需 Space-Track 账号**：

```bash
git clone https://github.com/你的用户名/tle-forecast-homework.git
cd tle-forecast-homework
pip install -r requirements.txt
python main.py
```

Windows PowerShell：

```powershell
git clone https://github.com/你的用户名/tle-forecast-homework.git
cd tle-forecast-homework
pip install -r requirements.txt
python main.py
```

程序将读取 `data/tle/` 中 15 颗卫星（LEO / MEO / GEO 各 5 颗）过去 90 天的 TLE，完成数据清洗、误差计算与完整统计流程，结果写入 `data/processed/` 与 `output/`。

## 统计方法链

```
TLE 半长轴 3σ 清洗 → 差分法预报误差（1 / 3 / 7 天）
  → 协方差传播 Σ(Δt)=ΦΣ₀Φᵀ+Q → 增长律拟合
  → Bootstrap 置信区间（Percentile / BCa）→ Gaussian 对比
  → ANOVA / Kruskal-Wallis 轨道类型差异检验
```

## 项目结构

```
tle_forecast_homework/
├── main.py              # 主入口
├── config.py            # 卫星 NORAD ID、预报天数等参数
├── download.py          # Space-Track 下载与本地缓存
├── tle_utils.py         # TLE 解析、清洗、SGP4 差分误差
├── statistics.py        # 协方差、Bootstrap、拟合、假设检验
├── plots.py             # 图表
├── requirements.txt
├── data/
│   ├── tle/             # 已下载 TLE
│   └── processed/       # 运行后生成的 CSV
└── output/              # 图表与 report.json
```

## 运行方式

### 默认：使用本地缓存

```bash
python main.py
```

优先读取 `data/tle/` 中的 `.tle` 文件。仓库已包含 15 颗卫星数据，适合直接复现分析结果。

### 可选：从 Space-Track 重新下载

需要 [Space-Track.org](https://www.space-track.org) 账号，通过环境变量传入凭证：

Linux / macOS：

```bash
export SPACETRACK_USER=你的用户名
export SPACETRACK_PASSWORD=你的密码
python main.py --download
```

Windows PowerShell：

```powershell
$env:SPACETRACK_USER = "你的用户名"
$env:SPACETRACK_PASSWORD = "你的密码"
python main.py --download
```

`--download` 会从 Space-Track 拉取最新 90 天 TLE 并覆盖 `data/tle/`。若未设置凭证或下载失败，程序会回退到已有本地缓存。

## 依赖

Python 3.10+，主要包见 `requirements.txt`：

- `numpy`, `pandas`, `scipy`, `matplotlib`
- `sgp4`（TLE 传播）
- `requests`（Space-Track 下载）

## 输出文件

| 路径 | 说明 |
|------|------|
| `data/processed/errors_all.csv` | 全部差分预报误差样本 |
| `data/processed/summary.csv` | 按轨道类型与预报天数的汇总统计 |
| `data/processed/tle_filter_stats.csv` | TLE 清洗留存统计 |
| `data/processed/confidence_intervals.csv` | Gaussian / Bootstrap 置信区间 |
| `data/processed/hypothesis_tests.csv` | ANOVA 与 Kruskal-Wallis 检验结果 |
| `output/error_growth.png` | 误差增长律拟合图 |
| `output/error_boxplot.png` | 三轨道误差箱线图 |
| `output/tle_filter_usage.png` | TLE 清洗留存比例 |
| `output/liu2025_benchmark.png` | 与文献 RMSE 对比 |
| `output/ci_comparison_*.png` | 各轨道、各预报天数的 CI 对比图（共 9 张） |
| `output/report.json` | 结构化运行报告 |

终端会同步打印误差中位数、协方差传播、增长律、LEO 7 天 Bootstrap 对比及假设检验摘要。

