# result

## 当前基准结果（2026-10-21整理）

| 目录 | 内容 | 证据编号 |
| --- | --- | --- |
| `ieee14_layer2_accuracy/`、`ieee30_layer2_accuracy/`、`ieee57_layer2_accuracy/`、`ieee118_layer2_accuracy/` | 计入PV机组无功限值时的完整计算结果：`accuracy_by_quantity.csv`（按电气量汇总）、`distribution_errors_all_objects.csv`（逐对象误差）、`fit_status_*.csv`、`fitted_models.mat`、`summary.json`、`run_manifest.json`、`pdf/`与`pdf_selected/`概率密度诊断图 | E-PF-PV-20261021、E-ACC-L2-20261021 |
| `ieee14_layer2/`、`ieee30_layer2/`、`ieee57_layer2/`、`ieee118_layer2/` | 上述计算所用的20000组检验场景及其交流Monte Carlo参照解（`states.mat`），以及当时的批量潮流核对与运行方式记录（`summary.json`、`mode_log_test.csv`）。参照解为MATPOWER 8.0牛顿-拉夫逊法，`pf.enforce_q_lims=1`、平衡机组不设无功限值，可由精度计算入口在`ReferenceRunDir`为空时重新生成 | E-PF-PV-20261021 |
| `ieee14_reactive_limit_check/` | IEEE14上8个确定性无功限值算例 | E-QLIM-IEEE14-20261021 |
| `case33bw/`、`case123/` | 模板示例目录，未用于本项目 | — |

不计无功限值的结果、云端复算包及中间计算结果已于2026-10-21删除，主要数值保留在`01_IDEA/research_trace.yaml`和`TMM_Transmission_changelog.md`。

## 模板说明

Generated numerical evidence stays under one case folder per runnable case.

```text
result/
├── case33bw/
│   ├── summary.txt
│   ├── run_manifest.json
│   └── figures/
└── case123/
    ├── summary.txt
    ├── run_manifest.json
    └── figures/
```

Each formal case run writes `run_manifest.json` during finalization. The
manifest records the entry command, Git state, MATLAB/platform information,
solver visibility, case configuration, input files, cache snapshot, output
files, formal figure manifest paths, and evidence-registration reminder.

Formal paper figures must be written to `result/<case>/figures/` through
`TMM_Transmission_utils.plotting.save_figure`. That helper exports the MATLAB source
figure, preview image, vector image, plot data CSV when provided, a JSONL
manifest, and a check report.

Before a figure is treated as formal evidence, its metadata must include the
scientific/engineering question, data files, field/unit/dimension description,
global visual encoding, target layout size, command, key parameters, and random
seed or `not_applicable`.

Do not alter, smooth, crop, or relabel data for visual appeal. If a diagnostic
figure excludes points or changes log-axis zero handling, the original full-data
figure should still be kept and the diagnostic choice must be declared in the
manifest and check report.

Do not use screenshots, copy/paste, Python/Matplotlib, Excel, Origin, or
MATLAB `subplot`/`tiledlayout` composites as the formal figure export path.

Do not add extra generated-output subfolders unless a specific case needs them
and the author chooses that layout.

Do not commit generated result files unless the project owner explicitly decides
that a specific artifact is curated evidence.
