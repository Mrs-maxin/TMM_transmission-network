# tests

本目录保留给具体项目的临时验证、调参检查和本地调试入口。模板不附带默认 MATLAB 测试函数或统一测试入口，部署和复制模板时不得自动执行 MATLAB。测试结果只有通过已说明的验证入口生成，并登记到 `01_IDEA/evidence_map.md` 后，才可作为论文证据。

当前保留三个计算入口（2026-10-21整理）。依赖MATLAB统计工具箱、MATPOWER 8.0（`4-PHD/00codezip/matpower8.0`）和本地AQPPF_TMM（`[3]博三/[1]ACPPF_TMM`，数值代码同`0d5c999d`），路径在各入口的配置段中设置。

| 入口 | 内容 | 结果目录 |
| --- | --- | --- |
| `run_tmm_transmission_accuracy.m(networkName, layer, Initialization=, ComponentCounts=, ReferenceRunDir=)` | IEEE14/30/57/118的完整计算：源荷场景、批量潮流（层2计入PV机组无功限值）、GMM/标准TMM/分块TMM、电压幅值相角与线路功率分布、与交流Monte Carlo比较、概率密度图。`ReferenceRunDir`为空时重新计算20000组交流Monte Carlo | `result/<系统>_layer<层>_accuracy/<时间>/` |
| `check_reactive_limit_modes.m` | IEEE14上8个确定性无功限值算例，比较批量潮流与牛顿-拉夫逊法的最终运行方式、电压和机组无功 | `result/ieee14_reactive_limit_check/` |
| `replot_accuracy_pdf.m(runDir, objects)` | 由已完成精度计算的保存模型重画指定节点或线路的概率密度图 | `<runDir>/pdf_selected/` |

当前基准结果由`run_tmm_transmission_accuracy(<系统>, 2, Initialization="balanced-physical", ComponentCounts=..., ReferenceRunDir=...)`生成，分量数候选IEEE14/30为3–20，IEEE57/118为3–12。单个系统用时约10–30分钟，复用交流Monte Carlo时不含参照计算时间。
