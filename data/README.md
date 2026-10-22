# data

`networks/ieee14/case14.m`：MATPOWER 8.0 原始case14，未修改，SHA256 `2ffc4e1b734ae6c5e92dbe68b4e36010ed695a4bbcc4d065c74c4fbc39fcf3c1`。`networks/ieee30/case_ieee30.m`、`networks/ieee57/case57.m`、`networks/ieee118/case118.m`：MATPOWER 8.0 原始数据，未修改，SHA256依次为`4dd5c382…3269`、`2218325a…d28b`、`bc2e6f22…f4b5f7`（完整值见各运行的run_manifest.json）。`networks/pglib118/pglib_opf_case118_ieee.m`：PGLib-OPF的IEEE118数据，含线路额定容量，供后续线路过载风险计算使用，SHA256 `b1af0833…cec78`。MATPOWER与PGLib许可证见`03_REFERENCE/MATPOWER_LICENSE.txt`、`PGLIB_LICENSE.txt`。负荷样本由`TMM_Transmission_utils.scenarios.build_correlated_load_samples`按随机抽样设置重建，不另存输入文件。

`rts_gmlc/`：RTS-GMLC 2020年风电、集中式光伏与区域负荷的日前预测和实时序列及`gen.csv`、`bus.csv`，未修改，来源提交与SHA256见`rts_gmlc/SOURCE.md`，供云端阶段一建立厚尾新能源预测误差输入。

Put raw input data, case files, and small template examples here.

This directory is the complete reproducibility input set. A submitted project should
be able to start its full formal computation from the code plus files in
`data/`, after `../cache/` has been removed or renamed.

Do not make files under `cache/` required startup inputs. If a computation needs
an artifact before it can begin, keep that artifact here or make the code rebuild
it from files kept here.

Keep large or private datasets out of Git unless the project author explicitly
decides they are publishable and lightweight enough for the repository.
When public Git omits large or private inputs, the final submission or handoff
package should include the complete `data/` tree needed for reproduction.

By default, generic `.mat` files outside `data/` are ignored, while files under
`data/` are allowed so curated reproduction inputs can travel with a specific
research project. Add explicit ignore rules only for data files that must stay
private or intentionally omitted from a public Git mirror.
