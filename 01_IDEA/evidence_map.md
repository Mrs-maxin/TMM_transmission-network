# Evidence Map

> 2026-10-21整理后只保留与现存结果对应的三项证据。删除的证据（云端IEEE14核查、受控重拟合、不计无功限值的第一层结果等）的主要数值见`research_trace.yaml`节点N21—N28与N31。

本文件记录可验证研究结论与原始结果之间的对应关系。它回答一个问题：每个结论编号由哪个入口、哪份结果、哪组指标支撑。

## 证据登记表

| 证据编号 | 支撑结论编号 | 生成入口 | 原始结果文件 | 证据类型 | 关键指标/数值 | 求解器与配置 | 状态 | 备注 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E-PF-PV-20261021 | C0：含PV节点批量潮流与牛顿-拉夫逊法求解同一交流潮流方程，计入机组无功限值时运行方式与状态一致 | tests/run_tmm_transmission_accuracy.m（层2）；参照解见result/ieee*_layer2/ | result/ieee{14,30,57,118}_layer2_accuracy/<时间>/summary.json中pairedFastVersusMatpower；result/ieee{14,30,57,118}_layer2/<时间>/summary.json、mode_log_test.csv、states.mat | 20000组检验场景逐场景比较 | 四个系统运行方式20000/20000一致，复电压最大误差3.7e-9、1.4e-8、2.3e-8、4.8e-9 p.u.；运行方式数2、5、6、37，IEEE30中98.6%场景越限、最多4轮；建立展开点的牛顿-拉夫逊计算2、9、7、53次，无逐点补算；批量潮流0.3–16 s，逐点计算401–1350 s | MATLAB R2026a、MATPOWER 8.0（case14、case_ieee30、case57、case118原始数据），pf.tol=1e-10，pf.enforce_q_lims=1，平衡机组不设无功限值；ε_y=1e-8 p.u.，ε_Q=5e-6 Mvar | supported | 仅正常拓扑；负荷为有界合成波动 |
| E-QLIM-IEEE14-20261021 | C0：无功限值处理覆盖单机越上限、越下限、多机同时越限及多轮越限 | tests/check_reactive_limit_modes.m | result/ieee14_reactive_limit_check/reactive_limit_modes.csv、run_manifest.json | 确定性算例 | 8个算例最终运行方式全部一致，节点4、5、9无功负荷增加时4台机组经3轮到上限；电压最大误差2.2e-11 p.u.，机组无功最大误差1.7e-8 Mvar | 同上 | supported | 无功负荷改变量为人为设定的检验工况，不是概率场景 |
| E-ACC-L2-20261021 | C0：四个系统计入PV机组无功限值时节点电压幅值、相角与线路有功、无功分布的准确度和精度；GMM、标准TMM、分块TMM精度相当 | tests/run_tmm_transmission_accuracy.m（层2，等人数分组初值）；tests/replot_accuracy_pdf.m | result/ieee{14,30,57,118}_layer2_accuracy/<时间>/accuracy_by_quantity.csv、distribution_errors_all_objects.csv、fit_status_*.csv、fitted_models.mat、pdf/、pdf_selected/、run_manifest.json | 独立交流Monte Carlo（20000组）分布比较；概率密度诊断图 | TMM均值误差：电压幅值≤7.5e-5 p.u.、相角≤0.03°、线路有功≤0.41 MW、无功≤0.08 Mvar；相角KS≤0.025；线路功率中一般分布的线路端（72%–86%）KS中位数0.012–0.014，平顶分布的末端负荷支路（14%–24%）0.027–0.038，固定注入支路（2%–4%）约0.5；IEEE118中紧邻PV机组的16个PQ节点KS最高0.37（物理偏差≤1e-3 p.u.）；分量数14/30为20（GMM与分块TMM在IEEE30为12），57为9（GMM为6），118为6 | derivation第3、4节；建模场景5256组，定阶场景1752组，积分点50000 | supported | 当前负荷输入有界且无厚尾，分块TMM的尾部作用未检验 |

## 图表证据包

每个正式算例先写出 `result/<case>/run_manifest.json`，记录入口命令、Git、MATLAB/平台、求解器可见性、配置、输入文件、缓存快照、summary、图表证据路径和证据登记提示。图表证据包记录单张图的科学问题、数据和导出质量；运行 manifest 记录整个 case 的复现输入、运行配置和输出文件。

正式图在导出前先在 `figure_plan.md` 完成算例故事：绑定结论编号（导出 metadata 字段为 `claimId`），写清具体工程问题、本文方法改变的物理状态或工程决策、代表性对比方法、读者能够直接看到的工程效果，再确定数据、物理量和图形类型。`evidenceRole` 只用于审计分类，不构成固定出图顺序。

正式图按 case 保存到 `result/<case>/figures/`。每张图对应一个完整证据包：

| 文件 | 作用 |
| --- | --- |
| `.fig` | MATLAB 源图，便于复查和微调 |
| `.png` | 预览图 |
| `.svg` | 论文排版用矢量图 |
| `*_plot_data.csv` | 与图中曲线、柱状或散点对应的绘图数据 |
| `figure_manifest.jsonl` | 科学问题、数据文件、字段/单位/维度、视觉编码、目标尺寸、命令、参数和随机抽样设置 |
| `figure_check_report.md` | 自动检查项和人工复核项 |

## 结果解释

| 证据编号 | 支持的具体判断 | 工程场景和证据口径 | 论文中可使用的图表或数字 | 需要同步更新的文件 |
| --- | --- | --- | --- | --- |

## 状态约定

| 状态 | 含义 |
| --- | --- |
| scaffold | 模板占位，还没有真实计算 |
| running | 正在计算或等待补充输出 |
| supported | 当前证据支撑对应结论 |
| weakened | 证据只支撑较弱版本的结论 |
| revised | 证据促使结论表述、工程场景或证据状态被修正 |
| archived | 曾经有用，但不再支撑当前论文主线 |
