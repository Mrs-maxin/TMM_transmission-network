# TMM_Transmission：含PV节点的交流概率潮流扩展及分布精度验证

**当前版本：2026-10-20，IEEE14基础测试。计算理论已成文；数值实现和验证未完成。**

先读[计算理论文档](01_IDEA/derivation.md)：从PV节点电压控制、机组Q限值和参考切换，说明快速潮流、普通TMM拟合、电压/线路功率CDF与PDF，以及独立完整交流MC误差验证。

## 当前研究范围

只采用标准IEEE14正常拓扑和相关负荷P/Q波动。保留原非参考Pg、Vset、Q限值、变比、充电与并联元件。当前不安排停运、区域压力、新增新能源或IEEE39/118概率算例；原其他网络文件仅保留来源参考。

初始参考1、PV2/3/6/8，9个PQ节点；13个非参考节点实虚电压共26未知量，对应13行P、9行Q和4行V²。PV的Qg为输出，触限时解除电压设定并重算。

## 当前材料与结果状态

| 内容 | 状态 |
| --- | --- |
| 计算理论、符号、近似及MC验证协议 | 已成文 |
| 官方IEEE14原文、字面JSON、来源SHA及单算例配置 | 已保存，AC尚未运行 |
| 材料/环境检查Python入口 | 已编写，未执行；本次只同步study_id |
| 混合潮流、Q控制、普通TMM、半解析查询 | 尚未实现与验证 |
| MC状态、实际PDF/CDF、误差、速度与费用 | 尚无本轮结果 |

原case14的1号机组Qg记录与Q限值不一致，需要真实复算核查参考切换。6/8号Vset高于通用Vmax需记录。全部支路rateA为0，本轮只报告MW/MVAr/MVA分布，不补造额定值或过载结论。

## 文档入口

- [研究目的](01_IDEA/story.md)、[理论正文](01_IDEA/derivation.md)、[符号与单位](01_IDEA/symbols.md)、[来源及控制依据](03_REFERENCE/ieee14_theory_sources.md)。
- [研究笔记](TMM_Transmission_note.md)、[待验证命题](01_IDEA/claims.md)、[证据状态](01_IDEA/evidence_map.md)。
- [IEEE14算例](01_IDEA/case_plan.md)、[验证标准](01_IDEA/validation_plan.md)、[分布图计划](01_IDEA/figure_plan.md)、[实现接口](01_IDEA/implementation_contract.md)。
- [算例注册](data/protocol/case_registry.json)、[网络来源SHA](data/networks/source_registry.json)、[当前状态](result/project_status.json)、[活动论文配置](02_PAPER/config/paper-profile.yaml)。

## 计算和分布比较

先做标准确定性潮流及26维方程检查，然后同输入核查电压、5机组Pg/Qg和20支路双端功率。关闭Q切换的固定模式只作诊断，正式概率比较采用匹配Q控制。成功全节点电压统一相角，保留旧参考幅值，拟合共同尺度普通多元t混合，初始固定3分量。

正式图预先选择V4/V9/V14及1–2的1端、4–9的9端、9–14的14端P/Q/MVA，比较独立MC直方图/经验CDF与TMM半解析曲线。完整表覆盖14节点和20支路两端。机组输出用于同输入校核和AC-MC统计，当前不声称仅电压TMM得到全部机组PDF。

同时报告均值、标准差、分位数、CDF误差、MC区间、PV点质量和模型泄漏。前向、拟合、积分与MC统计误差分别说明。不能只看曲线或均值宣布完成。

## 执行与保存

当前执行器离线，GitHub读写不是云端潮流计算。只在实际OpenAI托管执行器运行研究，用户Mac不运行研究计算。旧恢复已经停止，不依赖原工作目录。新执行交接见[START_NEW_CLOUD_TASK.md](03_REFERENCE/START_NEW_CLOUD_TASK.md)。

TMM_Transmission.py仅做材料检查/环境记录，尚未执行；--check-materials和--preflight不是数值研究入口。候选依赖见requirements.in，实际安装后才生成锁定版本。模板MATLAB仅供参考，未完成本地化。

每阶段在云端保存代码、输入、状态、失败、模型、种子、日志、manifest与SHA。用户明确说“Push／推送”时，把当时测试内容和结果整体同步GitHub并核对实际数组/LFS对象。此次用户明确要求补充仓库理论文档，因此保存本次文档和IEEE14对应输入配置；后续测试不自动推送。

result/materials中的旧重建/读回记录只适用于各自指定提交，不是方法验证。目标期刊暂沿用Applied Energy，当前没有投稿放行结论。
