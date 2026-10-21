# 电压越限与支路过载并发风险：文献、工程信息与可行性

日期：2026-10-19。研究对象：给定拓扑、控制规则和源荷统计条件，同一随机稳态下某节点电压越限与指定支路端过载的交集事件。当前为文献与模型分析，无研究代码和数值结果。

定位更新：本记录形成时并发为候选；随后用户已确认纳入主线。下文保留检索时的判断和读取边界，当前研究目的及算例以 story.md 和研究笔记为准；原推导附录 B 已提升为 derivation 第 9 节，沿用 B1–B19。

## 1. 检索与证据边界

本轮使用 PowerLit 文献检索技能，检索 TPWRS/TSG 索引的题名与正文开头，并补充中英文网络检索。主题包括 joint probabilistic load flow、joint distribution、joint chance constraints、voltage violation、line overload、current、simultaneous/concurrent，以及联合概率、并发越限、交集。

本地索引构建于 2026-07-03，原始 JSON 目录当前不可访问。除技能入口检索外，对索引进行了以下 AND 条件核查：

> ("joint probability" OR "joint distribution") AND voltage AND (overload OR "line flow" OR current)

索引正文通常只保存论文开头，不能用片段未出现某算法证明全文不涉及。下表注明实际读取范围；不以检索未命中证明无人研究。

## 2. 直接近邻与可核查差异

| 文献与来源 | 本轮实际读取 | 已有能力与当前对象的差异 |
| --- | --- | --- |
| Wang, Shen, Liu, Gao，Analytical Expressions for Joint Distributions in Probabilistic Load Flow，TPWRS 32(3):2473–2474，2017，[DOI](https://doi.org/10.1109/TPWRS.2016.2612881) | PowerLit 片段覆盖引言、式(1)–(7)及 GMM 联合积分开头；另核对清华作者页面的书目信息 | 已明确以多线路同时过载为动机，用 GMM、线性潮流和多重积分获得支路有功联合 PDF/CDF。联合概率、同时越限、非高斯相关输入并非新概念；本文对象可细化为非线性幅值约束下的电压与电流/MVA 交集 |
| T. Wang 等，An Improved Analytical Methodology for Joint Distribution in Probabilistic Load Flow，AECE 20(1):49–56，2020，[出版社摘要](https://www.aece.ro/abstractplus.php?article=7&number=1&year=2020)，10.4316/AECE.2020.01007 | 出版社摘要与元数据 | 改进 GMM 参数辨识并用线性潮流计算线路功率联合分布，测试 39/118 节点。本文不能把“联合分布用于输电网”视为充分增量 |
| Baker & Bernstein，Joint Chance Constraints in AC Optimal Power Flow: Improving Bounds Through Learning，TSG 10(6)，2019，[DOI](https://doi.org/10.1109/TSG.2019.2903767)，[已读作者预印本](https://arxiv.org/html/1809.04153v3) | 预印本全文，重点读交集概率与改进 Boole 界部分 | 已估计成对及高阶违约事件交集以改进联合机会约束界；示例聚焦配电网电压。不能说联合机会约束研究都忽略交集。本文目标是指定电压与过载交集的概率和后果，而非构造全局安全概率界或优化控制 |
| Jia 等，Frequency-Control-Aware Probabilistic Load Flow: An Analytical Method，TPWRS 38(6):5170–5187，2023，[DOI](https://doi.org/10.1109/TPWRS.2022.3223884) | PowerLit 正式摘要；[作者预印本全文](https://arxiv.org/html/2203.00045v1)题名为 Control-aware Probabilistic Load Flow for Transmission Systems | 已有输电网联合状态概率潮流，预印本基于控制感知的分段线性模型与校正。不能声称输电网不存在联合状态分布；本文拟补具体非线性风险区域的计算和解释。预印本细节不自动归于最终版本 |
| Hong 等，Thermal Overloading Risk Mitigation With a Semi-Analytical Probabilistic Model on Branch Current，TPWRS 36(4):3384–3393，2021，[DOI](https://doi.org/10.1109/TPWRS.2020.3048302) | PowerLit 摘要/引言；作者实验室页面核对年份与卷页，未读完整算法 | 已在线性直角坐标潮流与 GMM 下半解析计算电流越限概率。电流幅值过载的半解析计算不是空白；待补的是与节点电压事件结合的具体联合评估 |
| Gao 等，Analytical Probabilistic Power Flow and Global Sensitivity Analysis of Distribution Systems Based on Gaussian Mixture Model of Input-Output Variables，TPWRS 39(3)，2024，[DOI](https://doi.org/10.1109/TPWRS.2023.3329197) | PowerLit 摘要、引言开头及作者公开介绍，未读完整算法 | 已有非高斯相关输入、分段线性潮流及输入输出联合 PDF。联合表示具有可扩展能力，不能因其使用配电网算例就认定其原则上无法计算本研究事件 |
| Tan 等，Computationally efficient tail distribution-aware large-scale power system overloading risk assessment，Nature Communications 17:1521，2026，[原文](https://www.nature.com/articles/s41467-025-68241-y) | 前轮已核对出版社全文，重点是支路有功及尾部重采样 | 已有交流求解器代理与尾部校正。本文不能以“其他方法不保留交流非线性或不能算尾部”笼统区分；应比较具体风险对象与计算路径 |

统计依据：[Ding，On the Conditional Distribution of the Multivariate t Distribution，2016](https://doi.org/10.1080/00031305.2016.1164756)，本轮读[作者预印本](https://arxiv.org/pdf/1604.00561)的条件分布公式与证明。线性投影、条件 t 分布和混合条件权重都是既有性质，不作为新理论。

## 3. 并发风险增加什么工程信息

单项概率只说明电压和支路各自的越限频率。并发概率补充“是否在同一随机状态越限”，由此区分：无所定义越限、仅电压越限、仅过载、两类同时越限。相同边际概率可以对应不同的四类状态比例。

条件概率“过载已经发生时电压也越限的概率”可用于同时监视对象与联合告警信息解释。它描述同一状态的关联，不表示随后发生电压越限、因果传递或级联故障。

并发时的电压偏差和过载幅度分别统计，可以识别“偶尔一起越限但幅度小”和“并发时两类偏差均较大”。不把不同单位后果任意加权为综合损失。

依赖增强可能同时提高交集概率、降低至少一项越限概率，因此不说“交集大就代表全部意义上的系统风险更高”。只有线性加总单项损失的期望时，边际信息本已足够；并发概率的新增价值来自状态组成、条件关联和联合条件后果，而非必然改变这个线性期望。

## 4. 当前能支持的研究缺口

1. 部分联合分布方法针对线性化的支路有功输出，不能直接把其线性区域积分结果等同于电压幅值与电流/MVA 的非线性交集。需保留实际事件几何关系。
2. 联合机会约束中的交集估计可服务概率界，但概率界和约束满足保证不直接给出选定并发事件的完整概率及条件幅值后果。若已有方法已经输出精确联合分布，则原则上可以计算本文事件，差异转为效率、近似与输出解释。
3. 大规模采样、代理模型和联合分布方法均可作为可行替代。本文需核查联合尾部精度和完整成本，不能只展示边际分布吻合或把输入相关性等同于输出独立性。
4. 本轮未找到与“电压幅值及同端电流/MVA 并发事件 + 条件 TMM 复用二次型尾部计算”完全相同的专门方案。这是有限检索结论，不是首创证明。

## 5. 拟采用路线与可行性

选定节点复电压和支路端复电流，四个实坐标均为全网电压实向量的线性投影，因此可直接得到至多四维 TMM，无需重新拟合此局部概率对象。已有方法必须提供一致的联合参数；仅有各节点/支路边际分布不能完成本步骤。

对每个 t 分量，给定电压的两个实坐标后，电流仍为二维 t，条件自由度为原自由度加 2。其电流平方可调用既有单二次型 CDF，外层在电压越限区域作二维积分。此为半解析条件积分方案，不是已验证的闭式或加速结论。

同端 MVA 事件在标幺制下为 |V||I| 超过限值，给定 V 后只需把电流平方阈值换为 MVA 限值平方除以电压幅值平方。该路线不把四次式说成二次型，而是通过条件化复用既有计算。若监视节点不是该 MVA 支路端节点，通常需加入该端电压，维度至多增加到六维；首轮限定同端对象。

过压事件的外层区域无界，需控制尾部截断；固定 PV 幅值、参考电压或近奇异局部分布需按实际自由度处理。跨拓扑仍要更新状态联合分布，局部积分维数不随网架节点数增加，不代表全部前处理成本不随网络规模增加。

## 6. 创新性与完成判断

| 内容 | 判断 |
| --- | --- |
| 首次考虑联合概率、同时越限或非高斯相关性 | 不支持，已有直接近邻 |
| t 的投影/条件分布及四维表示是新数学 | 不支持，均为标准性质 |
| 以已发表 TMM 计算电压与电流/MVA 混合事件及条件后果 | 有工程与计算模型扩展潜力；需和上述近邻按实际事件比较 |
| 用正常和选定故障断面证明边际风险遗漏的并发信息 | 有工程贡献潜力；需真实可解释差异，而非仅报告新指标 |
| 数学可行性 | 路线清楚，假设与参数可检查 |
| 数值可行性和效率 | 尚未验证；窄尾部、相关近退化、嵌套积分与状态表示误差是主要待核查内容 |

建议先限定两个输电系统、少量预选节点与支路端组合，正常及少量非孤岛 N−1 断面。先验证低维积分，再用独立完整交流参照验证最终并发概率与后果。首轮不增加全网任意组合、动态级联、主导约束、灵敏度或优化任务。

结论：值得作为已发表 TMM 的工程评估扩展继续推进，当前不能据此认定创新充分或数值验证已完成。
