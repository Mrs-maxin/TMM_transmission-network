# 输电网风险解析评估：近邻文献与候选工程扩展

日期：2026-10-19。目的：界定已有研究覆盖及本文可完成的评估增量，不作首创或投稿充分性判定。

定位更新：下文“分级首选、并发候选”为当时筛选记录；用户后续已将并发风险纳入主线。当前目的与算例见 story.md 和研究笔记，本记录保留原文献比较事实。

## 检索范围与读取边界

- 使用 PowerLit 文献检索技能，核查方法 canon，检索本地 TPWRS/TSG 分片索引的题名与正文开头。索引构建于 2026-07-03；当前配置的原始 JSON 目录不可访问，因此索引片段不冒充全文。
- 联网核对作者机构记录和出版社原文。以下明确区分全文页面、摘要/索引片段及题名线索。
- 本文初版未专门查新并发事件；后续专题检索已归档至 [并发风险专题](concurrent_voltage_overload_review_20261019.md)，核查了联合分布及交集估计近邻。有限检索不构成领域空白或首创证明。

## 直接相关研究

| 文献 | 本轮实际读取 | 可支持的比较事实 | 本文定位 |
| --- | --- | --- | --- |
| Li 等，Transmission Line Overload Risk Assessment for Power Systems With Wind and Load-Power Generation Correlation，TSG 6(3):1233–1242，2015，[DOI](https://doi.org/10.1109/TSG.2014.2387281) | PowerLit 开头片段与[作者机构摘要/元数据](https://researchwith.stevens.edu/en/publications/transmission-line-overload-risk-assessment-for-power-systems-with/) | 点估计评估过载概率，严重度函数表征后果，并考虑源荷及风资源相关性 | 概率与后果结合、相关输入均不能笼统作为新问题 |
| Xu 等，Risk Assessment of Rare Events in Probabilistic Power Flow via Hybrid Multi-Surrogate Method，TSG 11(2)，2020，[DOI](https://doi.org/10.1109/TSG.2019.2940928) | 本轮 PowerLit 摘要与引言开头；原文服务器本轮读取失败 | PCE 多代理与尾部高保真重评估处理罕见事件、线路退出 | 尾部偏差与故障拓扑已有专门研究；本文比较概率积分途径与所选工程输出 |
| Tan 等，Computationally efficient tail distribution-aware large-scale power system overloading risk assessment，Nature Communications 17:1521，2026，[出版社原文](https://www.nature.com/articles/s41467-025-68241-y) | 读取出版社全文 HTML，重点核查式(1)–(7)、输出与尾部重采样段落 | 已有交流求解器代理、预设故障、相关非高斯输入和尾部校正，主要过载量为支路有功 | 本文可相对这篇工作扩展电压及电流/实际 MVA 风险；不泛化为所有既有研究遗漏这些量 |
| Elazab 等，Novel risk index integrating practical operation limits enhances probabilistic contingency ranking for large-scale photovoltaic plant planning，Scientific Reports 14:10009，2024，[出版社原文](https://www.nature.com/articles/s41598-024-60024-7) | 读取出版社摘要及实际运行限值段落 | 已涉及电压、线路等级及允许纠正时间 | 分级本身并非空白；本研究候选增量是具体区间概率与幅值后果的解析计算 |
| Ding 等，Transmission Line Overload Risk Assessment Considering Dynamic Line Rating Mechanism in a High-Wind-Penetrated Power System: A Data-Driven Approach，2022，[IEEE 页面](https://ieeexplore.ieee.org/document/9714814/)，10.1109/TSTE.2022.3151921 | 题名、DOI 及可获取摘要开头；未读完整算法 | DLR 与风电条件下的过载风险是已有方向 | 天气条件下限值变化可作应用候选，但不能称首次考虑降额或 DLR |

## 可提出的工程计算问题

给定正常或故障后断面和源荷联合分布，按真实电压、支路限值计算越限概率与幅值后果，并分清轻微与严重越限。与上述工作相比，本文是采用已发表 TMM 的风险评估应用与工程扩展，不是首次将概率潮流用于输电网风险。

二次型与尾部准确性属于继承能力；本文所需证据是工程阈值下风险结果是否可靠、风险输出是否补充有用信息，以及完整计算成本是否可接受。状态分布构建成本与所有必要交流校验均应计入。

## 候选研究点及取舍

| 候选 | 可补充的工程信息 | 新增工作 | 当前建议 |
| --- | --- | --- | --- |
| 实际限值与分级风险 | 有功代理是否遗漏电流/MVA 过载，轻微与严重越限是否被总概率混淆 | 固定工程等级的区间概率、分项后果、匹配额定值的参照比较；MVA 时须补联合计算 | 低成本首选，不把分级本身作为新理论 |
| 选定节点与支路并发风险 | 低电压与过载是否在同一随机状态发生 | 少量对象对的联合概率与条件概率、低维积分及独立参照 | 增强候选；专题近邻和条件积分模型已整理，数值与工程效果待核查 |
| 天气条件源荷与线路限值共同影响 | 评估同一天气下输入压力和限值变化 | 天气条件化数据与有依据的限值模型 | 数据满足后再考虑，不列首轮必需 |
| 正常与选定 N−1 对照 | 拓扑变化后的风险位置与后果 | 对少量拓扑更新状态分布，按条件概率报告 | 验证设置，不单独包装创新 |

候选点均未实现或验证；不同时承诺全部开展。分级限值应有工程含义及适用条件，不是任意阈值灵敏度扫描。并发风险是统计同时事件，不等同于因果、级联或动态失稳。

2024 年分级近邻进一步表明：分级适合补足工程评估输出，单独划分等级不足以承载新意。需要更明确新增风险对象时，优先核查少量节点与支路并发事件；当前不宣称其首次被研究。
