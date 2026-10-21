# 输电网尾部风险评估框架、方法比较与投稿定位

> 最新修订：用户现选定 Applied Energy，并提供 AQPPF_TMM。源库当前前向限非参考 PQ；默认多尺度 t 与普通 t 分开。后续以研究笔记第14节和 03_REFERENCE/aqppf_tmm_transmission_extension_20261019.md 为当前继承与实施口径，本文保留原检索及前序判断。
日期：2026-10-19。研究设计与文献核查记录，没有研究代码、交流计算或数值结论。用户确认：英文 SCI，JCR Q1/Q2；Sustainable Energy 指 IEEE Transactions on Sustainable Energy；排除 TPWRS。

## 1. 当前论文要回答的问题

在给定源荷联合分布、预测时窗、拓扑、控制和限值下，电压越限与支路过载的尾部概率及期望越限程度能否准确、经济地计算？已有 TMM 提供联合表示、二次型和条件尾部内核，本文完成风险事件及后果计算、低维流程和输电网端到端比较。

并发风险保留为联合尾部对象和工程输出；传播矩阵只解释共同波动。二者的任务规模与它们对精度主线的作用相称，不要求为投稿另造分布定理、控制模型或因果发现。

## 2. 常规输电网风险评估体系

流程为：评估时窗和计划运行方式 → 源荷随机输入及故障集合 → 各状态交流潮流与既定控制 → 电压/负载约束和严重程度 → 条件风险 → 有依据的概率加权汇总 → 运行监视或后续措施评估。

1. 确定性 N−1 回答预定状态和单个故障是否违反要求。
2. 概率潮流给出源荷波动经交流关系传播后的随机状态。
3. 运行风险计算将这些状态与工程限值、严重程度函数连接。
4. 有故障概率时才进行相应加权；供电可靠性中的 EENS、动态失稳和级联风险需要额外模型，不能由本研究的静态越限指标直接代替。

Ni 等 2003 的在线风险评估已经结合故障可能性和严重程度，用于输电过载、低电压等安全监视。Li 等 2015 已对风电、负荷相关性下的输电过载建立概率与严重程度模型。因而“概率乘严重程度”“考虑相关性”属于成熟基础；本研究的价值需通过具体尾部风险计算能力及端到端比较体现。来源见第 7 节。

给定天气/源荷条件 w 及控制 u 时，可把每项风险写为对故障 c 和源荷 ξ 的联合概率加权严重程度。若 c 与 ξ 在给定 w 后仍相关，应保留条件联合分布。仅有三个故障候选时逐个报告条件风险；不能将三个条件概率相加后称全系统风险，也不能将等权故障平均称实际发生概率。数学表达归档于 derivation 第 13 节 E1–E5。

## 3. 指标顺序

| 层次 | 主输出 | 工程含义 |
| --- | --- | --- |
| 单项概率 | pV、pO，欠压/过压和不同支路端分开 | 限值被超过的可能性 |
| 单项期望严重程度 | RV=E[DV]、RO=E[DO] | 越限的概率和幅度共同作用 |
| 单项条件严重程度 | RV/pV、RO/pO | 已经越限时平均偏差，分母足够可靠才解释 |
| 并发 | pVO、JV、JO | 同时越限及其对两项期望后果的分别贡献 |
| 条件与状态组成 | P(A\|B)、P(B\|A)、四类状态 | 同时监视和条件后果信息 |
| 结构解释 | 成对传播矩阵与平方量协变 | 解释共同波动；不替代尾部概率 |

DV、DO 使用现有模型中与事件匹配的相对幅值偏差。相对严重程度并非社会经济损失，不增加任意权重总分。固定限值期望越限是主严重程度指标，首轮不强制加入 VaR、CVaR 等新模块。

## 4. 代表方法与公平比较

| 方法 | 作用与最低实现范围 | 比较边界 |
| --- | --- | --- |
| 独立完整交流 Monte Carlo | 相同输入、控制及物理事件的概率与严重程度参照，给区间和交流调用数 | 有限样本不是精确真值；零命中不代表概率零 |
| 累积量＋Cornish–Fisher 或 Gram–Charlier 重构 | 常规解析类，至少单项电压与实际端 MVA 事件 | 先检查所需矩；超范围概率、负密度或非单调 CDF 如实报告，不靠静默裁剪修饰 |
| PCE/HMS | 代理与高精度交流尾部校正相结合，匹配单项及联合事件和后果 | 优先复现 Xu2020 的必要校正思想；正式版本需精读全文和确认接口；不能只做弱化的裸 PCE |
| 已发表 TMM＋本文风险计算 | 同端低维联合事件与分项严重程度 | 给出构建、局部投影、内外层积分及更新总成本 |
| GMM 风险方法，可选 | 补充已有联合概率表示下的风险输出比较 | 不展开分量数理论或无限调参；线性 P-only 方法不能冒充完整 AC/MVA 联合参照 |
| 点估计，可选 | 均值、方差和已明确重构方式下的单项风险 | 少量矩本身不能唯一决定稀有概率，不额外强制作为主比较 |

核心为两类快速比较方法加交流参照，避免同时重做多套方法论文。缺少联合输出的累积量实现，在并发列记“未提供”，不计作错误。GMM 若扩展了原文事件或控制规则，明确称适配实现，不将其表现归于原文已报告结果。

统一输入和额定值；各方法允许使用适合自身的构建策略并计全成本，不强制所有方法采用相同训练样本数。可以共用独立验证样本降低比较方差，但验证集不能反向用于选参数、拟合或校正。代理的尾部修正交流调用计入方法成本，与最终交流参照分开。

旧方案 ν=6 是受控非高斯特例，但某些高阶输出矩可能不存在。例如自由度 ν 的非退化 t 状态，二次型的四阶矩通常需要 ν>8；直接对 MVA 等二次量使用四阶累积量时应检查相应条件。输入分布满足矩条件也不自动证明所有交流输出矩存在。主比较优先考虑 ν=12 或有物理边界的统一输入，ν=6 降为可选适用边界；截断/限幅、分布变换及由此产生的 TMM 近似均如实记录。

## 5. 尾部精度与成本验证

分别检查两个问题：

- 给定 TMM 参数时，与同一 TMM 独立采样比较，确认求积和事件处理正确。
- 从统一源荷输入计算独立完整交流参照，检查模型近似和风险计算共同产生的最终误差。不能把前一种验证当成后一种。

先报告固定工程阈值下 pV、pO、pVO、RV、RO、JV、JO，给绝对误差和参照充分时的相对误差。再用固定状态下的诊断阈值曲线检查尾部；阈值曲线是分布验证，不是运行参数灵敏度或风险边界搜索。条件概率和条件严重程度还受分母误差影响，联合命中不足时不作稳定比较。

尾概率拟从约 10⁻² 检查到 10⁻³；10⁻⁴ 是有条件的延伸，不预先承诺精度。对普通独立 Monte Carlo，小概率 p 的相对标准误差约为 1/√(Np)。p=10⁻³、N=10⁵ 时约有 100 次命中，95% 近似相对半宽约 20%；验证相对半宽约 10% 需接近 3.84×10⁵ 个样本。p=10⁻⁴、N=3×10⁵ 时约 30 次命中，相对半宽仍约 36%。这些是设计近似，正式概率采用适当区间，不对严重程度套二项公式。

沿用先一万试算、每断面十万正式参照、最多两个断面三十万的初始预算。预算不足时收缩可声明范围，或仅对少数对象采用预先约定的重要抽样并用完整交流响应验证。重要抽样需覆盖目标事件、给正确似然比及方差/重复试验；不能套普通二项区间。若采用子集模拟，应另处理相关样本误差；首轮不同时引入两套稀有事件工具。

严重程度的尾部积分对大偏差更敏感，概率积分截断质量小不代表严重程度截断误差小。对模型内积分分别声明概率和加权后果收敛要求。状态不可评估质量不静默剔除；相对于可评估状态的条件化风险须明确分母，不收敛不自动称电压崩溃。

成本列出状态/代理构建、全部交流训练/校正、参数拟合、拓扑更新、单次与多对象概率查询、严重程度积分、验证外支出和失败调用。独立参照单列，方法自身校准必须计入其总成本。同等精度比较时间，同等预算比较误差；解析性本身不保证加速。

## 6. 算例与英文 SCI 投稿

现有 IEEE39 五个计划断面和 PGLib118 三个计划断面作为主干，先用一个断面检查尾部与成本，再展开。单项对象可用同一批交流样本全量统计；并发只保留原四组同端对象。边际匹配算例降为可选，传播矩阵仅一张辅助表，不增加样本批次。

当前纯负荷压力设置能支撑通用输电风险计算。若面向 IEEE TSTE，应将有限个风电或光伏注入及源荷联合数据纳入同一评估体系，并说明预测误差、物理边界、Q 控制和平衡分摊。优先以这些输入替换部分既有随机负荷设置，不额外叠加整套场景矩阵。源荷数据集、注入位置和容量尚未选定，不在本轮编造。真实历史/预测误差数据可以提升工程可信度，但若只有合成分布，应如实称压力测试，不称真实极端天气验证。

| 期刊 | 本轮查得的范围依据 | 本项目适配判断 |
| --- | --- | --- |
| Applied Energy | 官方列能源网络、智能电网、电力系统规划运行及能源相关建模与决策 | 候选；重点证明风险尾部精度、完整成本及可再生能源相关工程用途。无需强加 OPF 才能进入该范围 |
| IEEE Transactions on Sustainable Energy | 官方强调可持续能源的发电、输配电及并网系统分析 | 候选；需要可再生能源不确定性成为实质对象，纯负荷噪声案例不够 |
| Energy | 官方当前优先跨部门综合能源系统，对纯电力研究建议 Smart Energy | 暂不默认首投，范围适配弱于前两者；不为了期刊强加综合能源模块 |
| IJEPES，补充候选 | 官方明确可靠性、概率建模及风险分析 | 若保留纯输电风险计算，方法对象更直接匹配；本轮不自动替换用户候选 |

公开大学期刊记录对前三刊列出的 2025 行均为 JCR-JIF Q1；它同时另列 SJR，两者不能混用。该页面未提供本轮 Clarivate 原始学科类别导出，因此只作公开交叉核验，正式投稿时复核当期 JCR 类别与分区。不能从放弃 TPWRS 推断 Applied Energy 或 TSTE 容易录用，也不承诺本方案一定达到其接收标准。

当前可完成的主张是风险事件与流程、风险输出精度及完整成本的可重复比较。若对照后 TMM 只在部分尾段或查询规模下更优，论文按相应范围陈述；未出现优势也保留并据实重评，不追加场景制造优越。

## 7. 来源与核查范围

PowerLit：沿用 literature-intelligence 技能，原始 JSON 目录不可读；技能附带 SQLite 索引可读，构建日期 2026-07-03，不能代表截至本日完整覆盖。定向检索 energy、applied_energy，返回 PPF/POPF 相邻记录较少；已结合在线检索核对，不能以未检出证明创新。索引中的 DOI 年与印刷年需区分，例如 122106 为 2023 DOI、2024 卷期。下面所列来源按实际可读范围使用。

- Ni 等，Online risk-based security assessment，TPWRS 18(1), 258–265, 2003：[IEEE 官方摘要](https://ieeexplore.ieee.org/document/1178806/)，[DOI](https://doi.org/10.1109/TPWRS.2002.807091)。摘要核查，作者站点 PDF 被拒绝；不声称本轮全文复读。
- Li 等，Transmission Line Overload Risk Assessment for Power Systems with Wind and Load-Power Generation Correlation，TSG 6(3), 1233–1242, 2015：[作者机构摘要及元数据](https://researchwith.stevens.edu/en/publications/transmission-line-overload-risk-assessment-for-power-systems-with/)，[DOI](https://doi.org/10.1109/TSG.2014.2387281)。概率、严重程度、相关性和点估计已有，非本文首次。
- Usaola，Probabilistic load flow with wind production uncertainty using cumulants and Cornish–Fisher expansion，IJEPES 31(9), 474–481, 2009：[出版社页面](https://www.sciencedirect.com/science/article/pii/S0142061509000258)，[DOI](https://doi.org/10.1016/j.ijepes.2009.02.003)。摘要与引言片段核查，不据此断言其所有尾部失准。
- Xu 等，Risk Assessment of Rare Events in Probabilistic Power Flow via Hybrid Multi-Surrogate Method，TSG 11(2), 1593–1603, 2020：[官方资助存档](https://par.nsf.gov/servlets/purl/10157204)，[DOI](https://doi.org/10.1109/TSG.2019.2940928)。搜索索引可读摘要；本轮全文打开超时，正式复现仍需精读全文。
- Tan 等，Computationally efficient tail distribution-aware large-scale power system overloading risk assessment，Nature Communications 17, 1521, 2026：[出版社全文](https://www.nature.com/articles/s41467-025-68241-y)。交流代理结合阈值附近交流重算；文中主要风险量为支路有功越限，不能直接当本文 MVA/电压事件基准。该文证明尾部专门校正已有研究，不能泛称现有方法均忽略尾部。
- Gao 等，Bayesian deep neural networks for spatio-temporal probabilistic optimal power flow with multi-source renewable energy，Applied Energy 353, 122106, 2024：[出版社摘要](https://www.sciencedirect.com/science/article/pii/S0306261923014708)，[DOI](https://doi.org/10.1016/j.apenergy.2023.122106)。PowerLit 记录及出版社摘要，仅用于相邻主题适配，POPF 不是本研究必须新增任务。
- [Applied Energy 官方范围](https://shop.elsevier.com/journals/applied-energy/0306-2619)；[Energy 官方范围](https://shop.elsevier.com/journals/energy/0360-5442)；[IEEE TSTE 官方范围](https://ieee-pes.org/publications/transactions-on-sustainable-energy/)；[IJEPES 官方范围](https://shop.elsevier.com/journals/international-journal-of-electrical-power-and-energy-systems/0142-0615)。ScienceDirect 的两刊范围页访问 403，改读 Elsevier 官方 Shop 描述。
- 大学公开 JCR-JIF 交叉核验：[Energy](https://www.iit.comillas.edu/publicacion/info_revista/es/219/Energy)、[Applied Energy](https://www.iit.comillas.edu/publicacion/info_revista/es/226/Applied_Energy)、[IEEE TSTE](https://www.iit.comillas.edu/publicacion/info_revista/en/352/IEEE_Transactions_on_Sustainable_Energy)。非 Clarivate 原始导出；不用于断言具体学科排名。
