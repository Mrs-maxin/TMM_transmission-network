# 待验证命题与理论依据

版本：2026-10-20，IEEE14基础测试。理论文件已经成文；数值命题均未验证。已有TMM、PV方程与普通t性质不认领为首次提出。

| claimId | 待验证命题 | 工程检查与证据位置 | evidenceState |
| --- | --- | --- | --- |
| N0 | 混合PV/PQ快速潮流可在相同负荷下得到与独立完整交流一致的电压和支路功率 | 26维P/Q/V²、Qg、参考Pg、20支路双端；拟用E14-PHYS | unverified |
| N1 | Q限值与参考切换后仍保留正确电压、机组输出与模式 | 单上/下、多PV、换参考、无平衡能力；拟用E14-QCTRL | unverified |
| N2 | 普通TMM能够在所测IEEE14负荷条件下描述电压及线路功率完整分布，并给出可复核误差 | PDF/CDF、均值/标准差/分位数、MC区间、PV泄漏；拟用E14-DIST | unverified |
| N3 | 所测计算流程可以保存完整输入与状态，并从代码/data复算，记录全部成本 | 空cache与实际归档恢复；拟用E14-REPRO | unverified |

混合方程维度、固定模式二次关系、普通t仿射投影和条件分布具有理论依据，见derivation.md；这不构成N0–N3数值通过证据。首版可能有大量回退、固定PV泄漏或分布误差超标，均需保留。

统一绑定：primaryInnovationAxis=technical-object；noveltyMode=one-to-hundred；technicalObjectType=model-formulation；manuscriptUse=conditional；lifecycleStage=zero-to-one-explore。每项evidenceBoundary均为“待实现及独立IEEE14数值验证；当前无新结果”。上表evidenceId为预留编号，不是已存在run。

同输入前向、模型内积分与独立交流分布分别验证。仅电压TMM不声称获得带随机母线负荷的全部机组Pg/Qg PDF。不得以最终MC选分量、改输入、删不利样本或把无效输入记为安全。
