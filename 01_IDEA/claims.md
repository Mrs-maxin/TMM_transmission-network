# 待验证研究命题

更新：2026-10-20。研究目的与2026-10-19用户基线一致，原基线保存于03_REFERENCE/user_protocol_baseline_20261019/。本次只迁移模板和可从对话恢复的材料，原研究代码和数组未取得。

| 编号 | 工程命题 | 当前记录 | 当前证据状态 |
| --- | --- | --- | --- |
| C0 | 混合PV/PQ潮流能准确计算电压、Qg、平衡机组Pg及双端潮流，并正确执行Q限值与换参考 | 历史报告有IEEE39的40组和118的20组同输入检查；原文件待恢复 | 结果待复核 |
| C1 | 普通TMM能准确表示节点电压和同端实际MVA的单项/并发概率及分项后果 | 已报告偏差；完整PDF和全分布误差未完成 | 当前证据不足 |
| C2 | 同时越限状态和分项后果补充单项概率信息 | 定义与协议保留，原数值证据未迁入 | 尚未验证 |
| C3 | 在匹配输入、机组控制和设备限值下，能确定与代表方法比较的精度和完整成本 | 累积量、PCE/HMS适配范围和费用记录仍有缺口 | 当前证据不足 |
| C4 | 矩存在时，传播矩阵可辅助解释电压与线路电流的共同波动 | 原传播表和模型未取得 | 尚未验证 |

当前优先检查完整电压、相角和线路P/Q/MVA分布，再讨论尾部。有限同输入误差只能支持潮流求解检查；模型内抽样只能核查积分；最终分布精度需要独立完整交流参照。

TMM已作为已有方法基础，PV方程和普通t局部性质不重新认领为首次。普通t分支已确认，不切换默认多尺度t。无权重的拓扑条件结果不汇总年度风险。

绑定：primaryInnovationAxis=technical-object，technicalObjectType=model-formulation，noveltyMode=one-to-hundred，caseContract=new-variable-scenario，lifecycleStage=zero-to-one-explore，manuscriptUse=conditional。正式数值evidenceId和历史run_manifest待原文件恢复后登记。
