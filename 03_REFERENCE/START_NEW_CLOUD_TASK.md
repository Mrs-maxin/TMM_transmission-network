# 新执行器交接：IEEE14正常拓扑

当前study_id=ieee14-baseline-v1，目标是含PV节点的普通TMM交流概率潮流扩展及完整分布精度验证。只用标准IEEE14正常拓扑与相关负荷P/Q波动，旧八算例停用，旧目录恢复已停止，旧数字不作新证据。

先读AGENTS.md、rule05、研究笔记、01_IDEA/derivation.md及case_plan/validation_plan/implementation_contract，读取data/protocol与data/networks/source_registry.json。理论已成文，数值模块待实现，科学测试/MC/拟合/图为零。Python材料入口未运行。

第一项真实命令须确认工具位于OpenAI托管执行器并记录OS、目录、UTC、Python/包版本、CPU/线程与空间。脚本不能自行证明运行商；不能在用户Mac、GitHub CI或Sites运行后称为OpenAI研究计算。当前离线状态不能靠文档提交解除。

先复算原case14，审核源文件1号Qg与限值及参考切换；保留原非参考Pg、Vset和Q限值，不重新优化。所有rateA为0，不给过载风险。再实现26维混合方程、固定模式对照、完整Q控制和全节点重建，之后训练普通TMM与独立MC。

代表图为V4/V9/V14及1–2的1端、4–9的9端、9–14的14端P/Q/MVA；全表14节点、5机组、20支路双端。机组仅同输入校核与AC-MC统计，不声称电压TMM获得全部机组PDF。PV点质量、幅值泄漏和失败质量单列。

每阶段在云端保存完整检查点；用户明确“Push／推送”时将当时测试与结果整体保存GitHub并同步相关文档、核查实际数组/SHA。直接编辑文档请求只授权该次范围。不得恢复旧自动推送政策。

没有数值证据就不报告精度、曲线贴合、回退率、加速或完成；失败和未覆盖同样保存。依赖锁、真实入口和后端版本须在实际安装/运行后登记。
