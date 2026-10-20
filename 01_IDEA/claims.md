# 待验证命题

版本restart-v1。以下命题均尚未完成新数值验证。主研究对象是含PV节点的交流概率潮流，不以风险尾部计算作为完成条件。

| 编号 | 工程命题 | 所需证据 | 当前状态 |
| --- | --- | --- | --- |
| N0 | 固定模式快速求解能与独立完整AC一致计算电压、相角、Qg、平衡Pg和线路两端P/Q/MVA | 同输入误差与混合方程残差，涵盖导纳、变比/相移/充电和索引 | 尚未验证 |
| N1 | Q限值与换平衡节点能给出一致的最终节点类型、机组限值等号及全节点状态 | 单上/下限、多PV、旧参考触限、能力不足及非孤岛退出诊断 | 尚未验证 |
| N2 | 普通TMM能在声明误差内表示节点电压和线路功率的完整分布 | 独立MC与TMM的PDF/CDF、均值/标准差/分位数误差，PV点质量/泄漏 | 尚未验证 |
| N3 | 新流程在压力和选定线路退出及118节点系统中具有可说明的适用范围和完整成本 | 同协议的迁移结果、回退率、失败输入、耗时和交流调用 | 尚未验证 |

N0通过不能自动放行N2；模型内采样只检验TMM概率计算，独立源荷AC才检验实际分布。固定电压的点质量与连续密度误差分开，不借归零或投影隐藏偏差。

绑定：primaryInnovationAxis=technical-object；technicalObjectType=model-formulation；noveltyMode=one-to-hundred；caseContract=new-variable-scenario；lifecycleStage=zero-to-one-explore；manuscriptUse=conditional。新evidenceId只有实际运行后才登记。
