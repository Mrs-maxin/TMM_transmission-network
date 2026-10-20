# 新研究符号与接口

版本restart-v1。代码映射为拟采用字段，数值模块尚未实现；不把命名登记当成实现证据。

| 工程量 | 符号/含义 | 单位 | 拟采用字段 |
| --- | --- | --- | --- |
| 节点与参考 | i、原相角基准r0、当前参考r | 外部母线号 | bus_id, original_reference_bus, reference_bus |
| 节点电压 | V=e+jf、幅值Vm、对齐相角Va | p.u.、rad | voltage_complex, Vm_pu, Va_rad_aligned |
| 负荷与新能源 | PD/QD、PR/QR | MW、MVAr | PD_MW, QD_MVAr, RES_P_MW, RES_Q_MVAr |
| 同步机计划与输出 | Pg_plan、实际Pg/Qg | MW、MVAr | Pg_plan_MW, Pg_MW, Qg_MVAr |
| 调压与无功限额 | Vset、QMIN/QMAX | p.u.、MVAr | Vset_pu, Qmin_MVAr, Qmax_MVAr |
| 节点导纳与支路导纳 | Ybus、Yf、Yt | p.u. | Ybus, Yfrom, Yto |
| 线路两端功率 | Pf/Qf/Sf、Pt/Qt/St | MW、MVAr、MVA | P_from_MW, Q_from_MVAr, S_from_MVA, P_to_MW, Q_to_MVAr, S_to_MVA |
| 混合方程与目标 | F、b | P/Q为p.u.，V²为p.u.² | mixed_equations, mixed_targets |
| 全节点电压坐标 | v、v0、TV、z | p.u.及白化坐标 | voltage_center, voltage_transform, coordinates |
| 普通t混合 | 权重π、位置μ、散布Σ、自由度ν | 对应坐标单位 | weights, location, scatter, df, family=elliptical |
| 局部端量 | x=[e,f,iR,iI] | p.u. | local_location, local_scatter, local_factor |
| 无功模式 | 调压、上限、下限和当前参考 | 状态标签 | bus_type_final, Q_limit_direction, reference_bus |
| 输入与评估状态 | 源荷样本、可评估/失败 | input_id、状态 | input_id, status, failure_reason |
| 分布误差 | 均值/标准差差、分位数差、CDF差 | 原单位或概率 | mean_abs_error, std_abs_error, quantile_abs_error, CDF_grid_max_error |

散布矩阵不能当协方差；角度统一为rad，外部源文件Va的度数在适配时显式转换；baseMVA只用于功率标幺换算，不给电压幅值乘基准。线路正号表示从该端节点注入支路。

完整状态字段与文件规范见data/protocol/output_contract.json；正式代码实现后补充真实模块路径和行级对应。
