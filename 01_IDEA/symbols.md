# 符号、单位与实现状态

版本：2026-10-20。对应[计算理论](derivation.md)。下列代码字段均为拟定接口，尚无已运行数值实现，不把字段名当作已存在函数。

| 符号 | 电力系统含义 | 单位/维度 | 拟定字段 |
| --- | --- | --- | --- |
| n、r_m、r_0 | 节点数、当前参考、共同相角基准 | IEEE14为14；初始参考1 | n_bus, reference_bus, original_reference |
| m | PV/PQ及Q限值控制模式 | 每输入独立重置 | control_mode |
| ξ | 负荷P/Q联合输入 | MW/MVAr | PD_MW, QD_MVAr |
| Pg_plan、Vset | 非参考有功计划、调压电压设定 | MW、p.u. | Pg_plan_MW, Vset_pu |
| Pg、Qg | 机组有功、无功输出 | MW/MVAr；方程中p.u. | Pg_MW, Qg_MVAr |
| Qmin、Qmax | 真实机组无功上下限 | MVAr或统一p.u. | Qmin_MVAr, Qmax_MVAr |
| τQ、δQ | Q触发容差、近边界回退保护带 | 同一Q单位；用途不同 | Q_tolerance_pu, Q_guard_MVAr |
| B0 | 三相功率基准 | 100 MVA | baseMVA |
| Y、Yf、Yt | 完整节点导纳、支路首/末端电流映射 | 14×14、20×14复矩阵，p.u. | Ybus, Yf, Yt |
| e、f、V | 节点电压实部、虚部、复电压 | 14维，p.u. | V_complex |
| iR、iI | 节点或指定端电流实虚部 | 对应p.u. | current_real, current_imag |
| y、F_m、b_m | 非参考电压、P/Q/V²计算与目标 | 初始26维，分行功率或电压平方 | fixed_mode_coordinates, residuals |
| J_a、d1、d2 | 混合雅可比、两次逆初值修正 | 26×26、26维 | anchor_factorization, first_step, second_step |
| Vbar、v | 对齐到原1号相角的全节点电压 | 复14维、实28维，p.u. | voltage_aligned, voltage_real_full |
| v0、T_V、z | 训练中心、仿射变换、白化电压坐标 | 28维、28×d、d维 | voltage_center, voltage_transform, coordinates |
| k、πk、νk | 统计分量、权重、t自由度 | k=1..3；不同于控制模式m | weights, df |
| μk、Σk | 坐标位置与散布 | d维、d×d；不是协方差 | locations, scatter |
| U、g | 普通t的共同卡方尺度、高斯扰动 | 正标量、有效维向量 | common_chi_square, normal_vector |
| H_a、x | 全电压到同端电压/电流的投影、局部实量 | 4×28、4维 | local_projection, local_state |
| m_a,k、Λ_a,k、F | 局部位置、散布、有效生成因子 | 4维、4×4、4×r；r≤4 | local_location, local_scatter, local_factor |
| T、η | 条件积分标量、剩余联合t坐标 | 1维、r−1维；与T_V不同 | conditional_coordinate, outer_coordinates |
| F_V、u | 局部因子的电压两行、电压不变单位方向 | 2×r、r维 | voltage_factor, invariant_direction |
| P_a、Q_a、S_a | 指定支路端有功、无功、视在功率幅值 | MW、MVAr、MVA | P_end_MW, Q_end_MVAr, S_end_MVA |
| F_X、f_X | 物理量X的CDF、连续部分PDF | 概率、物理单位倒数 | CDF_grid, PDF_grid |
| ε | 不可评估输入概率质量 | [0,1]，有限MC还需区间 | unevaluable_mass |

初始维度26来自13行非参考P、9行PQ的Q、4行PV的V²；触限后行类别变化，维度不变。全节点统计坐标保留旧参考幅值；相角单位在表中采用rad，作图可转deg并标注。

功率从p.u.换回MW/MVAr/MVA乘B0。Q容差10⁻⁸ p.u.对应本例10⁻⁶ MVAr，不能把二者数值直接混用。电压原子质量是概率，不是密度高度。运行限值、数值容差和Q保护带分别登记。
