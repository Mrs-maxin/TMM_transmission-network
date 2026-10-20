# 标准IEEE14输入

官方MATPOWER 8.1固定提交的case14原文与数值矩阵JSON，来源/SHA见../source_registry.json。未进行交流求解。原始节点编号、非参考Pg计划、Vset、Q限值、变比和充电均保留。

14节点、5机组、20支路，100 MVA；初始参考1、PV2/3/6/8。1号机组源Qg为−16.9 MVAr而Qmin/Qmax为0/10 MVAr，实际Q控制需重新验证。6/8号Vset高于通用Vmax，记录而不擅自修改。全部rateA为0，不能据此给出真实额定负载或过载风险。

network.json仅转换baseMVA、version和必需数值矩阵，不是求解结果；case14.m保留完整来源、名称与转换警告。详见../../../03_REFERENCE/ieee14_theory_sources.md。
