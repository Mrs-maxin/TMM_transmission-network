# 正式输入

当前唯一活动网络为networks/ieee14/case14.m与network.json，来源和SHA见networks/source_registry.json。它们是输入，不是本轮潮流结果。

protocol保存IEEE14正常拓扑、负荷P/Q、独立样本流、固定3分量ordinary-t、Q控制、验收与输出协议。所有正式计算从代码与data出发，cache不得成为唯一输入。

已有IEEE39/PGLib118数据保留为inactive_reference_networks，不在当前概率算例中运行。case14线路rateA均为0，额定负载及过载指标不可用；不得补造额定值。真实随机输入、全部尝试样本和失败记录待执行后保存。
