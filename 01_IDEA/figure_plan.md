# IEEE14电压与线路功率分布图

版本：2026-10-20。当前没有可用数值数据，图均为计划，不绘制示意曲线冒充结果。用户当前要求先看完整潮流分布。

| figurePlanId | 读者需要判断的电力系统问题 | 对象和图形 | 比较与指标 | claimId |
| --- | --- | --- | --- | --- |
| F14-V | PQ负荷节点电压波动幅度、偏斜和曲线贴合程度 | V4/V9/V14分别给PDF与CDF | MC直方图/经验CDF对TMM半解析曲线，均值/标准差/分位数/CDF误差 | N2 |
| F14-P | 线路有功方向与波动是否一致 | 1–2的1端、4–9的9端、9–14的14端P，单位MW | 有符号P的完整分布，禁止先取绝对值 | N2 |
| F14-Q | 充电与无功交换分布是否一致 | 同三端Q，单位MVAr | 同上并保留正负 | N2 |
| F14-S | 同端实际MVA分布是否一致 | 同三端MVA | 电压乘电流计算，说明积分及PDF方法 | N2 |
| F14-PV | PV电压控制和Q限值对分布有什么影响 | PV2/3/6/8与原参考1的CDF、调压点质量、限Q部分 | 不把固定电压平滑为PDF；报告TMM设定偏离及容差带质量 | N1,N2 |

每图使用独立子图、清楚单位及颜色/线型双重区分。图上曲线必须来自真实MC和冻结TMM。机组Pg/Qg提供AC-MC统计与同输入误差；当前没有全部机组的电压TMM概率模型，不添加虚构TMM机组曲线。全14节点和20支路两端误差表不能由代表图替代。

## PDF和CDF保存

- MC直方图保存全部箱边界、计数、样本数、箱宽及区间外质量；KDE如使用需明确标注及带宽。MC经验CDF不依赖密度平滑。
- TMM曲线保存模型SHA、查询网格、半解析积分设置、随机化重复及CDF/PDF数值方法。模型抽样曲线另标，不叫半解析。
- 固定/原子输出保存位置和质量，连续密度只画连续部分并明确其质量。原子附近的容差带质量不等于严格原子。
- 网格和显示区间在最终MC前按训练/验证及冻结模型登记，不为改善视觉贴合选范围。仍输出完整分位数及区间外质量。
- 每图保存PNG、SVG或PDF、绘图CSV、figure_manifest.jsonl、命令、种子、单位、源数据/模型SHA和误差表。Python云端出图已获授权。

正式图的模板字段包括contractVersion、figurePlanId、claimId、evidenceRole、sciQuestion、physicsReproduction、metric、dataFiles、visualEncoding，以及readerQuestion、engineeringObject、visibleResult、comparedMethods、storyContribution、dataDescription。当前为physical-reproduction比较，不为填分类增加其他场景。PowerLit在本会话不可用，正式出图前仍需核查近邻原始文献；这不阻碍本次理论文档成文。

预期曲线形状和优势目前未知。MC与TMM若出现偏移、宽度不一致、多峰遗漏或PV泄漏，保留并解释；没有运行结果就不写“贴合良好”。
