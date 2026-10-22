# 阶段K：OP1场景组重新生成与一致性核对

计算入口：`python cloud/tests/run_stage_g_identity.py --result-stamp 20261022T125454 --op OP1`。场景组由`python cloud/tests/run_stage_d_montecarlo.py OP1 --result-stamp 20261021T202200 --workers 4`重新生成。

- 随机抽样设置与阶段一相同，SeedSequence熵2026102111。
- 2026-10-22 16:03–16:08 UTC（17:03–17:08 BST）完成，用时306 s，计算输出见`regeneration_console_output.txt`。
- 场景块位于阶段一目录`stage_d_montecarlo/OP1/`，不进入Git。被覆盖的阶段一summary和manifest已用`git checkout`恢复。

## 结果

1. 建模、定阶、厚尾检验、轻尾检验四组场景的以下各项与Git中的`stage_d_summary_OP1.json`全部一致，127008个场景全部收敛：
   - 场景数、收敛数、失败场景；
   - 平均与最大迭代次数；
   - 转换轮数分布。
2. 厚尾检验场景的逐对象统计（电压幅值，支路首末端有功、无功的均值、标准差、超额峰度、0.5%/99.5%分位数）与阶段一阶段E的`tail_statistics_OP1.csv`最大差2.9×10⁻¹³。
3. 越限事件计数与阶段一阶段E分级统计一致：
   - 任一支路负载率 > 1.00的场景27745个；
   - 任一支路负载率 > 1.30的场景30个，全部在支路94–100；
   - 另记录：任一支路 > 1.05的场景17285个，λ < −325.75 MW的场景0个。
4. 场景块SHA256记录在`stage_g_identity.json`。

## 判断

重新生成的OP1场景与阶段一原始场景一致，可用于阶段K（已由数值实验支持）。
