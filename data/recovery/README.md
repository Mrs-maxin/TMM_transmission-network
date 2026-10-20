# 原计算成果恢复清单

本次原工作目录扫描失败：执行器返回No usable execution environment is available。原文件是否仍存在、具体路径、字节数和SHA256尚未重新核查。

原工作目录：/workspace/scratch/e880dc0ac187/TMM_Transmission。
原归档目录：/workspace/scratch/e880dc0ac187/release。
目标仓库：Mrs-maxin/TMM_transmission-network。

## 尚未迁入的原件

- Python研究源码、测试、实际运行入口和完整依赖锁。
- IEEE39/PGLib118网络、源荷数据、机组计划、限额、拓扑和配置。
- 训练/验证/校准/最终随机输入、冻结模型、全节点电压、机组输出和状态标记。
- 原始结果表、日志、成本账本、源码快照、run_manifest及图表数据。
- RELEASE_INDEX.json、原SHA256SUMS.txt、归档实物及恢复验证记录。
- 原项目MATLAB历史骨架、完整推导历史章节和此前预审/源库对照文件。

reported_archive_inventory.json只保存先前对话提到的归档角色与大小。未知文件名和SHA保留为空，不用新计算的文档SHA替代原件SHA。

## 恢复后的保存要求

1. 只在OpenAI托管执行环境读取原目录或用户重新附上的归档；先保存只读清单与SHA256。
2. 保留原包结构、依赖、输入和结果的原件，再按migration_path_map.csv确认路径。不能重写算法后称为原源码恢复。
3. data保存完整正式输入；cache只放可再生中间结果；result按断面保存数组、模型、日志、原manifest和图表。原manifest的源码SHA、种子与失败状态不得改写。
4. 数组和归档通过已配置的Git LFS上传真实内容。当前没有LFS对象或替代指针上传。保留空间、存储费用和流量账本，未知费用不填零。
5. 核对原RELEASE_INDEX与SHA清单，确认缺失、重复和排除项。不能因Git忽略规则、文件大小或缓存目录名静默遗漏。
6. 从目标仓库克隆并取得真实LFS对象，逐文件核查SHA，再空cache复算代表流程。完成后才能修改complete_migration状态。
7. 提取节点电压、相角和线路双端P/Q/MVA全分布，生成MC与冻结TMM的真实PDF/CDF、统计表及误差，不能用当前对话摘要代替。

当前待补资源是原云端执行器及目录的可读取性，或完整原归档在本任务中的可访问副本。重新授权GitHub写入不能补足原文件字节。
