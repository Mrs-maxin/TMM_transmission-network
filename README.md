# TMM_Transmission

含PV节点的输电网交流概率潮流研究。节点电压、机组无功和线路两端P/Q/MVA是直接检查对象；统计模型采用用户确认的普通多元t混合分支。

**当前是部分迁移，原计算成果尚未完整保存到本仓库。** 本次已读取模板并恢复本对话中的研究文件；原云端执行器离线，Python研究源码、固定依赖、正式网络输入、模型、MC原始数组、历史日志和运行manifest尚未读取或上传。模板MATLAB代码只作示例，不是本项目历史源码或MATLAB本地化成果。

## 当前可查看内容

- [研究笔记](TMM_Transmission_note.md)：研究目的、潮流计算逻辑、历史报告结果和证据边界。
- [研究目的与方案](01_IDEA/story.md)、[待验证命题](01_IDEA/claims.md)、[扩展协议第15节](01_IDEA/derivation.md)和[符号表](01_IDEA/symbols.md)。
- [完整分布图表要求](01_IDEA/figure_plan.md)：节点电压及线路P/Q/MVA与MC的全范围PDF/CDF对照，当前尚未生成。
- [迁移状态](result/migration_20261020/migration_manifest.json)及[原文件恢复清单](data/recovery/README.md)。
- [对话记录中的断面计算摘要](result/migration_20261020/reported_case_summary.csv)：属于待原件复核的记录，不是本次重算结果。
- [模板来源与路径清单](03_REFERENCE/template_provenance.json)：固定模板提交及原始Git blob SHA；完整模板快照保留在03_REFERENCE/template_snapshot/。

## 项目文件职责

| 路径 | 保存内容 | 本次状态 |
| --- | --- | --- |
| 01_IDEA/ | 研究目的、命题、推导、符号、图表要求及证据索引 | 已按对话文本重建 |
| 02_PAPER/ | Applied Energy配置及论文材料 | 配置已重建；未放行投稿结论 |
| 03_REFERENCE/ | 模板、源库提交说明、用户协议基线 | 可读取部分已保存 |
| data/ | 网络、源荷、配置及完整复算输入 | 原研究输入待恢复；恢复清单已保存 |
| cache/ | 可再生中间缓存 | 未迁移；不能替代正式输入 |
| result/ | 断面结果、模型、原始数组、图表和真实运行manifest | 仅迁移记录与对话摘要；原成果待恢复 |
| tests/ | 原项目验证代码与结果 | 原测试代码待恢复 |
| +TMM_Transmission_core/、+TMM_Transmission_utils/、+TMM_Transmission_sota/ | 模板MATLAB代码分类 | 模板示例已改名前缀；未运行 |
| 03_REFERENCE/template_examples/ | 原模板33节点/123节点示例入口 | 与IEEE39/PGLib118研究结果分开 |

仓库根目录是新项目目录，项目简称为TMM_Transmission。原Python包恢复后先保留其导入和相对路径关系，确认目录映射后再整理可复用代码；当前没有虚构可运行的Python入口。

## 已记录的研究状态

历史对话报告IEEE39五个断面和PGLib118三个断面各完成10万输入的交流参照计算，并报告了有限条件下快速潮流与独立交流的一致性。原结果未恢复，本次没有重新验证这些数值。

完整电压和线路功率分布的TMM—MC叠图及误差表尚未完成。已有记录不能支持“分布贴合良好”的结论。IEEE39压力均值工况存在机组有功越限；PGLib118部分条件全部采用完整交流回退，不能声明快速计算成本优势。

## 复算和大文件

当前仓库不能从零复算，因为原研究源码和正式输入尚未迁入。依赖锁、实际运行命令和历史run_manifest须从原文件恢复后原样登记；不以模板入口或本次迁移记录替代。

已配置Git LFS追踪数组和归档格式，尚未上传LFS数据对象。只有在文件字节、SHA256、LFS对象、克隆恢复和代表流程均核查后，才能将迁移状态改为完整。详见data/recovery/README.md。

## 来源与计算位置

模板：[wuhanichina/project-template-lite](https://github.com/wuhanichina/project-template-lite)，提交861ad6f573aa64605dee9b9970b4fe3dcc0bc314，模板版本2026.08.15。

方法基础：[wuhanichina/AQPPF_TMM](https://github.com/wuhanichina/AQPPF_TMM)，固定提交0d5c999db4cbee2ea7055844d31b729d9694d747。本次未从该仓库重新取得或公开其源码。

云端Python是用户已授权路线。研究计算限定OpenAI托管执行环境；本次GitHub文件整理不构成交流潮流重算，也不代表执行器已恢复。
