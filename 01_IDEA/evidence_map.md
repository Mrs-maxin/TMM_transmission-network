# 当前理论材料与待产生的数值证据

版本：2026-10-20。已有理论文档和IEEE14源数据；本轮科学测试、MC状态、TMM模型和PDF/CDF均为零。

| evidenceId | 对应内容 | 文件/预留结果位置 | 状态 |
| --- | --- | --- | --- |
| DOC14-THEORY | 全链计算理论、近似与适用条件 | 01_IDEA/derivation.md；03_REFERENCE/ieee14_theory_sources.md | 文档已成文，不是数值通过 |
| SRC14-NET | 官方case14原文与字面JSON转换、来源SHA | data/networks/ieee14；data/networks/source_registry.json | 已取得，AC未运行 |
| E14-PHYS | 方程/导数/同输入潮流 | 拟用result/C14_N/physics | unverified；尚无run |
| E14-QCTRL | Q限值、参考切换与保护审计 | 拟用result/C14_N/control | unverified；尚无run |
| E14-DIST | 冻结TMM与独立AC-MC完整分布 | 拟用result/C14_N/distributions及figures | unverified；尚无run |
| E14-REPRO | 真实环境、空cache复算及实际归档恢复 | 拟用result/C14_N/reproduction | unverified；尚无run |

配置、材料校验、GitHub写入/读回SHA不是潮流或分布验证。result/materials中重建及旧提交校验记录仅证明各自记载版本的文件，不证明当前方法通过。

新实现没有复用完成的源码，没有独立交流实际后端、依赖锁或速度结论。固定AQPPF提交只作已知方法依据，正式论文引文待登记。旧成果留存分支不进入当前证据。
