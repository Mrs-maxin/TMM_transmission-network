# TMM_Transmission：含PV节点的交流概率潮流扩展及分布精度验证

**新研究基线：2026-10-20，版本restart-v1。新数值结果为零，旧工作恢复已停止。**

研究首先回答节点电压、机组无功和线路功率能否算准，再回答TMM能否准确表示这些量的完整分布。独立完整交流MC为参照；主要展示PDF、CDF、均值、标准差和分位数误差。

## 已建立的材料

- [研究目的与范围](01_IDEA/story.md)、[待验证命题](01_IDEA/claims.md)、[新计算协议](01_IDEA/derivation.md)、[统一符号](01_IDEA/symbols.md)。
- [预期算例](01_IDEA/case_plan.md)及[data/protocol/case_registry.json](data/protocol/case_registry.json)。
- [分布比较与验收标准](01_IDEA/validation_plan.md)、[图表计划](01_IDEA/figure_plan.md)和[实现接口](01_IDEA/implementation_contract.md)。
- 官方IEEE39与PGLib118原始网络文件、数值矩阵JSON及[来源SHA清单](data/networks/source_registry.json)。网络文件已重新取得，尚未做AC求解。
- [研究笔记](TMM_Transmission_note.md)、[活动论文配置](02_PAPER/config/paper-profile.yaml)及[当前状态](result/project_status.json)。

## 已做、待做与结果

| 项目 | 当前状态 |
| --- | --- |
| 研究材料、算例配置、输入来源与文件归档 | 已建立 |
| 网络数值矩阵读取、字段与编号检查 | 已静态检查，不是潮流计算 |
| 材料检查及环境记录入口 | 已编写，Python尚未运行 |
| 快速混合PV/PQ潮流、Q控制与独立交流适配 | 待重新实现和测试 |
| 普通t混合拟合、局部投影和半解析分布查询 | 待重新实现和测试 |
| MC样本、真实PDF/CDF、误差表、耗时与费用 | 尚未产生 |

当前云端执行器不可用；GitHub文件操作不是研究计算。项目不依赖旧工作目录，也不等待旧源码或数组恢复。旧状态保留在[archive/before-restart-20261020](https://github.com/Mrs-maxin/TMM_transmission-network/tree/archive/before-restart-20261020)，不参与新结果、参数选择或通过率。

## 计算与展示顺序

1. 建立并冻结满足声明要求的均值断面，检查机组Pg/Qg、节点电压和线路负载。
2. 同输入比较快速与独立完整交流，覆盖Q限值和换平衡节点。
3. 独立划分训练、验证、校准、试算和最终MC；固定普通多元t/elliptical，不做分量数扫描。
4. 全节点/全支路输出统计表；代表对象给MC直方图与TMM半解析PDF/CDF叠图。
5. 分开报告前向误差、统计模型误差、积分误差和MC不确定性。固定PV电压的点质量单独表示。

## 预期算例

先完成IEEE39常规断面；其后为区域压力、压力下5–8/6–11/14–15分别退出，以及PGLib118常规、压力、24–70退出，共八个预期概率断面。每个断面先过均值与控制检查，才执行概率计算。小规模Q控制诊断不增加概率场景。

首轮每断面1万MC试算，正式最多10万MC；本轮不自动增加至30万。不以增加场景、分量或样本改善结论。分布偏差即使超标也保留。

## Python与复算入口

TMM_Transmission.py当前只负责材料校验和环境记录，不运行潮流。以下是未来在OpenAI托管执行器上的命令，尚未执行：

```bash
python TMM_Transmission.py --check-materials
python TMM_Transmission.py --preflight
```

真实研究入口、依赖锁与执行命令在实现并通过检查后登记；不以检查入口代替概率潮流程序。候选依赖列在requirements.in，不伪称已锁定或安装。

保留模板的01_IDEA、02_PAPER、03_REFERENCE、data、cache、result和tests职责。模板MATLAB源仅保存在03_REFERENCE/template_snapshot，不作为研究实现；Python代码拟放在tmm_transmission按物理模型、统计模型、分布查询和参照交流分别组织。云端Python已授权，用户Mac不运行研究计算。

## 保存要求

正式输入、模型、随机种子、原始状态、失败输入、日志、真实run_manifest、图表数据和SHA256随阶段保存。大数组使用Git LFS，并从克隆后实际下载的对象复核，不能只上传指针。每个计算阶段完成后立即提交代码、配置及可上传结果；不等全部研究结束才备份。

目标期刊暂沿用Applied Energy，投稿证据尚未形成。普通t、PV方程及已有TMM基础不作为首次提出。

机组Pg/Qg保留同输入快算校核和完整AC-MC统计；仅电压TMM不包含随机母线负荷的联合模型，因此本轮不声称给出TMM机组Pg/Qg概率预测。
