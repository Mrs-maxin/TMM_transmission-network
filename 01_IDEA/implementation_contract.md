# 数值实现与数据接口

版本ieee14-baseline-v1。以下模块为重新实现的约定，当前尚未实现；TMM_Transmission.py仅是材料和环境检查入口。

| 拟采用模块 | 输入/输出 | 必要核查 |
| --- | --- | --- |
| network/admittance | 原始外部编号网络→内部映射、Ybus/Yf/Yt | 独立后端、变比/相移/充电、正常拓扑与编号 |
| power_flow/mixed_equations | 当前模式、复电压、指定输入→P/Q/V²与雅可比 | 维度、有限差分、全PQ退化 |
| power_flow/fast_fixed_mode | 冻结展开点、样本→全电压、残差、修正次数 | 同输入独立AC及不可靠回退 |
| reference/pypower_adapter | 同一网络与输入→独立极坐标AC状态 | 实际包版本、原始branch/gen输出核查 |
| power_flow/q_control | 初始状态、Q策略→最终模式/参考/限值状态 | 单上/下、多PV、换参考及无能力 |
| statistics/voltage_coordinates | 对齐全电压→训练中心/变换/坐标 | 旧参考幅值、弱差压和重建误差 |
| statistics/ordinary_t_mixture | 训练/验证坐标→冻结普通t混合 | 共同尺度、散布、df和正则化 |
| distributions/local_projection | 冻结模型、支路端→有效局部因子 | 秩、确定方向、真实同端MVA |
| distributions/conditional_t_cdf | 标量事件、局部参数→CDF与积分诊断 | 已知极限、同模型抽样、收敛 |
| outputs/whole_distribution | 原始状态/冻结模型→全量统计及PDF/CDF | 单位、端点、点质量和误差 |
| workflow/run_case | 注册算例→阶段产物与真实manifest | 输入独立、失败质量、费用和种子 |
| workflow/checkpoint | 当前阶段→云端检查点；用户Push后远端实物和SHA核验 | 数组完整、空cache复算、克隆核对 |

模块位于tmm_transmission/，按物理/统计/分布/参照职责分开；上表路径相对于该包。当前只登记接口，不创建冒充实现的空求解函数。

状态对象须保留input_id、全部节点电压、逐机组Pg/Qg、最终节点类型、Q触限方向、参考节点、P/Q/V²残差、求解/回退状态、内外完整交流调用和错误原因。未评估行保留，用状态及显式缺失值表示，不用零电压或零功率冒充成功状态。

原始输入、状态、模型与全部表/图字段见data/protocol/output_contract.json。固定电压和固定机组功率的原子质量标记随表保存。

快速网络实现不能只复用参照的同一套有缺陷导纳然后互相印证。允许共享原始输入和明确控制规则，但混合方程、网络输出和控制诊断至少各有独立来源核查。

参考后端若不可得，先记录阻塞，不用快速求解器自采样冒充独立完整交流。没有实际运行MATPOWER就不写“MATPOWER参照”。正式代码实现之后再更新模块状态和可复算命令。

当前只有IEEE14正常拓扑。控制器与实际参照后端尚未冻结，所有数值模块仍待实现。理论依据与近似界限见derivation.md。不得自动运行旧的压力、停运或IEEE39/118算例。
