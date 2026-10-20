# IEEE14理论文档的数据与控制依据

版本：2026-10-20。本文区分源数据、代码阅读和实际数值执行；当前没有运行MATPOWER或其他交流求解器。

## 模板要求

[project-template-lite](https://github.com/wuhanichina/project-template-lite/tree/861ad6f573aa64605dee9b9970b4fe3dcc0bc314)，固定提交861ad6f573aa64605dee9b9970b4fe3dcc0bc314。已核对README、AGENTS及现有共享规则：story保存唯一研究目的；derivation记录变量/假设/公式/近似与输入；symbols保存单位/维度/字段；claims与evidence_map区分理论与数值证据；研究笔记、活动配置和changelog同步。

本轮在上述现有目录补充MD，不增加与问题无关的实验。模板MATLAB约定由用户明确授权的云端Python实现/出图覆盖；没有完成MATLAB本地化。

## 标准IEEE14数据

- 来源：[MATPOWER 8.1 case14.m](https://github.com/MATPOWER/matpower/blob/1a828c7af590714499284e36ee9c81273388c594/data/case14.m)。
- 固定提交：1a828c7af590714499284e36ee9c81273388c594；Git blob：a67682c6654077a0dd71f0eab936de2ec830903d。
- 原文SHA256：2ffc4e1b734ae6c5e92dbe68b4e36010ed695a4bbcc4d065c74c4fbc39fcf3c1。
- 原文位置：data/networks/ieee14/case14.m；字面数值转换：data/networks/ieee14/network.json。外部节点编号和所有必需物理列保留；原文还保留母线名称及转换警告。
- 已检查条目为14节点、5机组、20支路、100 MVA；参考1、PV2/3/6/8。这里是文件读取与字段核查，不是本项目计算结果。

原文转换注释说明1号机组原Qmax=Qmin时将Qmax改为Qmin+10，并说明各支路MVA限值缺失后置0。本轮保留该源版本，不把其Q限额或零额定值称为完整设备实测能力。1号Qg=−16.9与限值0/10 MVAr不一致、6/8号Vset高于通用Vmax，均须在真实复算中解释，不能通过放宽限额消除。

MATPOWER软件许可证文件已保存，但其软件授权声明不覆盖外部case数据集；不把case14重新标注为本项目原创或由软件许可证自动授予。本项目保留原文件全部出处与说明。

## Q限值和独立参照依据

已读取相同固定提交的下列官方代码：

- [lib/runpf.m](https://github.com/MATPOWER/matpower/blob/1a828c7af590714499284e36ee9c81273388c594/lib/runpf.m)：legacy交流路径的Q限值判断、同轮同时/逐个处理、PV转PQ、参考更新及相角还原。
- [lib/bustypes.m](https://github.com/MATPOWER/matpower/blob/1a828c7af590714499284e36ee9c81273388c594/lib/bustypes.m)：内部编号与新参考选取。没有PV可选时需由适配层明确记录无可用平衡状态，不能把索引错误当成物理崩溃。
- [lib/newtonpf.m](https://github.com/MATPOWER/matpower/blob/1a828c7af590714499284e36ee9c81273388c594/lib/newtonpf.m)：独立极坐标功率平衡Newton的结构依据。

runpf读取原始QG/QMAX/QMIN与opf.violation；必须核查MVAr单位并与本项目p.u.容差转换一致。计划控制采用同时切换、限Q后本次不恢复PV、首个剩余PV接任参考。实际后端若用PYPOWER及项目控制器，必须记录适配范围和差别；代码阅读尚不能证明数值行为一致。测试前冻结版本和控制规则。

## TMM基础与本文推导边界

用户授权参考私有[AQPPF_TMM](https://github.com/wuhanichina/AQPPF_TMM/tree/0d5c999db4cbee2ea7055844d31b729d9694d747)，固定提交0d5c999db4cbee2ea7055844d31b729d9694d747。按本会话已确认核查，原前向仅非参考PQ且发电限参考，默认multiple-scaled-t，另有ordinary multivariate t / elliptical。局部投影、单项电压/MVA事件和条件积分已有基础。

本轮明确选ordinary t，沿用用户提供的第15节扩展协议及已确认模型关系，补写IEEE14计算链、普通t条件积分、MVA适用方向和MC误差协议。这里没有重新导入、移植或运行源库具体函数，不声称新实现已经复用其源码。正式已发表TMM论文元数据仍待登记，不从仓库名猜测引用。

混合方程、雅可比、t仿射/条件公式是计算依据，不作为新发现。PV电压支撑近似、全局混合控制模式、数值逆初值和积分均须用真实计算验证。文件可核查不等于方法已验证。
