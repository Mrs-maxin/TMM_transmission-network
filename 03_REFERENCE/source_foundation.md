# TMM方法基础与继承边界

源库为用户授权只读的wuhanichina/AQPPF_TMM，固定提交0d5c999db4cbee2ea7055844d31b729d9694d747。正式已发表方法论文元数据待登记，不能从仓库名推断。

按本会话已确认的源库核查，原前向限非参考PQ且在线发电仅在参考节点；默认multiple-scaled-t，另有ordinary multivariate t / elliptical。局部普通t投影、单项电压与双端MVA事件、条件标量积分已有基础。

当前明确选ordinary t，采用用户提供的第15节扩展协议，完整新计算逻辑见../01_IDEA/derivation.md。IEEE14通过非参考P、PQ的Q、PV的V²适配，增加Qg/限值控制、参考变化和全节点重建。PV方程、普通t性质和既有局部积分不称新理论。

本轮尚未导入、移植或运行源库具体数值函数，不声称已经源码复用完成。官方IEEE14输入重新保存，MATPOWER控制代码已读取；代码阅读和数值运行分别记录。尚无本轮MC、模型、分布图或误差。
