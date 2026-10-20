# 新云端任务交接

目标：含PV节点的TMM交流概率潮流扩展及完整分布精度验证。仓库为Mrs-maxin/TMM_transmission-network，当前研究study_id=restart-v1。

旧恢复已停止，不需要旧目录、旧模型或旧数组。旧数字不作证据。本轮官方IEEE39/PGLib118网络和全部协议在data，原始来源和SHA在source_registry。

先读取AGENTS.md、rule05、研究笔记、case_plan、validation_plan、implementation_contract和data/protocol/study.json。用户已授权云端Python及GitHub保存，用户Mac不做研究计算。

第一步必须在实际OpenAI托管执行器上运行真实命令，记录工具执行位置、OS、工作目录、Python/包、CPU/线程、UTC时间和可用空间。TMM_Transmission.py可做材料校验及环境记录，provider字段需要工具上下文确认，脚本不能自行证明运行商。

当前数值模块尚未实现、科学测试为零；检查入口未运行。实现顺序：独立网络/参照交流→均值断面→混合方程/快速修正/Q控制→独立输入和MC→普通TMM→局部CDF/PDF→全量统计与图表。

从C39_N开始，先生成真实节点电压和线路双端P/Q/MVA完整结果。机组Qg/Pg做同输入校核与AC-MC统计，不虚构仅电压TMM的机组概率预测。固定PV点质量显式展示；半解析查询和模型抽样区分。

八个预期断面、样本配额和种子流已登记。参数仅用训练/验证/校准确定，最终MC不用于选分量、改源荷或删样本。分布不达标照实报告。

每阶段立即保存、上传并核对源码/输入/模型/状态/失败/日志/费用/真实run_manifest与SHA；数组要有LFS实物和克隆验证。不等所有场景跑完才备份。目标期刊沿用Applied Energy，但无投稿放行结论。
