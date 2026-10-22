# cloud

本目录存放云端Python基础测试的计算程序与结果（作者2026-10-21要求设立）。云端测试用于在作者电脑关机时继续计算；经检验无误的计算再改写为MATLAB，进入根目录的`+TMM_Transmission_core/`、`tests/`与`result/`主线。本目录内容不直接作为论文正式证据，需按主线流程复算或登记后才能使用。

## 目录

```text
cloud/
├── README.md
├── requirements.txt      # Python依赖
├── tmm_cloud/            # 可复用计算程序（潮流、源荷场景、概率模型、指标）
├── tasks/                # 云端任务说明，启动云端代理时输入“按 cloud/tasks/<任务>.md 执行”
├── tests/                # 云端计算入口
└── results/<任务>/<日期时间>/   # 每次计算的结果、计算记录run_manifest.json与说明README.md
```

## 工作方式

- 云端代理在自己的分支`cloud/<任务名>`上工作，每完成一个阶段即提交并推送该分支，便于及时查看；是否并入`main`由作者决定。
- 每个结果目录必须有`README.md`（计算内容、主要数值、结论与证据状态）和`run_manifest.json`（入口命令、代码提交号、Python与依赖版本、数据来源与SHA256、关键参数、随机抽样设置）。
- 单个文件不超过50 MB，不使用Git LFS。Monte Carlo场景结果（按块保存的`*_chunk*.npz`）只留在云端工作目录，不提交Git；需要时由计算入口按`run_manifest.json`记录的随机抽样设置重新生成。提交Git的结果只含README、JSON、CSV、PNG与小于5 MB的运行点、误差模型文件。
- 原始外部数据存入`data/`对应子目录并登记来源与SHA256；`cloud/`只放由其派生的结果。
- 诊断图可用matplotlib导出PNG，不作为论文正式图。
- 文字说明使用电力系统概率潮流用语，禁用词见`.cursor/rules/lexicon.md`。
