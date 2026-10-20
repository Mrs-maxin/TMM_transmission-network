# 模板MATLAB代码说明

本目录中的.m文件来自project-template-lite模板，已经替换项目简称，保留作组织示例。它们不是从原TMM_Transmission工作目录恢复的研究源码，也不代表MATLAB本地化完成。本次没有运行这些文件。

原Python研究实现待恢复；先保留其实际包结构和依赖，再按核心潮流/统计模型、支撑流程、对比方法的职责整理。

# +TMM_Transmission_sota

Put reusable comparison-method and baseline implementations here.

Organize comparison code into semantic subpackages, usually one subpackage per
algorithm family, baseline, or reproduced external method. Subpackage and file
names should identify the method or comparison procedure, for example
`+dcopf_baseline/`, `+saa_monte_carlo/`, `+linearized_power_flow/`, or
`solve_reference_baseline.m`.

Keep fairness controls, shared metrics, and data loading aligned with the
proposed-method case entries before using a comparison result as paper evidence.
