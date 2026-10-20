# 模板MATLAB代码说明

本目录中的.m文件来自project-template-lite模板，已经替换项目简称，保留作组织示例。它们不是从原TMM_Transmission工作目录恢复的研究源码，也不代表MATLAB本地化完成。本次没有运行这些文件。

原Python研究实现待恢复；先保留其实际包结构和依赖，再按核心潮流/统计模型、支撑流程、对比方法的职责整理。

# +TMM_Transmission_core

Put reusable model formulation, contribution-defining algorithms, solver
wrappers, and paper metrics here.

Organize implementation code into semantic subpackages instead of leaving many
loose functions in this package root. Subpackage names should describe the
physical object, mathematical component, or solution stage they own, for example
`+network_model/`, `+power_flow/`, `+optimization/`, `+risk_metrics/`, or
`+methods/`.

File names should say what the function physically or procedurally does, such as
`build_power_flow_model.m`, `solve_dispatch_subproblem.m`, or
`evaluate_voltage_risk.m`. Avoid generic names such as `helper.m`, `utils.m`,
`misc.m`, `temp.m`, or `new_method.m` for formal code.

Keep case orchestration in the root `TMM_Transmission_case*.m` files and keep
plotting, IO, cache, and Monte Carlo helpers in `+TMM_Transmission_utils/`.
