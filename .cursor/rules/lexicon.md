# Writing Lexicon

This rule file defines wording checks for AI-generated research notes,
manuscript drafts, README text, code comments, review responses, and figure
captions in this template.

## Preferred Style

- State the technical conclusion directly.
- Write innovation points as positive statements about the minimum identifiable
  research object, concrete engineering object, contribution-appropriate
  relation, theoretical significance, engineering value, mathematical model,
  engineering setting, measurable result, and evidence source. The relation
  may be physical causality, model mapping, computational property, control
  response, validity condition, observability, or an engineering decision link.
- Use complete declarative sentences.
- Start with the system object, model term, constraint, data source, metric, or
  numerical result.
- Name the contribution-appropriate relation, affected object, case set,
  metric, and evidence state. Use "physical mechanism" only when the evidence
  establishes a physical causal chain and a testable prediction.
- In Chinese power-system prose, name the physical object first: bus, branch,
  feeder, source/load process, operating scenario, voltage, power flow, current,
  constraint, probability distribution, violation event, or recovery action.
- In `TMM_Transmission_note.md`, `01_IDEA/story.md`, `01_IDEA/claims.md`, and
  innovation anchors, use field-recognized terms, engineering objects, physical
  quantities, model relations, and verifiable case names in body text.
- Use plain technical verbs such as compute, estimate, solve, bound, compare,
  validate, export, register, update, and test.
- Attach an intensifier to a number, tolerance, scale, or explicit evidence
  source.

## Rewrite Checks

Before releasing AI-generated text, scan for these patterns:

- The technical statement is missing the minimum identifiable research object,
  model relation, measurable result, or evidence source. Add those elements
  directly.
- Human-facing notes, story files, and manuscript prose should organize
  innovation points as minimum identifiable research object, theoretical
  significance, engineering value, model, result, engineering setting, evidence
  source, and next experiment.
- Negative buildup and fragmentary emphasis. Replace them with one complete
  sentence.
- Rhetorical question setup. Replace it with the technical condition or result.
- Prompting phrases that introduce an explanation. Start with the object or
  formulation.
- Meta commentary about the section, future explanation, or authorial process.
- Business jargon and promotional verbs.
- Over-mathematized or programming-shaped carriers. If the reader-facing point
  is about a power-system object, rewrite the sentence around the physical
  object and its equation, constraint, state, or metric.
- Vague consequence statements. Name the concrete implication.
- Scope-free extreme terms. State the case set, time window, model class, solver,
  dataset, or assumption range.
- Hard-banned over-mathematized or programming-shaped labels. If generated text
  contains any token listed in `Hard Banned Terms`, rewrite it before release
  around the power-system object, equation, constraint, state, metric,
  calculation procedure, or model parameter.
- Over-mathematized terms listed in `Use With Care Terms`. Prefer a
  power-system object, equation, constraint, state, metric, calculation
  procedure, or model parameter unless the term is mathematically necessary and
  defined in the derivation.
- Coined method names, internal experiment configuration aliases, temporary
  parameter-group names, branch names, or script abbreviations in body text.
  Move them to evidence ids, entry-file fields, result paths, run manifests,
  Git records, or configuration fields, and rewrite the sentence around the
  engineering object, physical relation, model name, metric, or case name.
- Em dashes. Use commas or periods.
- Empty intensifiers. Delete them or attach a number.

## Use With Care Terms

The following terms are not hard banned, but they should usually be rewritten in
reader-facing prose because they often make a power-system claim sound more
abstract than the engineering object requires.

| Type | Terms | Prefer |
| --- | --- | --- |
| Over-mathematized Chinese labels | 算子, 可微 | equations, mappings, constraints, states, physical relations, model parameters, or calculation steps |
| Over-mathematized English labels | operator, differentiable | equations, mappings, constraints, states, physical relations, model parameters, or calculation steps |

## Hard Banned Terms

The following terms may appear only in this lexicon or in a quoted source being
criticized. They must not appear in `TMM_Transmission_note.md`, `01_IDEA/story.md`,
`01_IDEA/claims.md`, innovation anchors, introductions, abstracts,
contribution lists, captions, or reader-facing README prose.

| Type | Forbidden tokens | Rewrite around |
| --- | --- | --- |
| Over-mathematized Chinese labels | 证书, 零空间 | power-system equations, constraints, states, feasible checks, measurements, physical relations, or model parameters |
| Programming-shaped Chinese labels | 范式 | methods, models, calculation procedures, input-output mappings, metrics, or engineering decisions |
| Over-mathematized English labels | certificate, nullspace | power-system equations, constraints, states, feasible checks, measurements, physical relations, or model parameters |
| Programming-shaped English labels | paradigm, estimator | methods, models, calculation procedures, input-output mappings, metrics, or engineering decisions |
| Machine-learning and software-shaped Chinese labels (author requirement, 2026-10-21) | 训练, 训练集, 训练样本, 训练数据, 验证集, 验证样本, 测试集, 测试样本, 标签, 种子, 对照, 回退, 锚点, 过拟合, 超参数, 进程, 后台, 日志, 脚本 | the probabilistic power-flow terms in the table below |

Probabilistic power-flow wording that replaces the machine-learning and
software-shaped labels in reader-facing prose and in conversation with the
author:

| Avoid | Write |
| --- | --- |
| 训练样本 / 训练数据 | 建模场景：用于建立TMM的源荷场景及其节点电压 |
| 验证样本 / 验证集 | 定阶场景：用于确定分量数的另一组独立源荷场景 |
| 测试样本 / 测试集 | 检验场景：独立交流Monte Carlo场景 |
| 种子 | 随机抽样设置，或直接写“同一组源荷场景” |
| 标签 | 所属分量、所属运行方式 |
| 对照 | 比较、参照（交流Monte Carlo参照） |
| 回退 | 改用逐点牛顿-拉夫逊潮流计算 |
| 锚点 | 展开点（建立快速映射所用的牛顿-拉夫逊潮流解） |
| 选阶、选K | 确定分量数 |
| 接收/拒绝拟合 | 满足/不满足分量有效性条件 |
| 过拟合 | 分量数过多使定阶场景对数似然下降 |
| 进程、后台、日志、脚本 | 计算、计算记录、计算入口文件 |

## Avoid Terms

Use these tokens only when quoting source material or documenting a writing rule:

- navigate
- unpack
- landscape
- deep dive
- game-changer
- at its core
- it is worth noting
- it turns out
- as we will see
- in this section
- when it comes to
- think about it
- here is what I mean
- significant
- structural
- fundamentally
- inherently
- simply
- really
- every
- always
- never
- all
