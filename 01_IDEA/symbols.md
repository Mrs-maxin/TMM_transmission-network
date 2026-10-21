# 符号表

版本2026-10-21，与`derivation.md`重写版一致。重写前的符号见`derivation_archive_20261021.md`附录。

## 1 潮流模型（derivation第1、2节）

| 符号 | 含义 | 单位/维度 | 代码变量 |
| --- | --- | --- | --- |
| $n$，$s$ | 节点数、平衡节点 | — | `control.busCount`、`control.slack` |
| $\mathcal N^{*}$，$\mathcal N_{\rm PV}$，$\mathcal N_{\rm PQ}$ | 非平衡、PV、PQ节点集合 | 节点编号 | `control.nonReference`、`control.pv`、`control.pq` |
| $e_i,f_i$ | 节点电压实部、虚部 | p.u. | `real(voltage)`、`imag(voltage)` |
| $G_{ij},B_{ij}$ | 节点导纳矩阵实部、虚部 | p.u. | `anchor.Ybus` |
| $P_{{\rm G},i},Q_{{\rm G},i}$ | 机组计划有功、无功出力（固定于限值时为限值） | MW、Mvar | `control.scheduledP`、`control.scheduledQ` |
| $P_{{\rm L},i},Q_{{\rm L},i}$ | 负荷有功、无功 | MW、Mvar | `loadP`、`loadQ` |
| $U_i^{\rm set}$ | PV节点及平衡节点电压设定值 | p.u. | `control.voltageSet` |
| $\boldsymbol x$ | 非平衡节点电压实部、虚部，实部在前 | $2(n-1)$维 | `state` |
| $\boldsymbol F,\boldsymbol y$ | 潮流方程与给定量（P、PQ节点Q、PV节点电压平方） | $2(n-1)$维 | `evaluate_mixed_power`、`build_mixed_power_targets` |
| $\boldsymbol J_0$ | 展开点雅可比矩阵 | $2(n-1)$阶 | `anchor.jacobian`，逆为`anchor.inverseMatrix` |
| $\boldsymbol h$ | 潮流方程二次项 | $2(n-1)$维 | 由`evaluate_mixed_power`差分得到 |
| $\alpha_m$ | 修正步长 | — | `refine_mixed_jacobian`中`step` |
| $\varepsilon_y,\varepsilon_x$ | 功率残差、状态残差估计容许值 | $10^{-8}$、$10^{-7}$ p.u. | `PowerTolerance`、`StateTolerance` |
| $\sigma_g$ | 机组运行状态（0调压，1上限，−1下限） | — | `modeState` |
| $\varepsilon_Q$ | 无功越限判别容差 | $5\times10^{-6}$ Mvar | `ReactiveTolerance` |

## 2 电压联合概率模型（derivation第3节）

| 符号 | 含义 | 单位/维度 | 代码变量 |
| --- | --- | --- | --- |
| $\bar{\boldsymbol x},\boldsymbol E_d,\boldsymbol D$ | 建模场景电压均值、保留特征向量、方向标准差 | p.u.；$2(n-1)\times d$；$d$阶 | `transform.center`、`transform.basis`、`transform.scale` |
| $\boldsymbol z$，$d$ | 标准化电压坐标及维数 | $d$维 | `x`、`transform.rank` |
| $K$，$\pi_k$ | 分量数、分量权重 | — | `componentCounts`、`model.weights` |
| $\boldsymbol\mu_k,\boldsymbol\Sigma_k,\nu_k$ | 位置向量、散布矩阵、自由度 | $d$维、$d$阶、— | `model.means`、`model.scatter`、`model.df` |
| $S_k$ | 共同波动强度 | — | 抽样中的`strength` |
| $\boldsymbol H,\boldsymbol L_{\mathcal W}$ | 共同强度驱动的电压变化方向、注入选择矩阵 | $2(n-1)\times n_w$ | `determine_common_strength_directions`中`sensitivity` |
| $r$，$\boldsymbol Q_u,\boldsymbol Q_w$ | 强度作用方向数及两块方向 | $d\times r$、$d\times(d-r)$ | `directions.rank`、`directions.Q1`、`directions.Q2` |
| $\boldsymbol m_k,\boldsymbol\Sigma_k^{u}$ | 分块TMM强度块位置与散布 | $r$维、$r$阶 | 分块模型`means`、`scatter` |
| $\boldsymbol a_k,\boldsymbol G_k,\boldsymbol\Omega_k$ | 其余块截距、回归系数、残差协方差 | $(d-r)$维等 | 分块模型`intercept`、`coefficient`、`residual` |
| $\gamma_{ik},\tau_{ik},\delta_{ik}$ | 后验概率、强度权重、t距离 | — | `membership`、`tau`、`distance` |
| $N$，$N_{\rm v}$ | 建模场景数、定阶场景数 | 5256、1752 | `modelingCount`、`determinationCount` |

## 3 电气量与评价指标（derivation第4、6节）

| 符号 | 含义 | 单位 | 代码变量 |
| --- | --- | --- | --- |
| $U_i,\theta_i$ | 节点电压幅值、相角（以平衡节点为基准） | p.u.、(°) | `voltage`、`voltageAngle` |
| $P_{ij},Q_{ij}$ | 支路首端（末端同理）有功、无功 | MW、Mvar | `linePfrom`、`lineQfrom`（`linePto`、`lineQto`） |
| $S_{\rm B}$ | 功率基准值 | 100 MVA | `mpc.baseMVA` |
| $N^{\rm s}$ | 分布积分点数 | 50000 | `modelDrawCount` |
| $c,\nu$ | 负荷变异系数、共同强度自由度 | 5%、12 | `CoefficientOfVariation`、`LatentDf` |
| $e_\mu,e_\sigma,e_p,D_{\rm KS}$ | 均值误差、标准差相对误差、分位误差、KS距离 | 与电气量同单位 | `meanAbsoluteError`、由方差计算、`q01/q99`差、`KS` |
