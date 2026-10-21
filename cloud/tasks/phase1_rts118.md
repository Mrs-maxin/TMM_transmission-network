# 云端任务：阶段一 RTS-GMLC 厚尾输入与 IEEE118 交流 Monte Carlo

启动方式：在 Cursor 云端代理中选择仓库 `Mrs-maxin/TMM_transmission-network`、基准分支 `main`，输入“按 `cloud/tasks/phase1_rts118.md` 执行”。以下为给云端代理的完整说明。

---

You are working in the GitHub repository Mrs-maxin/TMM_transmission-network (branch main is the base). It is a power-system research project: extending a published TMM (multivariate t mixture model) AC probabilistic power flow method to transmission networks with PV buses, and later assessing voltage/line violation risk. The main computation line is MATLAB (root packages +TMM_Transmission_core, tests/, result/) and runs only on the author's laptop. Your job is to run Python "basic tests" in the cloud so that work continues while the laptop is off.

READ FIRST: AGENTS.md, .cursor/rules/00-lightweight-research-coding.mdc (section "Git Sync With GitHub"), .cursor/rules/lexicon.md, cloud/README.md, 01_IDEA/derivation.md (method text, Chinese), TMM_Transmission_note.md section "当前基准". Follow them.

## Hard rules

- Work only inside cloud/ (Python package cloud/tmm_cloud/, entries cloud/tests/, results cloud/results/) and data/rts_gmlc/ (raw external data). Do not modify the MATLAB main line (+TMM_Transmission_core, +TMM_Transmission_utils, tests/, result/) or 01_IDEA/*.
- Work on a branch named cloud/phase1-rts118 (or the branch your environment gives you). After EACH completed stage below, commit and push that branch so results are visible promptly. Never push to main, never merge.
- Every result folder cloud/results/phase1_rts118/<YYYYmmddTHHMMSS>/ (one folder per stage or one shared folder with stage subfolders) must contain README.md (Chinese, power-system wording: what was computed, key numbers, judgment with evidence state such as 已由数值实验支持 / 尚未验证) and run_manifest.json (command, git commit, Python and package versions, data sources with SHA256, key parameters, random sampling settings).
- Single files must be < 50 MB. Store large Monte Carlo arrays as chunked compressed .npz (float32 allowed for voltages and powers of the large check set). No Git LFS.
- In all Chinese text use power-system probabilistic power flow terms. FORBIDDEN words in prose: 训练, 训练集, 训练样本, 验证集, 测试集, 测试样本, 标签, 种子, 对照, 回退, 锚点, 过拟合, 超参数. Use instead: 建模场景 (scenarios used to build a model), 定阶场景 (independent scenarios used to choose the number of mixture components), 检验场景 (independent AC Monte Carlo scenarios), 随机抽样设置 (random-number settings), 比较/参照, 改用逐点牛顿-拉夫逊潮流计算, 展开点.
- Diagnostic figures may be matplotlib PNGs (not formal paper figures).
- Do not silently drop failed power flows; record counts and reasons.

## Stage A: RTS-GMLC data

1. Download from https://github.com/GridMod/RTS-GMLC (record the repository commit SHA via the GitHub API): RTS_Data/timeseries_data_files/WIND/{DAY_AHEAD_wind.csv, REAL_TIME_wind.csv}, PV/{DAY_AHEAD_pv.csv, REAL_TIME_pv.csv}, Load/{DAY_AHEAD_regional_Load.csv, REAL_TIME_regional_Load.csv}, RTS_Data/SourceData/{gen.csv, bus.csv}. Save under data/rts_gmlc/ with data/rts_gmlc/SOURCE.md listing URL, repo commit, file SHA256, units and time conventions (2020, day-ahead hourly 8784 periods, real-time 5-minute 105408 periods). Add one line to data/README.md describing the folder.
2. Forecast errors: real-time averaged to hourly minus day-ahead, for the 4 wind farms (309_WIND_1, 317_WIND_1, 303_WIND_1, 122_WIND_1, total 2508 MW) and the 25 utility PV plants (daylight hours only for PV). Report per site: std (MW and % of capacity), excess kurtosis, 1%/99% quantiles; pairwise correlation of errors and of squared standardized errors (common fluctuation strength indicator); dependence of error spread on forecast level. A quick check already found wind excess kurtosis 2.8–4.1, mean error correlation 0.40, mean squared-error correlation 0.31, and negligible load forecast error. Push.

## Stage B: Network, renewable placement, operating points, balancing

1. Network: data/networks/pglib118/pglib_opf_case118_ieee.m (MATPOWER format, has branch ratings rateA). Parse it in Python (write a small .m parser or use pandapower/matpowercaseframes) and build a pandapower or PYPOWER model; check bus/gen/branch counts and total load against the file.
2. Partition IEEE118 into three electrically connected regions by a documented rule (for example spectral clustering on the admittance graph, or a literature three-area division if you can cite it) and record the bus lists. Map RTS-GMLC areas 1/2/3 to these regions. Place the 4 wind farms and a subset of UPV plants (aggregate plants of the same RTS area into 2–4 IEEE118 host buses per region) at non-radial load or PQ buses with adequate line capacity; record the rule and the final mapping table.
3. Scale renewable capacity so that total installed renewable capacity is 25% of the IEEE118 total base active load (record capacities). Renewables are PQ injections at unity power factor (Q=0). If the base operating points are infeasible or lines are overloaded in the base case, try 20% once, record the reason; do not keep searching for favourable settings.
4. Operating points: choose two representative hours from 2020 using day-ahead forecasts: OP1 high aggregate wind (night, aggregate wind forecast between the 75th and 90th percentile), OP2 midday high PV with moderate wind. Loads: scale IEEE118 base loads per region by the RTS regional day-ahead load at that hour divided by the RTS regional annual mean (constant power factor). Generator dispatch: obtain a feasible base dispatch at each operating point (net load = load minus renewable forecast) with an AC OPF (pandapower runopp or PYPOWER runopf) used only as initialization; record it.
5. Active-power balancing (author decision): the imbalance caused by renewable and load deviations is shared by all online conventional units in proportion to participation factors proportional to their Pmax (record them); the reference bus only fixes the angle. PV units enforce reactive limits (a unit at its limit becomes PQ with Q fixed at the limit; all violating units switch together; no switching back), the reference unit has no reactive limit. Use pandapower runpp(distributed_slack=True, enforce_q_lims=True) or an equivalent, AND write your own Newton-Raphson in cloud/tmm_cloud (rectangular coordinates, unknowns = non-reference bus voltages plus the scalar system imbalance; equations = active power at all buses incl. reference, reactive power at PQ buses, squared voltage magnitude at PV buses; same reactive-limit rule). Cross-check the two solvers on at least 200 random scenarios per operating point (max voltage difference, unit P/Q difference, identical final modes). Report base-case voltages vs bus Vmin/Vmax and branch loading vs rateA. Push.

## Stage C: Forecast-error model with common fluctuation strength

Author decision: fitted model, not plain resampling. For each operating point: select historical hours whose day-ahead forecasts are close to the operating point (k-nearest neighbours on the capacity-normalized forecast vector of the mapped sites, k about 500–1000, record k and distance). Use their errors. Model: empirical marginal per site (kernel-smoothed, respecting physical bounds 0 ≤ forecast+error ≤ capacity), dependence by a t copula whose degrees of freedom is estimated by maximum likelihood (this carries the common fluctuation strength); also fit a Gaussian copula with the same marginals as the light-tailed counterpart. Check that the fitted model reproduces per-site excess kurtosis, pairwise correlation and squared-error correlation of the historical errors (table). Loads: small independent fluctuation (coefficient of variation 1%, Gaussian, constant power factor) on top of the regional profile, documented as a secondary random source. Push.

## Stage D: Scenarios and AC Monte Carlo

For each operating point draw independent sets: modeling 5256, determination 1752, check 100000 (heavy-tailed t-copula model), plus a light-tailed check set of 20000 (Gaussian copula). Solve every scenario with the AC power flow of Stage B (your verified Newton or pandapower, parallelized). Save per scenario: complex bus voltages, unit P and Q, final unit modes, branch P/Q at both ends, system imbalance, convergence flag and iterations; record failures with reasons. Record computation times. Push after each operating point.

## Stage E: Does the heavy tail reach transmission voltages and line flows?

For every non-reference bus voltage magnitude and angle and every branch P/Q end: mean, std, excess kurtosis, 0.5%/1%/99%/99.5% quantiles, (q99.5−q0.5)/IQR, for the heavy-tailed and light-tailed check sets; list the objects with the largest tail difference; relate them to the voltage sensitivity to renewable injections (columns of the inverse power-flow Jacobian for the renewable buses at the operating point). Also estimate, with 95% binomial confidence intervals, voltage-limit violation probabilities (bus Vmin/Vmax) and branch overload probabilities (apparent power > rateA) under both input models, and the frequency of conventional units leaving [Pmin, Pmax] because of the shared imbalance. Write the stage README with a clear judgment: does the heavy tail propagate (yes/no/partly, where), and are violation probabilities in the 1e-2–1e-3 range present (needed for later tail checks). Push.

## Stage F (only if time remains): Python batch power flow

Follow 01_IDEA/derivation.md Sections 1–2 (expansion point, second-order initial value eq. (8), fixed-Jacobian corrections eq. (9)–(10), reactive-limit mode grouping Section 2.3) extended with the distributed-slack unknown; verify against your Newton on the check scenarios (identical modes, voltage error < 1e-7 p.u.) and report counts of expansion points and timing. Push.

## Final report

Write cloud/results/phase1_rts118/SUMMARY.md and push (Chinese, power-system wording, no forbidden words): branch name and last pushed commit; per stage what was done and the key numbers (data statistics, mapping and capacity, operating points, solver cross-check errors, error-model checks, Monte Carlo counts/failures/times, tail-propagation findings with the main objects, violation probabilities with intervals); every decision you had to make and why; problems or deviations from this plan; what remains. Be honest about anything not done or not working.
