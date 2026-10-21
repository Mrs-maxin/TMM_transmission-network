function summary = run_tmm_transmission_accuracy(networkName, layer, options)
%RUN_TMM_TRANSMISSION_ACCURACY Accuracy of TMM AC probabilistic power flow on PV/PQ grids.
% Claim C0. Source-load scenarios -> fast PQ/PV power flow -> node-voltage
% scenarios -> joint voltage model (GMM, standard TMM, block TMM of AQPPF_TMM
% derivation 2026-10-21) -> distributions of node voltage magnitude and
% angle and of branch active and reactive power -> comparison with an
% independent AC Monte Carlo. LAYER 1 applies no reactive limits; LAYER 2
% limits PV units (slack unlimited, MATPOWER pf.enforce_q_lims=1 rule).
%
% Steps:
% 1. Configuration and external dependencies.
% 2. Source-load scenarios (modeling, component-count determination, check).
% 3. Fast power flow of modeling and determination scenarios.
% 4. AC Monte Carlo reference (reused when the same check scenarios exist).
% 5. Fast power flow of the check scenarios paired with the reference.
% 6. GMM and standard TMM over component counts.
% 7. Common-strength voltage directions and block TMM over component counts.
% 8. Distributions and accuracy against the AC reference.
% 9. PDF comparison figures, tables, summary and calculation record.
arguments
    networkName (1,1) string {mustBeMember(networkName, ["ieee14", "ieee30", "ieee57", "ieee118"])}
    layer (1,1) double {mustBeMember(layer, [1, 2])}
    options.Initialization (1,1) string {mustBeMember(options.Initialization, ["kmeans", "balanced-physical"])} = "kmeans"
    options.ComponentCounts double = [3, 6, 9, 12, 15, 20]
    options.ReferenceRunDir (1,1) string = ""
end
projectRoot = string(fileparts(fileparts(mfilename("fullpath"))));
addpath(char(projectRoot));
config = local_config(projectRoot, networkName, layer, options);
tag = sprintf("[%s layer %d accuracy]", networkName, layer);
step = @(k, text) fprintf("%s step %d/9 %s (%s)\n", tag, k, text, string(datetime("now", "Format", "HH:mm:ss")));

step(1, "configuration and dependencies");
dependencies = TMM_Transmission_utils.io.setup_external_dependencies(config);
TMM_Transmission_utils.io.ensure_dir(config.resultDir);
mpc = loadcase(char(config.networkFile));
limited = layer == 2;
if limited
    mpc = TMM_Transmission_core.powerflow.release_slack_reactive_limits(mpc);
end
control = TMM_Transmission_core.powerflow.describe_bus_control(mpc);
nr = control.nonReference;
mapOptions = {"MaximumAnchorCount", config.maximumAnchorCount, "StateTolerance", config.stateTolerance, ...
    "PowerTolerance", config.powerTolerance, "MaximumCorrections", config.maximumCorrections, ...
    "PowerFlowTolerance", config.powerFlowTolerance};

step(2, "source-load scenarios");
modeling = TMM_Transmission_utils.scenarios.build_correlated_load_samples(mpc, config.modelingCount, config.seed + 1);
determination = TMM_Transmission_utils.scenarios.build_correlated_load_samples(mpc, config.determinationCount, config.seed + 2);
check = TMM_Transmission_utils.scenarios.build_correlated_load_samples(mpc, config.checkCount, config.seed + 3);

step(3, "fast power flow of modeling and determination scenarios");
modelingMap = local_map(mpc, modeling, limited, mapOptions, config);
determinationMap = local_map(mpc, determination, limited, mapOptions, config);

step(4, "AC Monte Carlo reference");
if strlength(config.referenceReusedFrom) > 0
    saved = load(fullfile(config.referenceReusedFrom, "states.mat"), "test", "reference");
    assert(isequal(saved.test.loadP, check.loadP) && isequal(saved.test.loadQ, check.loadQ), ...
        "Check scenarios in %s differ from the regenerated scenarios.", config.referenceReusedFrom);
    reference = saved.reference;
    clear saved
else
    reference = TMM_Transmission_core.powerflow.solve_reference_power_flow(mpc, check.loadP, check.loadQ, ...
        PowerFlowTolerance=config.powerFlowTolerance, EnforceReactiveLimits=limited);
end
assert(all(reference.success), "%d Monte Carlo power flows failed.", sum(~reference.success));
if ~isfield(reference, "modeState"), reference.modeState = zeros(size(reference.voltage, 1), size(mpc.gen, 1)); end

step(5, "fast power flow of the check scenarios paired with the reference");
checkMap = local_map(mpc, check, limited, mapOptions, config);
paired = struct("count", size(check.loadP, 1), "sameModeCount", sum(all(checkMap.modeState == reference.modeState, 2)), ...
    "maximumComplexVoltageErrorPu", max(abs(checkMap.voltage - reference.voltage), [], "all"), ...
    "modeCount", size(unique(reference.modeState, "rows"), 1), ...
    "fastSeconds", checkMap.seconds, "referenceSeconds", reference.seconds);
fprintf("  paired: same mode %d/%d, max voltage error %.3g p.u.\n", paired.sameModeCount, paired.count, ...
    paired.maximumComplexVoltageErrorPu);

step(6, "GMM and standard TMM over component counts");
fitConfig = config.fit;
fitConfig.resultDir = fullfile(config.resultDir, "fits");
fits = TMM_Transmission_core.statistics.fit_voltage_mixture_orders(modelingMap.voltage, ...
    determinationMap.voltage, control, fitConfig);

step(7, "common-strength voltage directions and block TMM");
[modes, ~, group] = unique(modelingMap.modeState, "rows");
[~, dominant] = max(accumarray(group, 1));
directions = TMM_Transmission_core.statistics.determine_common_strength_directions(mpc, ...
    modeling.loadP, modeling.loadQ, modes(dominant, :), fits.transform, PowerFlowTolerance=config.powerFlowTolerance);
directionCheck = local_direction_check(modelingMap.voltage, nr, fits.transform, directions);
blockFits = TMM_Transmission_core.statistics.fit_block_tmm_orders(fits, directions, modelingMap.voltage, ...
    determinationMap.voltage, control, config.fit);
fprintf("  directions r=%d of %d; block TMM K=%d\n", directions.rank, fits.rank, blockFits.selectedK);

step(8, "distributions and accuracy against the AC reference");
slackProfile = modelingMap.voltage(1, :);
observed = TMM_core.network.evaluate_outputs(reference.voltage, mpc);
observed.voltageAngle = rad2deg(angle(reference.voltage(:, nr)));
floorPrediction = TMM_core.network.evaluate_outputs(modelingMap.voltage, mpc);
floorPrediction.voltageAngle = rad2deg(angle(modelingMap.voltage(:, nr)));
floorPrediction.sampleWeights = ones(size(modelingMap.voltage, 1), 1);
floor = TMM_Transmission_core.metrics.compare_power_distributions(reference.voltage, floorPrediction, mpc, control, "建模场景");
errorTables = local_tag(floor.errors, 0);
accuracy = TMM_Transmission_core.metrics.summarize_accuracy_by_quantity(floor.errors, "建模场景", 0);
selected = struct("GMM", fits.selected.gmmK, "TMM", fits.selected.tmmK, "blockTMM", blockFits.selectedK);
plotted = struct([]);
for method = ["GMM", "TMM", "blockTMM"]
    k = selected.(method);
    index = find(fits.componentCounts == k);
    if method == "blockTMM", model = blockFits.models{index}; else, model = fits.pairs{index}.(lower(method)); end
    prediction = TMM_Transmission_core.distributions.propagate_voltage_model(model, fits.transform, mpc, slackProfile, config.distribution);
    label = local_label(method);
    comparison = TMM_Transmission_core.metrics.compare_power_distributions(reference.voltage, prediction, mpc, control, label);
    errorTables = [errorTables; local_tag(comparison.errors, k)]; %#ok<AGROW>
    accuracy = [accuracy, TMM_Transmission_core.metrics.summarize_accuracy_by_quantity(comparison.errors, label, k)]; %#ok<AGROW>
    entry = struct("name", label, "style", local_style(method), "K", k, "prediction", prediction);
    if isempty(plotted), plotted = entry; else, plotted(end+1) = entry; end %#ok<AGROW>
    fprintf("  %s K=%d done\n", method, k);
end

step(9, "PDF comparison figures, tables, summary and calculation record");
objects = local_plot_objects(observed, errorTables, mpc, nr);
figures = TMM_Transmission_utils.plotting.export_pdf_comparison(observed, plotted, objects, fullfile(config.resultDir, "pdf"));
accuracyTable = struct2table(accuracy);
writetable(accuracyTable, fullfile(config.resultDir, "accuracy_by_quantity.csv"));
writetable(errorTables, fullfile(config.resultDir, "distribution_errors_all_objects.csv"));
writetable(fits.status, fullfile(config.resultDir, "fit_status_gmm_tmm.csv"));
writetable(blockFits.status, fullfile(config.resultDir, "fit_status_block_tmm.csv"));
save(fullfile(config.resultDir, "fitted_models.mat"), "fits", "blockFits", "slackProfile", "-v7.3");
summary = struct("study", sprintf("%s TMM AC probabilistic power flow accuracy, layer %d", networkName, layer), ...
    "claim", "C0", "network", networkName, "layer", layer, "pairedFastVersusMatpower", paired, ...
    "voltageCoordinateCount", fits.rank, "commonStrengthDirections", directionCheck, ...
    "selectedK", selected, "accuracy", accuracy, "figures", figures);
local_write_json(fullfile(config.resultDir, "summary.json"), summary);
record = struct("entry", "tests/run_tmm_transmission_accuracy.m", "arguments", [networkName, string(layer)], ...
    "config", config, "dependencies", dependencies, "networkSha256", local_sha256(config.networkFile), ...
    "projectGit", "【Git 待提交】 local project is not a Git repository", ...
    "finishedAt", string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), "componentCountUsesReference", false);
local_write_json(fullfile(config.resultDir, "run_manifest.json"), record);
fprintf("%s done: %s\n", tag, config.resultDir);
end

function config = local_config(projectRoot, networkName, layer, options)
files = struct("ieee14", "case14.m", "ieee30", "case_ieee30.m", "ieee57", "case57.m", "ieee118", "case118.m");
stamp = string(datetime("now", "Format", "yyyyMMdd'T'HHmmss"));
yearRoot = fileparts(char(projectRoot));
phdRoot = fileparts(yearRoot);
config = struct();
config.caseName = networkName + "_layer" + layer + "_accuracy";
config.networkFile = fullfile(projectRoot, "data", "networks", networkName, files.(networkName));
config.matpowerRoot = string(fullfile(phdRoot, "00codezip", "matpower8.0"));
config.aqppfRoot = string(fullfile(yearRoot, "[1]ACPPF_TMM"));
config.aqppfReferenceCommit = "0d5c999db4cbee2ea7055844d31b729d9694d747";
config.blockTmmSource = "derivation Section 3.3-3.4; AQPPF_TMM 01_IDEA/derivation.md, section 2026-10-21 (B1)-(B16)";
config.resultDir = fullfile(projectRoot, "result", config.caseName, stamp);
config.reactiveLimitRule = "none";
if layer == 2
    config.reactiveLimitRule = "PV units limited by source Qmin/Qmax; slack unlimited; MATPOWER pf.enforce_q_lims=1";
end
config.referenceReusedFrom = options.ReferenceRunDir;
config.seed = 20261021;
config.modelingCount = 5256;
config.determinationCount = 1752;
config.checkCount = 20000;
config.maximumAnchorCount = 12;
config.stateTolerance = 1e-7;
config.powerTolerance = 1e-8;
config.maximumCorrections = 10;
config.powerFlowTolerance = 1e-10;
config.reactiveTolerance = 5e-6;
config.fit = struct();
config.fit.componentCounts = options.ComponentCounts;
config.fit.pcaVarianceTolerance = 1e-10;
config.fit.blockMaximumIterations = 300;
config.fit.fitter = struct("seed", config.seed, "initialization", options.Initialization, "selectionMode", "likelihood", ...
    "restarts", 3, "maximumIterations", 500, "likelihoodTolerance", 1e-6, "covarianceFloor", 1e-7, ...
    "projectInitialScale", true, "initialDf", [12, 200], "dfUpdateMethod", "latent-em", ...
    "dfBounds", [4.1, 1e5], "emAcceleration", "none", "saveCandidateModels", false);
config.distribution = struct("modelDrawCount", 50000, "seed", config.seed + 100, "propagationBatchSize", 5000);
end

function mapping = local_map(mpc, scenarios, limited, mapOptions, config)
if limited
    mapping = TMM_Transmission_core.powerflow.map_load_samples_with_reactive_limits(mpc, scenarios.loadP, scenarios.loadQ, ...
        mapOptions{:}, "ReactiveTolerance", config.reactiveTolerance);
    return
end
solved = TMM_Transmission_core.powerflow.map_load_samples_to_voltage(mpc, scenarios.loadP, scenarios.loadQ, mapOptions{:});
mapping = struct("voltage", solved.voltage, "modeState", zeros(size(solved.voltage, 1), size(mpc.gen, 1)), ...
    "seconds", solved.totalSeconds);
end

function check = local_direction_check(voltage, nr, transform, directions)
x = (([real(voltage(:, nr)), imag(voltage(:, nr))] - transform.center)*transform.basis)./transform.scale;
U = x*directions.Q1; W = x*directions.Q2;
physical = x.*transform.scale;
uShare = sum((physical*directions.Q1).^2, "all")/sum(physical.^2, "all");
check = struct("rank", directions.rank, "drivenInjectionCount", directions.drivenInjectionCount, ...
    "voltageCoordinateCount", transform.rank, "uExcessKurtosisMedian", median(kurtosis(U) - 3), ...
    "wExcessKurtosisMedian", median(kurtosis(W) - 3), "uVarianceShareInFittingCoordinates", uShare);
end

function objects = local_plot_objects(observed, errorTables, mpc, nr)
pqColumns = find(mpc.bus(nr, 2) == 1);
[~, v] = max(std(observed.voltage(:, pqColumns), 1, 1)); v = pqColumns(v);
[~, a] = max(std(observed.voltageAngle, 1, 1));
[~, p] = max(std(observed.linePfrom, 1, 1));
[~, q] = max(std(observed.lineQfrom, 1, 1));
branch = mpc.branch(mpc.branch(:, 11) > 0, 1:2);
objects = [ ...
    struct("quantity", "voltage", "column", v, "file", "V" + nr(v), "label", "节点" + nr(v) + "电压幅值 (p.u.)"), ...
    struct("quantity", "voltageAngle", "column", a, "file", "A" + nr(a), "label", "节点" + nr(a) + "电压相角 (°)"), ...
    struct("quantity", "linePfrom", "column", p, "file", "P" + branch(p, 1) + "-" + branch(p, 2), ...
        "label", branch(p, 1) + "–" + branch(p, 2) + "支路" + branch(p, 1) + "端有功 (MW)"), ...
    struct("quantity", "lineQfrom", "column", q, "file", "Q" + branch(q, 1) + "-" + branch(q, 2), ...
        "label", branch(q, 1) + "–" + branch(q, 2) + "支路" + branch(q, 1) + "端无功 (MVAr)")];
tmmRows = errorTables(errorTables.method == "TMM" & ismember(errorTables.quantity, ["linePfrom", "lineQfrom"]) ...
    & errorTables.referenceVariance > 0.05^2, :);
[~, worst] = max(tmmRows.KS);
quantity = tmmRows.quantity(worst); column = tmmRows.index(worst);
if ~any(arrayfun(@(o) o.quantity == quantity && o.column == column, objects))
    unit = "MW"; kind = "有功"; if quantity == "lineQfrom", unit = "MVAr"; kind = "无功"; end
    objects(end+1) = struct("quantity", quantity, "column", column, ...
        "file", "Worst_" + extractBefore(quantity, 6) + branch(column, 1) + "-" + branch(column, 2), ...
        "label", branch(column, 1) + "–" + branch(column, 2) + "支路" + branch(column, 1) + "端" + kind + "，TMM误差最大 (" + unit + ")");
end
end

function label = local_label(method)
label = method;
if method == "blockTMM", label = "分块TMM"; end
end

function style = local_style(method)
style = "proposed";
if method == "GMM", style = "gmm"; elseif method == "blockTMM", style = "block_tmm"; end
end

function tagged = local_tag(errors, k)
tagged = errors;
tagged.K = repmat(k, height(errors), 1);
end

function local_write_json(file, value)
fid = fopen(file, "w");
assert(fid > 0, "Cannot write %s.", file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(value, PrettyPrint=true));
end

function digest = local_sha256(file)
[status, text] = system(sprintf('shasum -a 256 "%s"', file));
digest = "unavailable";
if status == 0, digest = extractBefore(string(strtrim(text)), " "); end
end
