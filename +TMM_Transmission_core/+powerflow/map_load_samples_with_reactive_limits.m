function mapping = map_load_samples_with_reactive_limits(mpc, loadP, loadQ, options)
%MAP_LOAD_SAMPLES_WITH_REACTIVE_LIMITS Fast PQ/PV power flow with PV reactive limits.
% Every sample starts with all PV units regulating. After each batch solve,
% units outside their limits are fixed at the limit and their buses become
% PQ; samples are regrouped by operating mode and solved again until no new
% violation appears, without returning fixed units to regulation. This is
% the MATPOWER runpf rule (pf.enforce_q_lims=1). Each mode is solved by
% map_load_samples_to_voltage with its own anchors; the slack unit is
% unlimited (derivation Section 2.3).
arguments
    mpc (1,1) struct
    loadP (:,:) double
    loadQ (:,:) double
    options.MaximumAnchorCount (1,1) double {mustBeInteger, mustBePositive} = 12
    options.StateTolerance (1,1) double {mustBePositive} = 1e-7
    options.PowerTolerance (1,1) double {mustBePositive} = 1e-8
    options.MaximumCorrections (1,1) double {mustBeInteger, mustBeNonnegative} = 10
    options.PowerFlowTolerance (1,1) double {mustBePositive} = 1e-10
    options.ReactiveTolerance (1,1) double {mustBePositive} = 5e-6
end
timer = tic;
define_constants;
count = size(loadP, 1);
ng = size(mpc.gen, 1);
state = zeros(count, ng);
voltage = complex(nan(count, size(mpc.bus, 1)));
generatorQ = nan(count, ng);
rounds = zeros(count, 1);
active = (1:count).';
modeLog = struct([]);
mapOptions = {"MaximumAnchorCount", options.MaximumAnchorCount, "StateTolerance", options.StateTolerance, ...
    "PowerTolerance", options.PowerTolerance, "MaximumCorrections", options.MaximumCorrections, ...
    "PowerFlowTolerance", options.PowerFlowTolerance};
for pass = 1:ng + 1
    if isempty(active), break; end
    [modes, ~, group] = unique(state(active, :), "rows");
    stillActive = zeros(0, 1);
    for j = 1:size(modes, 1)
        rows = active(group == j);
        modeMpc = TMM_Transmission_core.powerflow.apply_reactive_limit_mode(mpc, modes(j, :));
        solved = TMM_Transmission_core.powerflow.map_load_samples_to_voltage(modeMpc, ...
            loadP(rows, :), loadQ(rows, :), mapOptions{:});
        voltage(rows, :) = solved.voltage;
        rounds(rows) = pass;
        [newState, q, violated] = TMM_Transmission_core.powerflow.find_reactive_limit_violations( ...
            modeMpc, modes(j, :), solved.voltage, loadQ(rows, :), options.ReactiveTolerance);
        generatorQ(rows, :) = q;
        state(rows, :) = newState;
        stillActive = [stillActive; rows(violated)]; %#ok<AGROW>
        entry = struct("pass", pass, "mode", mode_label(modes(j, :), mpc), "sampleCount", numel(rows), ...
            "anchorCount", solved.anchorCount, "fallbackCount", solved.fallbackCount, ...
            "exactAcCallCount", solved.exactAcCallCount, "newViolations", nnz(violated), ...
            "maximumFinalResidualPu", solved.maximumFinalResidualPu);
        if isempty(modeLog), modeLog = entry; else, modeLog(end+1) = entry; end %#ok<AGROW>
    end
    active = stillActive;
end
assert(isempty(active), "TMM_Transmission:powerflow:ModeIterationLimit", ...
    "%d samples still violate reactive limits after %d mode updates.", numel(active), ng + 1);
modeTable = struct2table(modeLog);
[finalModes, ~, finalGroup] = unique(state, "rows");
finalCounts = accumarray(finalGroup, 1);
finalLabels = strings(size(finalModes, 1), 1);
for j = 1:size(finalModes, 1), finalLabels(j) = mode_label(finalModes(j, :), mpc); end
mapping = struct("voltage", voltage, "modeState", state, "generatorQ", generatorQ, "rounds", rounds, ...
    "modeLog", modeTable, "finalModeLabels", finalLabels, "finalModeCounts", finalCounts, ...
    "exactAcCallCount", sum(modeTable.exactAcCallCount), "anchorCount", sum(modeTable.anchorCount), ...
    "fallbackCount", sum(modeTable.fallbackCount), "seconds", toc(timer));
end

function label = mode_label(state, mpc)
fixed = find(state ~= 0);
if isempty(fixed)
    label = "all PV regulating";
    return
end
parts = strings(1, numel(fixed));
for k = 1:numel(fixed)
    limitName = "Qmax";
    if state(fixed(k)) < 0, limitName = "Qmin"; end
    parts(k) = "G" + mpc.gen(fixed(k), 1) + "@" + limitName;
end
label = strjoin(parts, ", ");
end
