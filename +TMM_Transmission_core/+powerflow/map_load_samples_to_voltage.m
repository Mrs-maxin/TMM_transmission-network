function mapping = map_load_samples_to_voltage(mpc, loadP, loadQ, options)
%MAP_LOAD_SAMPLES_TO_VOLTAGE Load samples to node voltages without AC Monte Carlo.
% Transmission counterpart of AQPPF map_samples_with_atlas plus
% recover_training_voltage: adaptive anchors (one AC solve each), batch
% quadratic inverse (AQPPF map_power_samples_quadratic, unchanged), fixed
% Jacobian corrections, and exact AC only for samples still above tolerance.
arguments
    mpc (1,1) struct
    loadP (:,:) double
    loadQ (:,:) double
    options.MaximumAnchorCount (1,1) double {mustBeInteger, mustBePositive} = 12
    options.StateTolerance (1,1) double {mustBePositive} = 1e-7
    options.PowerTolerance (1,1) double {mustBePositive} = 1e-8
    options.MaximumCorrections (1,1) double {mustBeInteger, mustBeNonnegative} = 10
    options.PowerFlowTolerance (1,1) double {mustBePositive} = 1e-10
    options.VoltageMin (1,1) double = 0.5
    options.VoltageMax (1,1) double = 1.5
end
timer = tic;
control = TMM_Transmission_core.powerflow.describe_bus_control(mpc);
target = TMM_Transmission_core.powerflow.build_mixed_power_targets(control, loadP, loadQ);
count = size(target, 2);
[~, nextIndex] = min(vecnorm(target - mean(target, 2), 2, 1));
bestState = nan(size(target));
bestScore = inf(1, count);
bestAnchor = zeros(1, count);
accepted = false(1, count);
corrections = zeros(1, count);
anchors = {};
anchorRows = zeros(1, 0);
for anchorNumber = 1:options.MaximumAnchorCount
    anchor = TMM_Transmission_core.powerflow.build_mixed_anchor(mpc, control, ...
        loadP(nextIndex, :).', loadQ(nextIndex, :).', PowerFlowTolerance=options.PowerFlowTolerance);
    anchors{end+1} = anchor; %#ok<AGROW>
    anchorRows(end+1) = nextIndex; %#ok<AGROW>
    active = find(~accepted);
    inverseAction = @(rhs) anchor.inverseMatrix*rhs;
    powerMap = @(x) TMM_Transmission_core.powerflow.evaluate_mixed_power(x, anchor);
    mapped = TMM_core.inverse.map_power_samples_quadratic(anchor.state, anchor.power, ...
        target(:, active), inverseAction, powerMap);
    mapped.state(:, active == nextIndex) = anchor.state;
    refined = TMM_Transmission_core.powerflow.refine_mixed_jacobian(mapped.state, target(:, active), anchor, ...
        MaximumCorrections=options.MaximumCorrections, StateTolerance=options.StateTolerance, ...
        PowerTolerance=options.PowerTolerance, VoltageMin=options.VoltageMin, VoltageMax=options.VoltageMax);
    score = refined.stateResidualEstimate;
    score(~refined.voltageInside | ~isfinite(score)) = inf;
    improve = score < bestScore(active) | refined.accepted;
    selected = active(improve);
    bestState(:, selected) = refined.state(:, improve);
    bestScore(selected) = score(improve);
    bestAnchor(selected) = anchorNumber;
    corrections(selected) = refined.correctionCount(improve);
    accepted(selected) = refined.accepted(improve);
    if all(accepted), break; end
    ranking = bestScore;
    ranking(accepted) = -inf;
    ranking(anchorRows) = -inf;
    [worst, nextIndex] = max(ranking);
    if ~(worst > 0), break; end
end

fallbackRows = find(~accepted);
reference = TMM_Transmission_core.powerflow.solve_reference_power_flow(mpc, ...
    loadP(fallbackRows, :), loadQ(fallbackRows, :), PowerFlowTolerance=options.PowerFlowTolerance);
nr = control.nonReference; m = numel(nr);
voltage = repmat(anchors{1}.voltage.', count, 1);
voltage(:, nr) = complex(bestState(1:m, :), bestState(m+1:end, :)).';
voltage(fallbackRows, :) = reference.voltage;
assert(all(reference.success), "TMM_Transmission:powerflow:FallbackFailed", ...
    "Exact AC fallback failed for %d load samples.", sum(~reference.success));

finalResidual = target - TMM_Transmission_core.powerflow.evaluate_mixed_power( ...
    [real(voltage(:, nr)).'; imag(voltage(:, nr)).'], anchors{1});
maximumResidual = max(abs(finalResidual), [], "all");
assert(maximumResidual <= max(options.PowerTolerance, 10*options.PowerFlowTolerance), ...
    "TMM_Transmission:powerflow:FinalResidual", ...
    "Final PQ/PV residual %.3g p.u. exceeds tolerance.", maximumResidual);
mapping = struct("voltage", voltage, "control", control, "anchorCount", numel(anchors), ...
    "anchorRows", anchorRows, "assignedAnchor", bestAnchor, "acceptedByFastMapping", accepted, ...
    "fallbackRows", fallbackRows, "fallbackCount", numel(fallbackRows), ...
    "exactAcCallCount", numel(anchors) + numel(fallbackRows), ...
    "correctionCount", corrections, "maximumFinalResidualPu", maximumResidual, ...
    "fastMappingSeconds", toc(timer) - reference.seconds, "fallbackSeconds", reference.seconds, ...
    "totalSeconds", toc(timer));
mapping.anchors = anchors;
end
