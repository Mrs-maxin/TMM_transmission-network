function comparison = compare_power_distributions(referenceVoltage, prediction, mpc, control, method)
%COMPARE_POWER_DISTRIBUTIONS Model distributions against independent AC Monte Carlo.
% Distribution errors use AQPPF TMM_core.metrics.compare_distributions
% unchanged (KS, W1, mean, variance, 1%/99% quantiles). PV terminal voltages
% are held by control in the AC reference, so their model spread is reported
% separately as the share of draws away from the set point.
arguments
    referenceVoltage (:,:) double
    prediction (1,1) struct
    mpc (1,1) struct
    control (1,1) struct
    method (1,1) string
end
reference = TMM_core.network.evaluate_outputs(referenceVoltage, mpc);
fields = ["voltage", "nodeP", "nodeQ", "linePfrom", "lineQfrom", "linePto", "lineQto", "loss"];
if isfield(prediction, "voltageAngle")
    reference.voltageAngle = rad2deg(angle(referenceVoltage(:, control.nonReference)));
    fields = ["voltage", "voltageAngle", fields(2:end)];
end
model = struct("sampleWeights", prediction.sampleWeights);
observed = struct();
for name = fields
    model.(name) = prediction.(name);
    observed.(name) = reference.(name);
end
errors = TMM_core.metrics.compare_distributions(observed, model, method);
nr = control.nonReference;
busNumber = nr(:);
isVoltage = ismember(errors.quantity, ["voltage", "voltageAngle"]);
errors.bus = nan(height(errors), 1);
errors.bus(isVoltage) = busNumber(errors.index(isVoltage));
errors.busType = repmat("", height(errors), 1);
errors.busType(isVoltage & ismember(errors.bus, control.pv)) = "PV";
errors.busType(isVoltage & ismember(errors.bus, control.pq)) = "PQ";

[~, pvColumns] = ismember(control.pv, nr);
deviation = abs(prediction.voltage(:, pvColumns) - control.voltageSet(control.pv).');
referenceDeviation = abs(reference.voltage(:, pvColumns) - control.voltageSet(control.pv).');
weights = prediction.sampleWeights(:)/sum(prediction.sampleWeights);
pvLeakage = table2struct(array2table([control.pv(:), (weights.'*(deviation > 1e-4)).', ...
    (weights.'*(deviation > 1e-3)).', max(deviation, [], 1).', ...
    mean(referenceDeviation > 1e-4, 1).', max(referenceDeviation, [], 1).'], ...
    "VariableNames", ["bus", "shareAbove1e4pu", "shareAbove1e3pu", "maximumDeviationPu", ...
    "referenceShareAbove1e4pu", "referenceMaximumDeviationPu"]));
comparison = struct("errors", errors, "pvVoltageSpread", pvLeakage, "method", method);
end
