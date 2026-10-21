function fits = fit_voltage_mixture_orders(trainingVoltage, validationVoltage, control, config)
%FIT_VOLTAGE_MIXTURE_ORDERS Joint node-voltage GMM/TMM fits over component counts.
% Coordinates are the real and imaginary voltages of all non-slack buses, PV
% buses included. Each K is fitted by AQPPF TMM_core.mixture.fit_order_candidates
% unchanged (k-means restarts, GMM warm start, df starts, acceptance rule,
% validation-likelihood choice). K is then chosen by validation likelihood
% only; no AC Monte Carlo result enters the fit.
arguments
    trainingVoltage (:,:) double
    validationVoltage (:,:) double
    control (1,1) struct
    config (1,1) struct
end
nr = control.nonReference;
training = [real(trainingVoltage(:, nr)), imag(trainingVoltage(:, nr))];
validation = [real(validationVoltage(:, nr)), imag(validationVoltage(:, nr))];
[x, transform] = TMM_core.mixture.fit_coordinates(training, config.pcaVarianceTolerance);
z = ((validation - transform.center)*transform.basis)./transform.scale;
rows = struct([]);
pairs = cell(numel(config.componentCounts), 1);
for index = 1:numel(config.componentCounts)
    k = config.componentCounts(index);
    cfg = config.fitter;
    cfg.componentCount = k;
    cfg.resultDir = fullfile(config.resultDir, "K" + k);
    TMM_Transmission_utils.io.ensure_dir(cfg.resultDir);
    timer = tic;
    pair = TMM_core.mixture.fit_order_candidates(x, z, cfg, transform.scale);
    seconds = toc(timer);
    pair.transform = transform;
    pairs{index} = pair;
    for method = ["GMM", "TMM"]
        key = lower(method);
        row = struct("K", k, "method", method, "accepted", pair.(key + "Accepted"), ...
            "validationLogDensity", -Inf, "dfMin", NaN, "dfMax", NaN, "fitSecondsAllCandidates", seconds);
        if ~isempty(pair.(key))
            row.validationLogDensity = pair.audit.validationLogDensity(pair.(key + "SelectedRow"));
            row.dfMin = min(pair.(key).df);
            row.dfMax = max(pair.(key).df);
        end
        if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
    end
end
status = struct2table(rows);
selected = struct();
for method = ["GMM", "TMM"]
    eligible = find(status.method == method & status.accepted);
    assert(~isempty(eligible), "TMM_Transmission:statistics:NoAcceptedModel", ...
        "No accepted %s fit among the component counts.", method);
    [~, best] = max(status.validationLogDensity(eligible));
    selected.(lower(method) + "K") = status.K(eligible(best));
end
fits = struct("pairs", {pairs}, "componentCounts", config.componentCounts, "status", status, ...
    "transform", transform, "selected", selected, "trainingCount", size(x, 1), ...
    "validationCount", size(z, 1), "rank", transform.rank);
end
