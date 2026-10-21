function blockFits = fit_block_tmm_orders(fits, directions, trainingVoltage, validationVoltage, control, config)
%FIT_BLOCK_TMM_ORDERS Block TMM for every component count with an accepted GMM.
% Each K starts from that K's GMM with df 12 and 200; the converged candidate
% with the highest determination-scenario log likelihood is kept, and the
% component count is then chosen by the same likelihood (no Monte Carlo input).
arguments
    fits (1,1) struct
    directions (1,1) struct
    trainingVoltage (:,:) double
    validationVoltage (:,:) double
    control (1,1) struct
    config (1,1) struct
end
nr = control.nonReference;
transform = fits.transform;
x = (([real(trainingVoltage(:, nr)), imag(trainingVoltage(:, nr))] - transform.center)*transform.basis)./transform.scale;
z = (([real(validationVoltage(:, nr)), imag(validationVoltage(:, nr))] - transform.center)*transform.basis)./transform.scale;
U = x*directions.Q1; W = x*directions.Q2;
Uv = z*directions.Q1; Wv = z*directions.Q2;
models = cell(numel(fits.componentCounts), 1);
rows = struct([]);
for index = 1:numel(fits.componentCounts)
    k = fits.componentCounts(index);
    pair = fits.pairs{index};
    row = struct("K", k, "method", "block-TMM", "accepted", false, "validationLogDensity", -Inf, ...
        "dfMin", NaN, "dfMax", NaN, "status", "no accepted GMM start", "fitSeconds", 0);
    if ~isempty(pair.gmm) && pair.gmmAccepted
        best = [];
        for nu = config.fitter.initialDf(:).'
            [model, audit] = TMM_Transmission_core.statistics.fit_block_tmm(U, W, pair.gmm, nu, ...
                MaximumIterations=config.blockMaximumIterations, LikelihoodTolerance=config.fitter.likelihoodTolerance, ...
                CovarianceFloor=config.fitter.covarianceFloor, DfBounds=config.fitter.dfBounds(:).', ...
                Q1=directions.Q1, Q2=directions.Q2);
            row.fitSeconds = row.fitSeconds + audit.seconds;
            if audit.status == "collapsed-component", continue; end
            value = mean(TMM_Transmission_core.statistics.block_tmm_log_density(Uv, Wv, model));
            better = isempty(best) || (audit.status == "converged" && best.status ~= "converged") || ...
                ((audit.status == "converged") == (best.status == "converged") && value > best.value);
            if better
                best = struct("model", model, "value", value, "status", audit.status);
            end
        end
        if ~isempty(best)
            models{index} = best.model;
            row.accepted = best.status == "converged";
            row.validationLogDensity = best.value;
            row.dfMin = min(best.model.df); row.dfMax = max(best.model.df);
            row.status = best.status;
        end
    end
    if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
end
status = struct2table(rows);
eligible = find(status.accepted);
assert(~isempty(eligible), "TMM_Transmission:statistics:NoBlockModel", "No converged block TMM among the component counts.");
[~, best] = max(status.validationLogDensity(eligible));
blockFits = struct("models", {models}, "status", status, "selectedK", status.K(eligible(best)), ...
    "directions", directions);
end
