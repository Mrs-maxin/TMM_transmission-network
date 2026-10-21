function [model, audit] = fit_block_tmm(U, W, gmm, initialDf, options)
%FIT_BLOCK_TMM Expectation conditional maximization of the block TMM.
% Derivation eqs. (19)-(21), after AQPPF_TMM derivation 2026-10-21 (B12)-(B15): the U block updates
% location, scatter and df as a multivariate t; the W block updates intercept,
% regression on U and residual covariance by weighted least squares. The start
% is the same-K GMM split into the two blocks, scatter scaled by (nu-2)/nu.
% Floors, df bounds, iteration limit and stopping rule follow the standard TMM.
arguments
    U (:,:) double
    W (:,:) double
    gmm (1,1) struct
    initialDf (1,1) double {mustBePositive}
    options.MaximumIterations (1,1) double {mustBeInteger, mustBePositive} = 300
    options.LikelihoodTolerance (1,1) double {mustBePositive} = 1e-6
    options.CovarianceFloor (1,1) double {mustBePositive} = 1e-7
    options.DfBounds (1,2) double = [4.1, 1e5]
    options.Q1 double = []
    options.Q2 double = []
end
timer = tic;
[n, r] = size(U);
dw = size(W, 2);
k = numel(gmm.weights);
Q1 = options.Q1; Q2 = options.Q2;
model = struct("weights", gmm.weights(:).', "means", gmm.means*Q1, "scatter", zeros(r, r, k), ...
    "df", repmat(initialDf, 1, k), "intercept", zeros(k, dw), "coefficient", zeros(dw, r, k), ...
    "residual", zeros(dw, dw, k));
for j = 1:k
    C = gmm.scatter(:, :, j);
    Suu = Q1.'*C*Q1; Swu = Q2.'*C*Q1; Sww = Q2.'*C*Q2;
    G = Swu/Suu;
    model.scatter(:, :, j) = local_floor(Suu*(initialDf - 2)/initialDf, options.CovarianceFloor*(initialDf - 2)/initialDf);
    model.coefficient(:, :, j) = G;
    model.intercept(j, :) = gmm.means(j, :)*Q2 - (G*(gmm.means(j, :)*Q1).').';
    model.residual(:, :, j) = local_floor(Sww - G*Suu*G.', options.CovarianceFloor);
end
[logDensity, membership, distance] = TMM_Transmission_core.statistics.block_tmm_log_density(U, W, model);
history = mean(logDensity);
status = "maximum-iterations";
design = [ones(n, 1), U];
iteration = 0;
for iteration = 1:options.MaximumIterations
    counts = sum(membership, 1);
    if min(counts) <= dw + r + 2
        status = "collapsed-component";
        break
    end
    next = model;
    next.weights = counts/n;
    for j = 1:k
        nu = model.df(j);
        tau = (nu + r)./(nu + distance(:, j));
        weight = membership(:, j).*tau;
        next.means(j, :) = (weight.'*U)/sum(weight);
        deviation = U - next.means(j, :);
        next.scatter(:, :, j) = local_floor((deviation.'*(deviation.*weight))/counts(j), ...
            options.CovarianceFloor*(nu - 2)/nu);
        statistic = membership(:, j).'*(psi((nu + r)/2) - log((nu + distance(:, j))/2) - tau)/counts(j);
        score = @(v) log(v/2) - psi(v/2) + 1 + statistic;
        if score(options.DfBounds(1)) <= 0
            next.df(j) = options.DfBounds(1);
        elseif score(options.DfBounds(2)) >= 0
            next.df(j) = options.DfBounds(2);
        else
            next.df(j) = fzero(score, options.DfBounds);
        end
        weighted = design.*membership(:, j);
        b = (design.'*weighted)\(weighted.'*W);
        next.intercept(j, :) = b(1, :);
        next.coefficient(:, :, j) = b(2:end, :).';
        e = W - design*b;
        next.residual(:, :, j) = local_floor((e.'*(e.*membership(:, j)))/counts(j), options.CovarianceFloor);
    end
    [nextDensity, nextMembership, nextDistance] = TMM_Transmission_core.statistics.block_tmm_log_density(U, W, next);
    value = mean(nextDensity);
    gain = value - history(end);
    if gain < -1e-8
        status = "rejected-decreasing-update";
        break
    end
    model = next; membership = nextMembership; distance = nextDistance;
    history(end+1) = value; %#ok<AGROW>
    if abs(gain) < options.LikelihoodTolerance*max(1, abs(history(end-1)))
        status = "converged";
        break
    end
end
model.Q1 = Q1; model.Q2 = Q2; model.family = "block-tmm";
audit = struct("status", status, "iterations", iteration, "meanLogDensity", history(end), ...
    "df", model.df, "weights", model.weights, "effectiveCounts", sum(membership, 1), ...
    "initialDf", initialDf, "seconds", toc(timer));
end

function A = local_floor(A, floor)
A = (A + A.')/2;
[V, D] = eig(A);
A = V*diag(max(diag(D), floor))*V.';
A = (A + A.')/2;
end
