function [logDensity, membership, distance] = block_tmm_log_density(U, W, model)
%BLOCK_TMM_LOG_DENSITY Block TMM density, derivation eqs. (17)-(18).
% Within component k, U is multivariate t with location m_k, scatter Sigma_k
% and df nu_k; given U, W is Gaussian with mean a_k + G_k U and covariance
% Omega_k. Returns the mixture log density per scenario, the component
% membership probabilities and the t distances delta_ik.
arguments
    U (:,:) double
    W (:,:) double
    model (1,1) struct
end
[n, r] = size(U);
dw = size(W, 2);
k = numel(model.weights);
logs = zeros(n, k);
distance = zeros(n, k);
for j = 1:k
    L = chol(model.scatter(:, :, j), "lower");
    z = (L\(U - model.means(j, :)).').';
    delta = sum(z.^2, 2);
    nu = model.df(j);
    if isinf(nu)
        tLog = -0.5*(r*log(2*pi) + delta);
    else
        tLog = gammaln((nu + r)/2) - gammaln(nu/2) - r/2*log(nu*pi) - (nu + r)/2*log1p(delta/nu);
    end
    residualFactor = chol(model.residual(:, :, j), "lower");
    e = (residualFactor\(W - model.intercept(j, :) - U*model.coefficient(:, :, j).').').';
    gaussianLog = -0.5*(dw*log(2*pi) + sum(e.^2, 2)) - sum(log(diag(residualFactor)));
    logs(:, j) = log(model.weights(j)) + tLog - sum(log(diag(L))) + gaussianLog;
    distance(:, j) = delta;
end
peak = max(logs, [], 2);
logDensity = peak + log(sum(exp(logs - peak), 2));
membership = exp(logs - logDensity);
end
