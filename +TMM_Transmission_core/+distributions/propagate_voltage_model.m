function prediction = propagate_voltage_model(model, transform, mpc, slackProfile, config)
%PROPAGATE_VOLTAGE_MODEL Node-voltage magnitude/angle and power distributions of a voltage model.
% Component-stratified Sobol integration: GMM and standard TMM use AQPPF
% sample_stratified_coordinates unchanged; the block TMM draws its t block
% with one strength per point and the Gaussian block given it (derivation
% eq. (26)). Voltages are rebuilt on all non-slack buses and node and line
% powers follow from the exact AC relations (AQPPF evaluate_outputs).
arguments
    model (1,1) struct
    transform (1,1) struct
    mpc (1,1) struct
    slackProfile (1,:) double
    config (1,1) struct
end
timer = tic;
if isfield(model, "family") && model.family == "block-tmm"
    [x, weights] = local_sample_block(model, config.modelDrawCount, config.seed);
else
    [x, weights] = TMM_core.mixture.sample_stratified_coordinates(model, transform.scale, ...
        config.modelDrawCount, config.seed, Scheme="sobol");
end
pq = find(mpc.bus(:, 2) ~= 3);
m = numel(pq);
count = size(x, 1);
prediction = struct();
for first = 1:config.propagationBatchSize:count
    rows = first:min(count, first + config.propagationBatchSize - 1);
    rectangular = (x(rows, :).*transform.scale)*transform.basis.' + transform.center;
    voltage = repmat(slackProfile, numel(rows), 1);
    voltage(:, pq) = rectangular(:, 1:m) + 1i*rectangular(:, m+1:end);
    block = TMM_core.network.evaluate_outputs(voltage, mpc);
    block.voltageAngle = rad2deg(angle(voltage(:, pq)));
    for name = string(fieldnames(block)).'
        if first == 1, prediction.(name) = zeros(count, size(block.(name), 2)); end
        prediction.(name)(rows, :) = block.(name);
    end
end
prediction.sampleWeights = weights(:);
prediction.seconds = toc(timer);
end

function [x, weights] = local_sample_block(model, count, seed)
k = numel(model.weights);
r = size(model.means, 2);
dw = size(model.intercept, 2);
quota = count*model.weights;
allocation = max(1, floor(quota));
while sum(allocation) < count
    [~, j] = max(quota - allocation); allocation(j) = allocation(j) + 1;
end
while sum(allocation) > count
    candidates = find(allocation > 1);
    [~, j] = min(quota(candidates) - allocation(candidates)); allocation(candidates(j)) = allocation(candidates(j)) - 1;
end
x = zeros(count, size(model.Q1, 1));
weights = zeros(count, 1);
position = 0;
for j = 1:k
    nj = allocation(j);
    rng(seed + 1009*j, "twister");
    net = scramble(sobolset(r + dw + 1), "MatousekAffineOwen");
    u = net(1:nj, :);
    u = min(max(u, eps), 1 - eps);
    nu = model.df(j);
    if isinf(nu)
        strength = ones(nj, 1);
    else
        strength = 2*gammaincinv(u(:, end), nu/2)/nu;
    end
    U = model.means(j, :) + (norminv(u(:, 1:r))*chol(model.scatter(:, :, j), "lower").')./sqrt(strength);
    W = model.intercept(j, :) + U*model.coefficient(:, :, j).' + norminv(u(:, r+1:r+dw))*chol(model.residual(:, :, j), "lower").';
    rows = position + (1:nj);
    x(rows, :) = U*model.Q1.' + W*model.Q2.';
    weights(rows) = model.weights(j)/nj;
    position = position + nj;
end
end
