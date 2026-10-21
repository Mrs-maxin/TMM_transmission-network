function reference = solve_reference_power_flow(mpc, loadP, loadQ, options)
%SOLVE_REFERENCE_POWER_FLOW Sample-by-sample MATPOWER Newton power flow.
% Conventional AC power flow with the source slack bus. With
% EnforceReactiveLimits, MATPOWER pf.enforce_q_lims=1 fixes every violating
% unit at its limit and resolves; pass data from release_slack_reactive_limits
% so that only PV units are limited (derivation Section 2.3). Used for exact fallback
% and as the independent Monte Carlo reference.
arguments
    mpc (1,1) struct
    loadP (:,:) double
    loadQ (:,:) double
    options.PowerFlowTolerance (1,1) double {mustBePositive} = 1e-10
    options.EnforceReactiveLimits (1,1) logical = false
end
define_constants;
count = size(loadP, 1);
n = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
mpopt = mpoption("verbose", 0, "out.all", 0, "pf.tol", options.PowerFlowTolerance, ...
    "pf.enforce_q_lims", double(options.EnforceReactiveLimits));
originalType = mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE).';
voltage = complex(nan(count, n));
generatorP = nan(count, ng);
generatorQ = nan(count, ng);
modeState = zeros(count, ng);
success = false(count, 1);
iterations = zeros(count, 1);
timer = tic;
for j = 1:count
    sample = mpc;
    sample.bus(:, PD) = loadP(j, :).';
    sample.bus(:, QD) = loadQ(j, :).';
    [solved, ok] = runpf(sample, mpopt);
    success(j) = ok;
    if isfield(solved, "iterations"), iterations(j) = solved.iterations; end
    if ok
        voltage(j, :) = (solved.bus(:, VM).*exp(1i*deg2rad(solved.bus(:, VA)))).';
        generatorP(j, :) = solved.gen(:, PG).';
        generatorQ(j, :) = solved.gen(:, QG).';
        converted = originalType == PV & solved.bus(solved.gen(:, GEN_BUS), BUS_TYPE).' == PQ;
        atUpper = abs(solved.gen(:, QG) - solved.gen(:, QMAX)).' < 1e-9;
        modeState(j, converted & atUpper) = 1;
        modeState(j, converted & ~atUpper) = -1;
    end
end
reference = struct("voltage", voltage, "generatorP", generatorP, "generatorQ", generatorQ, ...
    "modeState", modeState, "success", success, "iterations", iterations, "seconds", toc(timer), ...
    "solver", "MATPOWER runpf Newton, pf.enforce_q_lims=" + double(options.EnforceReactiveLimits));
end
