function directions = determine_common_strength_directions(mpc, loadP, loadQ, modeState, transform, options)
%DETERMINE_COMMON_STRENGTH_DIRECTIONS Voltage directions driven by the common fluctuation strength.
% Derivation eqs. (15)-(16), after AQPPF_TMM derivation 2026-10-21 (B7)-(B9): at an
% operating point, the inverse PQ/PV power-flow Jacobian columns of the
% injections that share one fluctuation strength span the voltage changes
% they cause. Mapped to the fitting coordinates and orthonormalized, they give
% Q1 (r directions, r <= number of such injections); Q2 is the complement.
% Here every load shares one strength, so the driven injections are the
% active load at all non-slack buses and the reactive load at PQ buses of the
% operating mode; reactive load at a regulating PV bus is supplied by its unit.
% The operating point is the mean modeling-scenario load in the given mode.
arguments
    mpc (1,1) struct
    loadP (:,:) double
    loadQ (:,:) double
    modeState (1,:) double
    transform (1,1) struct
    options.RankTolerance (1,1) double {mustBePositive} = 1e-10
    options.PowerFlowTolerance (1,1) double {mustBePositive} = 1e-10
end
define_constants;
modeMpc = TMM_Transmission_core.powerflow.apply_reactive_limit_mode(mpc, modeState);
control = TMM_Transmission_core.powerflow.describe_bus_control(modeMpc);
anchor = TMM_Transmission_core.powerflow.build_mixed_anchor(modeMpc, control, mean(loadP, 1).', mean(loadQ, 1).', ...
    PowerFlowTolerance=options.PowerFlowTolerance);
nr = control.nonReference;
activeRows = find(modeMpc.bus(nr, PD) ~= 0);
reactiveRows = numel(nr) + find(modeMpc.bus(control.pq, QD) ~= 0);
drivenRows = [activeRows; reactiveRows];
sensitivity = anchor.inverseMatrix(:, drivenRows);
assert(size(transform.basis, 1) == size(sensitivity, 1), "TMM_Transmission:statistics:CoordinateSize", ...
    "Voltage coordinates of the fit and of the power-flow Jacobian differ.");
whitened = (transform.basis.'*sensitivity)./transform.scale(:);
[left, singular, ~] = svd(whitened, "econ");
singular = diag(singular);
r = sum(singular > singular(1)*options.RankTolerance);
Q1 = left(:, 1:r);
[complement, ~, ~] = svd(eye(transform.rank) - Q1*Q1.');
Q2 = complement(:, 1:transform.rank - r);
assert(norm(Q1.'*Q2, "fro") < 1e-8, "TMM_Transmission:statistics:Orthogonality", ...
    "Common-strength and remaining directions are not orthogonal.");
directions = struct("Q1", Q1, "Q2", Q2, "rank", r, "drivenInjectionCount", numel(drivenRows), ...
    "singularValues", singular.', "modeState", modeState, ...
    "rule", "inverse PQ/PV Jacobian columns of load injections at the mean modeling-scenario load");
end
