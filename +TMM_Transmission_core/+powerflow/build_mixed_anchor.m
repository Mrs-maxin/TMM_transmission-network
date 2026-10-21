function anchor = build_mixed_anchor(mpc, control, loadP, loadQ, options)
%BUILD_MIXED_ANCHOR One AC power flow and the inverse PQ/PV Jacobian.
% Transmission counterpart of AQPPF TMM_core.powerflow.build_anchor. The
% Jacobian is taken directly in rectangular voltages: P and Q rows from
% dSbus_dV, PV rows d|V|^2 = 2e de + 2f df (derivation eq. (7)).
arguments
    mpc (1,1) struct
    control (1,1) struct
    loadP (:,1) double
    loadQ (:,1) double
    options.PowerFlowTolerance (1,1) double {mustBePositive} = 1e-10
end
define_constants;
mpc.bus(:, PD) = loadP;
mpc.bus(:, QD) = loadQ;
mpopt = mpoption("verbose", 0, "out.all", 0, "pf.tol", options.PowerFlowTolerance, "pf.enforce_q_lims", 0);
[solved, success] = runpf(mpc, mpopt);
assert(success, "TMM_Transmission:powerflow:AnchorFailed", ...
    "AC power flow failed while building a transmission anchor.");
voltage = solved.bus(:, VM).*exp(1i*deg2rad(solved.bus(:, VA)));
Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
[dSdVr, dSdVi] = dSbus_dV(Ybus, voltage, 1);
nr = control.nonReference; m = numel(nr); npv = numel(control.pv);
voltageRows = zeros(npv, 2*m);
for k = 1:npv
    position = control.pvPosition(k);
    voltageRows(k, position) = 2*real(voltage(control.pv(k)));
    voltageRows(k, m + position) = 2*imag(voltage(control.pv(k)));
end
jacobian = full([real(dSdVr(nr, nr)), real(dSdVi(nr, nr)); ...
    imag(dSdVr(control.pq, nr)), imag(dSdVi(control.pq, nr)); voltageRows]);
inverseMatrix = jacobian\eye(2*m);
anchor = struct("control", control, "voltage", voltage, "Ybus", Ybus, ...
    "state", [real(voltage(nr)); imag(voltage(nr))], "jacobian", jacobian, ...
    "inverseMatrix", inverseMatrix, "inverseNorm", norm(inverseMatrix, 2), ...
    "admittanceNorm", norm(full(Ybus), 2), "reciprocalCondition", rcond(jacobian));
anchor.power = TMM_Transmission_core.powerflow.evaluate_mixed_power(anchor.state, anchor);
end
