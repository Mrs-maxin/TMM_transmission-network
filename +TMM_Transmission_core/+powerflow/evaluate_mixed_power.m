function power = evaluate_mixed_power(state, anchor)
%EVALUATE_MIXED_POWER Exact rectangular AC rows for PQ and PV buses in batch.
% STATE holds [real; imaginary] voltages of all non-slack buses per column.
% Output rows follow build_mixed_power_targets (derivation eq. (5)).
arguments
    state (:,:) double
    anchor (1,1) struct
end
m = numel(anchor.control.nonReference);
assert(size(state, 1) == 2*m, "TMM_Transmission:powerflow:StateDimension", ...
    "State must have 2*(n-1) rows.");
voltage = repmat(anchor.voltage, 1, size(state, 2));
voltage(anchor.control.nonReference, :) = complex(state(1:m, :), state(m+1:end, :));
injection = voltage.*conj(anchor.Ybus*voltage);
power = [real(injection(anchor.control.nonReference, :)); imag(injection(anchor.control.pq, :)); ...
    abs(voltage(anchor.control.pv, :)).^2];
end
