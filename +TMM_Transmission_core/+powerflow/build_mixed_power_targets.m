function targets = build_mixed_power_targets(control, loadP, loadQ)
%BUILD_MIXED_POWER_TARGETS Specified quantities of the PQ/PV power-flow rows.
% Rows are net active injection at every non-slack bus, net reactive
% injection at PQ buses and squared voltage set point at PV buses, in p.u.
% Columns are load samples (derivation eqs. (4)-(5)).
arguments
    control (1,1) struct
    loadP (:,:) double
    loadQ (:,:) double
end
assert(isequal(size(loadP), size(loadQ)) && size(loadP, 2) == control.busCount, ...
    "TMM_Transmission:powerflow:LoadShape", "Load samples must be count-by-bus in MW and MVAr.");
count = size(loadP, 1);
netP = (control.scheduledP.' - loadP)/control.baseMVA;
netQ = (control.scheduledQ.' - loadQ)/control.baseMVA;
targets = [netP(:, control.nonReference).'; netQ(:, control.pq).'; ...
    repmat(control.voltageSet(control.pv).^2, 1, count)];
end
