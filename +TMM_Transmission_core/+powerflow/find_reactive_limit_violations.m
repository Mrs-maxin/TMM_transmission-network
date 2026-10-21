function [newState, generatorQ, violated] = find_reactive_limit_violations(modeMpc, modeState, voltage, loadQ, tolerance)
%FIND_REACTIVE_LIMIT_VIOLATIONS Regulating PV units outside reactive limits.
% Unit reactive output is the bus reactive injection plus local reactive load
% (MATPOWER pfsoln, one unit per regulated bus). All violating units are
% fixed together, as MATPOWER runpf does with pf.enforce_q_lims=1.
arguments
    modeMpc (1,1) struct
    modeState (1,:) double
    voltage (:,:) double
    loadQ (:,:) double
    tolerance (1,1) double {mustBePositive}
end
define_constants;
Ybus = makeYbus(modeMpc.baseMVA, modeMpc.bus, modeMpc.branch);
busQ = imag(voltage.*conj((Ybus*voltage.').'))*modeMpc.baseMVA + loadQ;
count = size(voltage, 1);
ng = size(modeMpc.gen, 1);
generatorQ = repmat(modeMpc.gen(:, QG).', count, 1);
online = modeMpc.gen(:, GEN_STATUS) > 0;
busType = modeMpc.bus(modeMpc.gen(:, GEN_BUS), BUS_TYPE);
regulating = find(online & (busType == PV | busType == REF)).';
for g = regulating
    bus = modeMpc.gen(g, GEN_BUS);
    assert(nnz(online & modeMpc.gen(:, GEN_BUS) == bus) == 1, "TMM_Transmission:powerflow:SharedBus", ...
        "Bus %d has several online units; no reactive sharing rule is defined.", bus);
    generatorQ(:, g) = busQ(:, bus);
end
upper = false(count, ng);
lower = false(count, ng);
upper(:, regulating) = generatorQ(:, regulating) > modeMpc.gen(regulating, QMAX).' + tolerance;
lower(:, regulating) = generatorQ(:, regulating) < modeMpc.gen(regulating, QMIN).' - tolerance;
newState = repmat(modeState, count, 1);
newState(upper) = 1;
newState(lower) = -1;
violated = any(upper | lower, 2);
end
