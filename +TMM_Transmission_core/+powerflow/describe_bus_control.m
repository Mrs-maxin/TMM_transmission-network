function control = describe_bus_control(mpc)
%DESCRIBE_BUS_CONTROL Slack, PV and PQ sets with scheduled generation.
% The slack bus holds its voltage and balances active power and losses; PV
% buses hold scheduled active output and voltage set point; PQ buses carry
% fixed net injections (derivation eq. (4)). Units fixed at a reactive limit
% are converted to PQ by apply_reactive_limit_mode (derivation Section 2.3).
% Bus order of MATPOWER is kept.
arguments
    mpc (1,1) struct
end
define_constants;
n = size(mpc.bus, 1);
assert(isequal(mpc.bus(:, BUS_I), (1:n).'), "TMM_Transmission:powerflow:BusNumbering", ...
    "Contiguous MATPOWER bus numbers 1..n are required.");
online = mpc.gen(:, GEN_STATUS) > 0;
genBus = mpc.gen(online, GEN_BUS);
types = mpc.bus(:, BUS_TYPE);
slack = find(types == REF);
assert(isscalar(slack), "TMM_Transmission:powerflow:SlackBus", "Exactly one slack bus is required.");
pv = find(types == PV);
pq = find(types == PQ);
assert(all(ismember(pv, genBus)), "TMM_Transmission:powerflow:PvWithoutGenerator", ...
    "Every PV bus needs an online generator.");
nonReference = setdiff((1:n).', slack, "stable");
scheduledP = accumarray(genBus, mpc.gen(online, PG), [n, 1]);
scheduledQ = accumarray(genBus, mpc.gen(online, QG), [n, 1]);
voltageSet = mpc.bus(:, VM);
for bus = [slack; pv].'
    setPoints = mpc.gen(online & mpc.gen(:, GEN_BUS) == bus, VG);
    assert(all(abs(setPoints - setPoints(1)) < 1e-12), "TMM_Transmission:powerflow:VoltageSetPoint", ...
        "Generators at bus %d have different voltage set points.", bus);
    voltageSet(bus) = setPoints(1);
end
[~, pqPosition] = ismember(pq, nonReference);
[~, pvPosition] = ismember(pv, nonReference);
control = struct("busCount", n, "baseMVA", mpc.baseMVA, "slack", slack, "pv", pv, "pq", pq, ...
    "nonReference", nonReference, "pqPosition", pqPosition, "pvPosition", pvPosition, ...
    "scheduledP", scheduledP, "scheduledQ", scheduledQ, "voltageSet", voltageSet, ...
    "reactiveLimits", "not-applied");
end
