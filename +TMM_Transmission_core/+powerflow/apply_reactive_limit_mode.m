function mpc = apply_reactive_limit_mode(mpc, modeState)
%APPLY_REACTIVE_LIMIT_MODE Network data for one reactive-limit operating mode.
% modeState(g) = 0 keeps unit g regulating its terminal voltage, +1 fixes its
% reactive output at Qmax and -1 at Qmin; a fixed unit's bus becomes PQ and its
% voltage is released (MATPOWER runpf rule with pf.enforce_q_lims=1).
arguments
    mpc (1,1) struct
    modeState (1,:) double
end
define_constants;
assert(numel(modeState) == size(mpc.gen, 1) && all(ismember(modeState, [-1, 0, 1])), ...
    "TMM_Transmission:powerflow:ModeState", "modeState needs one entry in {-1,0,1} per generator.");
fixed = find(modeState ~= 0);
for g = fixed
    bus = mpc.gen(g, GEN_BUS);
    assert(mpc.bus(bus, BUS_TYPE) == PV, "TMM_Transmission:powerflow:FixedNonPv", ...
        "Generator %d at bus %d is not a PV unit and cannot be fixed at a reactive limit.", g, bus);
    sameBus = find(mpc.gen(:, GEN_STATUS) > 0 & mpc.gen(:, GEN_BUS) == bus);
    assert(all(modeState(sameBus) ~= 0), "TMM_Transmission:powerflow:SharedBus", ...
        "All online units at bus %d must be fixed together; no sharing rule is defined.", bus);
    if modeState(g) > 0
        mpc.gen(g, QG) = mpc.gen(g, QMAX);
    else
        mpc.gen(g, QG) = mpc.gen(g, QMIN);
    end
    mpc.bus(bus, BUS_TYPE) = PQ;
end
end
