function mpc = release_slack_reactive_limits(mpc)
%RELEASE_SLACK_REACTIVE_LIMITS Leave the slack unit without reactive limits.
% Conventional power flow keeps the slack bus as the voltage and angle
% reference and balances P and Q without limits; only PV units are checked
% against their reactive limits (derivation Section 2.3).
arguments
    mpc (1,1) struct
end
define_constants;
slackGenerators = mpc.gen(:, GEN_STATUS) > 0 & mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE) == REF;
mpc.gen(slackGenerators, QMAX) = Inf;
mpc.gen(slackGenerators, QMIN) = -Inf;
end
