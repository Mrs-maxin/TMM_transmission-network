% TEMPLATE SCAFFOLD: imported from project-template-lite; not recovered research implementation.
% Cloud Python and original project MATLAB files remain pending recovery.
function metrics = evaluate_case_metrics(caseConfig, proposedResults, baselineResults)
%EVALUATE_CASE_METRICS Scaffold for case metric aggregation.

metrics = struct();
metrics.caseName = string(caseConfig.name);
metrics.status = "scaffold";
metrics.proposedStatus = string(proposedResults.status);
metrics.baselineStatus = string(baselineResults.status);
metrics.metricNames = strings(0, 1);
metrics.files = struct();
metrics.warnings = "Replace +TMM_Transmission_core/+metrics/evaluate_case_metrics.m with project metrics.";
end
