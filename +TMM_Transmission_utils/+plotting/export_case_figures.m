% TEMPLATE SCAFFOLD: imported from project-template-lite; not recovered research implementation.
% Cloud Python and original project MATLAB files remain pending recovery.
function figureOutputs = export_case_figures(caseConfig, metrics)
%EXPORT_CASE_FIGURES Scaffold for formal case figure export.

TMM_Transmission_utils.io.ensure_dir(caseConfig.figureDir);

figureOutputs = struct();
figureOutputs.caseName = string(caseConfig.name);
figureOutputs.status = "scaffold";
figureOutputs.metricsStatus = string(metrics.status);
figureOutputs.figureDir = string(caseConfig.figureDir);
figureOutputs.files = strings(0, 1);
figureOutputs.warnings = "No formal figures are exported by the template scaffold.";
end
