function files = replot_accuracy_pdf(runDir, objects)
%REPLOT_ACCURACY_PDF Redraw PDF comparison figures of a finished accuracy run.
% Reloads the network, the fitted GMM/TMM/block-TMM models and the AC Monte
% Carlo reference recorded by run_tmm_transmission_accuracy, recomputes the
% selected models' distributions, and draws the requested objects into
% runDir/pdf_selected. objects: struct array with quantity, column, file, label.
arguments
    runDir (1,1) string
    objects (1,:) struct
end
projectRoot = string(fileparts(fileparts(mfilename("fullpath"))));
addpath(char(projectRoot));
record = jsondecode(fileread(fullfile(runDir, "run_manifest.json")));
config = record.config;
TMM_Transmission_utils.io.setup_external_dependencies(config);
mpc = loadcase(char(config.networkFile));
if string(config.reactiveLimitRule) ~= "none"
    mpc = TMM_Transmission_core.powerflow.release_slack_reactive_limits(mpc);
end
control = TMM_Transmission_core.powerflow.describe_bus_control(mpc);
saved = load(fullfile(runDir, "fitted_models.mat"), "fits", "blockFits", "slackProfile");
reference = load(fullfile(config.referenceReusedFrom, "states.mat"), "reference").reference;
observed = TMM_core.network.evaluate_outputs(reference.voltage, mpc);
observed.voltageAngle = rad2deg(angle(reference.voltage(:, control.nonReference)));
summary = jsondecode(fileread(fullfile(runDir, "summary.json")));
plotted = struct([]);
for method = ["GMM", "TMM", "blockTMM"]
    k = summary.selectedK.(method);
    index = find(saved.fits.componentCounts == k);
    if method == "blockTMM"
        model = saved.blockFits.models{index}; name = "分块TMM"; style = "block_tmm";
    else
        model = saved.fits.pairs{index}.(lower(method)); name = method; style = "proposed";
        if method == "GMM", style = "gmm"; end
    end
    prediction = TMM_Transmission_core.distributions.propagate_voltage_model(model, saved.fits.transform, mpc, ...
        saved.slackProfile, config.distribution);
    entry = struct("name", name, "style", style, "K", k, "prediction", prediction);
    if isempty(plotted), plotted = entry; else, plotted(end+1) = entry; end %#ok<AGROW>
end
files = TMM_Transmission_utils.plotting.export_pdf_comparison(observed, plotted, objects, fullfile(runDir, "pdf_selected"));
end
