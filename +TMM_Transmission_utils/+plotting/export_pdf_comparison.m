function files = export_pdf_comparison(observed, predictions, objects, outDir)
%EXPORT_PDF_COMPARISON Diagnostic PDF figures: AC Monte Carlo histogram with model densities.
% observed: AC Monte Carlo outputs (fields as evaluate_outputs plus voltageAngle).
% predictions: struct array with fields name, style, K, prediction.
% objects: struct array with fields quantity, column, label, file.
% Model curves are weighted Gaussian kernels with the AC-sample bandwidth, as
% in AQPPF export_aqppf_comparison; they are for display only. Diagnostic
% figures, not formal manuscript figures.
arguments
    observed (1,1) struct
    predictions (1,:) struct
    objects (1,:) struct
    outDir (1,1) string
end
TMM_Transmission_utils.io.ensure_dir(outDir);
files = strings(0, 1);
for object = objects
    a = observed.(object.quantity)(:, object.column);
    lo = min(a); hi = max(a);
    for j = 1:numel(predictions)
        values = predictions(j).prediction.(object.quantity)(:, object.column);
        lo = min(lo, prctile(values, 0.1)); hi = max(hi, prctile(values, 99.9));
    end
    pad = 0.08*(hi - lo);
    support = linspace(lo - pad, hi + pad, 400).';
    [~, ~, bandwidth] = ksdensity(a);
    bandwidth = max(bandwidth, (hi - lo)/200);
    fig = TMM_Transmission_utils.plotting.create_figure("single-column", Visible="off");
    ax = axes(fig); hold(ax, "on"); box(ax, "on"); grid(ax, "on");
    set(ax, "FontName", "PingFang SC", "FontSize", 9);
    h = histogram(ax, a, linspace(lo, hi, 41), "Normalization", "pdf", "FaceColor", [0.82 0.82 0.82], ...
        "EdgeColor", "none", "DisplayName", sprintf("交流MC（%d组）", numel(a)));
    curves = zeros(numel(support), numel(predictions));
    names = strings(1, numel(predictions));
    for j = 1:numel(predictions)
        p = predictions(j).prediction;
        values = p.(object.quantity)(:, object.column);
        w = p.sampleWeights(:)/sum(p.sampleWeights);
        curves(:, j) = ksdensity(values, support, "Bandwidth", bandwidth, "Weights", w);
        ks = TMM_core.metrics.compare_distributions(struct("value", a), ...
            struct("value", values, "sampleWeights", w), predictions(j).name).KS;
        TMM_Transmission_utils.plotting.plot_method_line(ax, support, curves(:, j), predictions(j).style, ...
            "DisplayName", sprintf("%s（K=%d），KS=%.3f", predictions(j).name, predictions(j).K, ks), "LineWidth", 1.2);
        names(j) = matlab.lang.makeValidName(predictions(j).name + "Density");
    end
    xlim(ax, [support(1), support(end)]);
    ylim(ax, [0, 1.5*max([h.Values(:); curves(:)])]);
    xlabel(ax, object.label); ylabel(ax, "概率密度");
    legend(ax, "Location", "northeast", "FontSize", 7, "FontName", "PingFang SC");
    base = fullfile(outDir, object.file + "_PDF");
    bins = table(h.BinEdges(1:end-1).', h.BinEdges(2:end).', h.Values.', ...
        VariableNames=["leftEdge", "rightEdge", "acDensity"]);
    exportgraphics(fig, base + ".png", "Resolution", 300);
    savefig(fig, base + ".fig");
    close(fig);
    writetable([table(support, VariableNames="value"), array2table(curves, VariableNames=names)], base + "_plot_data.csv");
    writetable(bins, base + "_histogram.csv");
    files(end+1) = base + ".png"; %#ok<AGROW>
end
end
