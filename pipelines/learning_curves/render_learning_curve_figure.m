function render_learning_curve_figure(runRoot)
%RENDER_LEARNING_CURVE_FIGURE Render Supplementary Figure S2.
%
% Same data as primary_learning_curve_positive_r, restyled for the
% supplement: no figure title and no in-figure note (both live in Word),
% abbreviation-only panel labels, and a physical figure size matching the
% 6.27 inch text width so 300 dpi output carries true print point sizes.

assert(nargin >= 1 && ~isempty(runRoot), ...
    'CPM:FigureS2:RunRootRequired', ...
    'Supply the learning-curve run directory.');
outputPath = fullfile(runRoot, 'output');
summaryCsv = fullfile(outputPath, 'primary_learning_curve_summary.csv');
assert(isfile(summaryCsv), 'CPM:FigureS2:MissingSummary', ...
    'Missing summary: %s', summaryCsv);

TEXT_WIDTH_INCHES = 6.73;   % tight-crop lands this on the 6.27in text width
FIGURE_HEIGHT_INCHES = 4.15;
RESOLUTION_DPI = 300;
OUTCOMES = ["GDE", "MEQ30", "VRS", "OBN", "BDE", "AED"];
lineColor = [0.85 0.33 0.10];

summary = readtable(summaryCsv, 'TextType', 'string');
assert(all(ismember(OUTCOMES, unique(summary.Outcome))), ...
    'CPM:FigureS2:SummaryRows', ...
    'The summary must contain all six outcomes.');

span = [min(summary.RPositiveP2_5), max(summary.RPositiveP97_5)];
padding = 0.05 * diff(span);
sharedYLimits = [floor((span(1) - padding) * 10) / 10, ...
    ceil((span(2) + padding) * 10) / 10];

figureHandle = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'inches', ...
    'Position', [1 1 TEXT_WIDTH_INCHES FIGURE_HEIGHT_INCHES]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for outcomeIndex = 1:numel(OUTCOMES)
    outcome = OUTCOMES(outcomeIndex);
    axisHandle = nexttile(layout);
    set(axisHandle, 'Color', 'white', 'XColor', 'black', ...
        'YColor', 'black', 'GridColor', [0.75 0.75 0.75], ...
        'GridAlpha', 0.45, 'FontSize', 7, 'LineWidth', 0.5);
    hold(axisHandle, 'on');
    grid(axisHandle, 'on');

    current = sortrows(summary(summary.Outcome == outcome, :), 'SampleN');
    x = current.SampleN';
    fill(axisHandle, [x fliplr(x)], ...
        [current.RPositiveP2_5' fliplr(current.RPositiveP97_5')], ...
        lineColor, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    yline(axisHandle, 0, ':', 'Color', [0.35 0.35 0.35]);
    plot(axisHandle, current.SampleN, current.MedianRPositive, ...
        'Color', lineColor, 'LineStyle', '-', 'LineWidth', 1.2, ...
        'Marker', 'o', 'MarkerSize', 3.2, 'MarkerFaceColor', 'white');

    title(axisHandle, outcome, 'FontWeight', 'bold', 'Color', 'black', ...
        'FontSize', 8.5);
    xlabel(axisHandle, 'Participants included (N)', 'Color', 'black', ...
        'FontSize', 7.5);
    ylabel(axisHandle, 'Prediction accuracy (r)', 'Color', 'black', ...
        'FontSize', 7.5);
    xlim(axisHandle, [18 69]);
    xticks(axisHandle, [20 30 40 50 60 67]);
    ylim(axisHandle, sharedYLimits);
end

pngPath = fullfile(outputPath, 'figure_s2_learning_curves.png');
pdfPath = fullfile(outputPath, 'figure_s2_learning_curves.pdf');
exportgraphics(figureHandle, pngPath, 'Resolution', RESOLUTION_DPI, ...
    'BackgroundColor', 'white');
exportgraphics(figureHandle, pdfPath, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
clear cleanup

info = imfinfo(pngPath);
fprintf('Figure S2 written: %s\n', pngPath);
fprintf('  pixels: %d x %d\n', info.Width, info.Height);
fprintf('  target print size: %.2f x %.2f inches at %d dpi\n', ...
    info.Width / RESOLUTION_DPI, info.Height / RESOLUTION_DPI, ...
    RESOLUTION_DPI);
fprintf('  shared y-limits: [%.1f %.1f]\n', sharedYLimits);
fprintf('  vector copy: %s\n', pdfPath);
end
