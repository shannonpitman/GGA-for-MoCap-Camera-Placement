function plotGA_FactorEffects(varargin)
% plotGA_FactorEffects  Side-by-side comparison of factor effects, with
% sample-size annotations and Mann-Whitney significance brackets.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Both target type (UAV/UGV) and grid discretisation (Uniform/Normal)
%   shift the cost-vs-camera-count relationship; the magnitudes
%   identify which factor matters more."
%
% Strengths
%   - Side-by-side panels enforce the same camera-count axis, making
%     direct comparison easy.
%   - Y-axes are linked so cross-panel magnitudes are visually comparable.
%   - Each box now carries its sample size n.
%   - Each pair of boxes (UAV vs UGV, Uniform vs Normal) carries a
%     Mann-Whitney U-test p-value and APA-style stars, so the visual
%     factor claim has statistical backing.
%
% Decisions taken to address prior examiner critiques
%   1. LINKED Y-AXES enforced in this version.
%   2. SIGNIFICANCE TESTING per camera count per panel (Mann-Whitney).
%   3. n= ANNOTATIONS above each box.
% Notes still belonging in the caption / discussion text
%   - State the cost function this plot was filtered to.
%   - Cite a 2- or 3-way ANOVA on cost ~ Cameras + TargetType +
%     GridMode if you want a single omnibus test rather than per-cam
%     pairwise tests.
%   - Interaction terms (TargetType × GridMode) are not shown here;
%     either argue they were tested and negligible, or add a third
%     interaction-plot panel.
% =====================================================================
%
%   plotGA_FactorEffects('Name', Value, ...)
%
%   Produces a 1x2 figure showing:
%     Left panel:  UAV vs UGV (target type effect)
%     Right panel: Uniform vs Normal grid (discretisation effect)
%
%   Each panel shows grouped box plots by camera count with n=
%   annotations and Mann-Whitney significance brackets.
%   Best used when filtered to a single cost function.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'ShowStats' - true to overlay Mann-Whitney brackets (default: true)
%     'SaveAs'    - Output filename without extension (default: auto)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'ShowStats', true, @islogical);
    addParameter(p, 'SaveAs',    '',   @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty   = gaPlotStyle();
    stats = gaStatsHelpers();

    uniqueCams = sort(unique([runs.NumCameras]));

    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    %% Panel 1: Target Type
    ax1 = subplot(1, 2, 1, 'Parent', fig);
    plotGroupedBoxes(ax1, runs, uniqueCams, 'TargetType', ...
        sty.TargetColors, sty.TargetNames, sty, stats, opts.ShowStats);
    xlabel(ax1, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax1, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);

    %% Panel 2: Grid Mode
    ax2 = subplot(1, 2, 2, 'Parent', fig);
    plotGroupedBoxes(ax2, runs, uniqueCams, 'GridMode', ...
        sty.GridColors, sty.GridNames, sty, stats, opts.ShowStats);
    xlabel(ax2, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax2, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);

    % Link y-axes so cross-panel magnitudes are visually comparable.
    linkaxes([ax1, ax2], 'y');

    if opts.ShowStats
        sgtitle(fig, ...
            'Factor effects  (significance: Mann–Whitney U; *<0.05, **<0.01, ***<0.001)', ...
            'FontSize', sty.FontSizeAxis + 1, 'FontWeight', 'normal', ...
            'FontName', sty.FontName);
    end

    % --- Apply thesis style: white bg, black text/axes/ticks/legend ---
    applyThesisStyle(fig);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_FactorEffects_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helper ----
function plotGroupedBoxes(ax, runs, uniqueCams, fieldName, colors, names, sty, stats, showStats)
% Draw simple grouped box plots on the given axes, with n= and optional
% Mann-Whitney significance brackets between two-level splits.

    hold(ax, 'on');
    nCams = length(uniqueCams);

    vals = [runs.(fieldName)];
    uniqueVals = sort(unique(vals(~isnan(vals))));
    nGroups = length(uniqueVals);
    if nGroups < 1
        hold(ax, 'off');
        return;
    end

    grpOffset = linspace(-0.2, 0.2, nGroups);
    boxHW = 0.12;

    legendHandles = gobjects(nGroups, 1);
    cellData = cell(nCams, nGroups);

    for c = 1:nCams
        cam = uniqueCams(c);
        for g = 1:nGroups
            mask = ([runs.NumCameras] == cam) & (vals == uniqueVals(g));
            costs = [runs(mask).BestCost];
            cellData{c,g} = costs;
            if isempty(costs), continue; end

            xc  = c + grpOffset(g);
            col = colors(min(g, size(colors,1)), :);

            hB = drawBoxPlotLocal(ax, xc, costs, col, boxHW);

            if isgraphics(hB) && ~isgraphics(legendHandles(g))
                legendHandles(g) = hB;
            end
        end
    end

    set(ax, 'XTick', 1:nCams, ...
        'XTickLabel', arrayfun(@(x) sprintf('%d', x), uniqueCams, 'Uni', false));
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    xlim(ax, [0.4, nCams+0.6]);
    grid(ax, 'on');
    drawnow;

    yLim = ylim(ax);
    yRange = yLim(2) - yLim(1);
    if yRange == 0, yRange = max(abs(yLim(2)),1); end

    % n= annotations
    for c = 1:nCams
        for g = 1:nGroups
            data = cellData{c,g};
            if isempty(data), continue; end
            xc = c + grpOffset(g);
            yTop = max(data) + 0.025*yRange;
            text(ax, xc, yTop, sprintf('n=%d', numel(data)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
                'Color', 'k');
        end
    end

    % Significance brackets (only for two-level splits)
    if showStats && nGroups == 2
        for c = 1:nCams
            d1 = cellData{c,1};
            d2 = cellData{c,2};
            if isempty(d1) || isempty(d2), continue; end
            pVal = stats.mannWhitney(d1, d2);
            lbl  = sprintf('%s (p=%.3f)', stats.sigStars(pVal), pVal);

            x1 = c + grpOffset(1);
            x2 = c + grpOffset(2);
            yLine = max([d1, d2]) + 0.10*yRange;
            stats.drawSigBracket(ax, x1, x2, yLine, lbl, sty.FontSizeAnnot);
        end
    end

    ylim(ax, [yLim(1), yLim(2) + 0.20*yRange]);

    validH = legendHandles(isgraphics(legendHandles));
    validN = names(isgraphics(legendHandles));
    if ~isempty(validH)
        legend(ax, validH, validN, 'Location', 'best', 'FontSize', sty.FontSizeLegend);
    end

    hold(ax, 'off');
end


function hBox = drawBoxPlotLocal(ax, xc, data, col, hw)
    hBox = gobjects(0);
    data = data(~isnan(data));
    if isempty(data), return; end

    if numel(data) == 1
        hBox = plot(ax, xc, data, 's', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 6);
        return;
    end

    q = quantile(data, [0.25, 0.50, 0.75]);
    iqr_val = q(3) - q(1);
    wLo = max(min(data), q(1) - 1.5*iqr_val);
    wHi = min(max(data), q(3) + 1.5*iqr_val);
    outliers = data(data < wLo | data > wHi);

    hBox = fill(ax, [xc-hw, xc+hw, xc+hw, xc-hw], ...
        [q(1), q(1), q(3), q(3)], col, ...
        'FaceAlpha', 0.35, 'EdgeColor', col, 'LineWidth', 0.8);
    plot(ax, [xc-hw, xc+hw], [q(2), q(2)], 'Color', col, 'LineWidth', 1.5);
    plot(ax, [xc, xc], [wLo, q(1)], 'Color', col, 'LineWidth', 0.6);
    plot(ax, [xc, xc], [q(3), wHi], 'Color', col, 'LineWidth', 0.6);
    plot(ax, [xc-hw*0.5, xc+hw*0.5], [wLo, wLo], 'Color', col, 'LineWidth', 0.6);
    plot(ax, [xc-hw*0.5, xc+hw*0.5], [wHi, wHi], 'Color', col, 'LineWidth', 0.6);
    if ~isempty(outliers)
        plot(ax, repmat(xc, size(outliers)), outliers, ...
            'o', 'MarkerSize', 3, 'MarkerEdgeColor', col, 'HandleVisibility', 'off');
    end
end
