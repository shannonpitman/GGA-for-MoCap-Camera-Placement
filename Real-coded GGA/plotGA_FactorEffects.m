function plotGA_FactorEffects(varargin)
% plotGA_FactorEffects  Side-by-side comparison of factor effects.
%   plotGA_FactorEffects('Name', Value, ...)
%
%   Produces a 1x2 figure showing:
%     Left panel:  UAV vs UGV (target type effect)
%     Right panel: Uniform vs Normal grid (discretisation effect)
%
%   Each panel shows grouped box plots by camera count.
%   Best used when filtered to a single cost function.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SaveAs' - Output filename without extension (default: auto)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SaveAs', '', @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty = gaPlotStyle();

    uniqueCams = sort(unique([runs.NumCameras]));
    nCams = length(uniqueCams);

    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', 'Color', 'w');

    %% Panel 1: Target Type
    ax1 = subplot(1, 2, 1, 'Parent', fig);
    plotGroupedBoxes(ax1, runs, uniqueCams, 'TargetType', ...
        sty.TargetColors, sty.TargetNames, sty);
    xlabel(ax1, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax1, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);

    %% Panel 2: Grid Mode
    ax2 = subplot(1, 2, 2, 'Parent', fig);
    plotGroupedBoxes(ax2, runs, uniqueCams, 'GridMode', ...
        sty.GridColors, sty.GridNames, sty);
    xlabel(ax2, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax2, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_FactorEffects_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helper ----
function plotGroupedBoxes(ax, runs, uniqueCams, fieldName, colors, names, sty)
% Draw simple grouped box plots on the given axes

    hold(ax, 'on');
    nCams = length(uniqueCams);

    vals = [runs.(fieldName)];
    uniqueVals = sort(unique(vals(~isnan(vals))));
    nGroups = length(uniqueVals);
    if nGroups < 1, return; end

    grpOffset = linspace(-0.2, 0.2, nGroups);
    boxHW = 0.12;

    legendHandles = gobjects(nGroups, 1);

    for c = 1:nCams
        cam = uniqueCams(c);
        for g = 1:nGroups
            mask = ([runs.NumCameras] == cam) & (vals == uniqueVals(g));
            costs = [runs(mask).BestCost];
            if isempty(costs), continue; end

            xc = c + grpOffset(g);
            col = colors(min(g, size(colors,1)), :);

            q = quantile(costs, [0.25, 0.50, 0.75]);
            iqr_val = q(3) - q(1);
            wLo = max(min(costs), q(1) - 1.5*iqr_val);
            wHi = min(max(costs), q(3) + 1.5*iqr_val);

            % Box
            hB = fill(ax, [xc-boxHW, xc+boxHW, xc+boxHW, xc-boxHW], ...
                [q(1), q(1), q(3), q(3)], col, ...
                'FaceAlpha', 0.35, 'EdgeColor', col, 'LineWidth', 0.8);
            % Median
            plot(ax, [xc-boxHW, xc+boxHW], [q(2), q(2)], 'Color', col, 'LineWidth', 1.5);
            % Whiskers
            plot(ax, [xc, xc], [wLo, q(1)], 'Color', col, 'LineWidth', 0.6);
            plot(ax, [xc, xc], [q(3), wHi], 'Color', col, 'LineWidth', 0.6);

            if ~isgraphics(legendHandles(g), 'matlab.graphics.primitive.Patch')
                legendHandles(g) = hB;
            end
        end
    end

    hold(ax, 'off');

    set(ax, 'XTick', 1:nCams, ...
        'XTickLabel', arrayfun(@(x) sprintf('%d', x), uniqueCams, 'Uni', false));
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    xlim(ax, [0.4, nCams+0.6]);
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.15);

    validH = legendHandles(isgraphics(legendHandles));
    validN = names(isgraphics(legendHandles));
    if ~isempty(validH)
        legend(ax, validH, validN, 'Location', 'best', 'FontSize', sty.FontSizeLegend);
    end
end