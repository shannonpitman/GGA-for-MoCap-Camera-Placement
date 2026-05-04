function plotGA_ComputationTime(varargin)
% plotGA_ComputationTime  Computation time grouped by camera count.
%   plotGA_ComputationTime('Name', Value, ...)
%
%   Bar chart showing mean computation time (in hours) per camera count,
%   split by target type. Error bars show ± 1 std. Individual run times
%   overlaid as scatter points.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SplitBy' - 'TargetType' (default), 'CostFunction', or 'none'
%     'SaveAs'  - Output filename without extension (default: auto)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SplitBy', 'TargetType', @ischar);
    addParameter(p, 'SaveAs',  '',           @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty = gaPlotStyle();

    uniqueCams = sort(unique([runs.NumCameras]));
    nCams = length(uniqueCams);

    %% Determine split groups
    switch lower(opts.SplitBy)
        case 'targettype'
            splitField = 'TargetType';
            splitVals  = sort(unique([runs.TargetType]));
            splitVals  = splitVals(~isnan(splitVals));
            splitNames = arrayfun(@(v) sty.TargetNames{v}, splitVals, 'Uni', false);
            splitColors = sty.TargetColors(splitVals,:);
        case 'costfunction'
            splitField = 'CostFunctionType';
            splitVals  = sort(unique([runs.CostFunctionType]));
            splitNames = arrayfun(@(v) sty.CostFuncShort{v}, splitVals, 'Uni', false);
            splitColors = sty.CostFuncColors(splitVals,:);
        otherwise
            splitField = '';
            splitVals  = 1;
            splitNames = {'All'};
            splitColors = [0.4 0.6 0.8];
    end
    nGroups = length(splitVals);

    %% Build data matrix [nCams x nGroups] for means and stds
    meanTimes = zeros(nCams, nGroups);
    stdTimes  = zeros(nCams, nGroups);
    rawTimes  = cell(nCams, nGroups);

    for c = 1:nCams
        cam = uniqueCams(c);
        for g = 1:nGroups
            if isempty(splitField)
                mask = [runs.NumCameras] == cam;
            else
                mask = ([runs.NumCameras] == cam) & ([runs.(splitField)] == splitVals(g));
            end
            times_h = [runs(mask).ElapsedTime] / 3600;  % seconds -> hours
            rawTimes{c,g}  = times_h;
            meanTimes(c,g) = mean(times_h);
            stdTimes(c,g)  = std(times_h);
        end
    end

    %% Plot
    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', 'Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    xPos = 1:nCams;
    totalWidth = 0.7;
    barW = totalWidth / nGroups;
    offsets = linspace(-totalWidth/2 + barW/2, totalWidth/2 - barW/2, nGroups);

    hBars = gobjects(nGroups, 1);
    for g = 1:nGroups
        xg = xPos + offsets(g);
        hBars(g) = bar(ax, xg, meanTimes(:,g), barW, ...
            'FaceColor', splitColors(g,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');

        errorbar(ax, xg, meanTimes(:,g), stdTimes(:,g), ...
            'k.', 'LineWidth', 0.8, 'CapSize', 4);

        % Individual points
        for c = 1:nCams
            nPts = length(rawTimes{c,g});
            jx = xg(c) + 0.03*(rand(1,nPts)-0.5);
            plot(ax, jx, rawTimes{c,g}, 'o', ...
                'MarkerSize', 2.5, 'MarkerEdgeColor', splitColors(g,:)*0.5, ...
                'MarkerFaceColor', splitColors(g,:)*0.5, 'HandleVisibility', 'off');
        end
    end

    hold(ax, 'off');

    %% Formatting
    set(ax, 'XTick', xPos, 'XTickLabel', ...
        arrayfun(@(x) sprintf('%d', x), uniqueCams, 'Uni', false));
    xlabel(ax, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax, 'Computation Time (hours)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.15);

    if nGroups > 1
        legend(hBars, splitNames, 'Location', 'northwest', ...
            'FontSize', sty.FontSizeLegend);
    end

    %% Print summary
    fprintf('Computation time summary (hours):\n');
    for c = 1:nCams
        for g = 1:nGroups
            fprintf('  %dC %s: %.1f +/- %.1f h  (range: %.1f - %.1f)\n', ...
                uniqueCams(c), splitNames{g}, ...
                meanTimes(c,g), stdTimes(c,g), ...
                min(rawTimes{c,g}), max(rawTimes{c,g}));
        end
    end

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_CompTime_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', outName);
end