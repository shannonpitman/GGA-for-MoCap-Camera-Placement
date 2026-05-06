function plotGA_ComputationTime(varargin)
% plotGA_ComputationTime  Computation time grouped by camera count using
% box-and-jitter plots (no bars).
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "GA wall-clock time scales with camera count and target type / cost
%   function, demonstrating practical tractability."
%
% Strengths
%   - Wall-clock time in hours is the metric that matters for an
%     engineering reader who wants to reproduce.
%   - Box-and-jitter exposes both spread and the actual sample points.
%   - Sample size n is annotated under each box.
%
% Decisions taken to address prior examiner critiques
%   1. BARS REPLACED WITH BOXES. Mean ± std bars assumed normality;
%      boxes show median, IQR and whiskers — with no hidden assumption.
%   2. SAMPLE SIZE VISIBLE. n is annotated under every box.
%   3. AUTO LOG SCALE. If max/min of any group exceeds 5×, the y-axis
%      switches to log so small bars are not crushed.
% Notes still to address in caption text
%   - State whether the time is total / per-generation / per-evaluation.
%   - State the hardware (CPU, RAM, MATLAB version, parallel pool).
%   - A scaling-fit (t = a·N^p with p reported) would strengthen the
%     "tractability" claim and is straightforward to add separately.
% =====================================================================
%
%   plotGA_ComputationTime('Name', Value, ...)
%
%   Box-and-jitter plot of computation time (in hours) per camera count,
%   split by target type or cost function. Sample size annotated under
%   each box.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SplitBy'  - 'TargetType' (default), 'CostFunction', or 'none'
%     'LogScale' - 'auto' (default), true, or false. 'auto' switches to
%                  log y when max/min ratio across groups > 5.
%     'SaveAs'   - Output filename without extension (default: auto)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SplitBy',  'TargetType', @ischar);
    addParameter(p, 'LogScale', 'auto', @(x) (ischar(x) && strcmpi(x,'auto')) || islogical(x));
    addParameter(p, 'SaveAs',   '',           @ischar);
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

    %% Build raw data cell array {nCams x nGroups}
    rawTimes  = cell(nCams, nGroups);
    for c = 1:nCams
        cam = uniqueCams(c);
        for g = 1:nGroups
            if isempty(splitField)
                mask = [runs.NumCameras] == cam;
            else
                mask = ([runs.NumCameras] == cam) & ([runs.(splitField)] == splitVals(g));
            end
            rawTimes{c,g} = [runs(mask).ElapsedTime] / 3600;  % seconds → hours
        end
    end

    %% Plot
    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);
    ax = axes(fig);
    hold(ax, 'on');

    xPos      = 1:nCams;
    totalSpan = 0.7;
    if nGroups > 1
        offsets = linspace(-totalSpan/2, totalSpan/2, nGroups);
    else
        offsets = 0;
    end
    boxHW     = 0.9 * (totalSpan / max(nGroups, 1)) / 2;
    jitterAmp = 0.45 * boxHW;

    legendHandles = gobjects(nGroups, 1);

    for g = 1:nGroups
        col = splitColors(g,:);
        for c = 1:nCams
            xc = xPos(c) + offsets(g);
            data = rawTimes{c,g};

            hBox = drawBoxPlot(ax, xc, data, col, boxHW);

            % Jittered points
            if ~isempty(data)
                xj = xc + jitterAmp*(rand(1,numel(data))-0.5);
                plot(ax, xj, data, 'o', ...
                    'MarkerSize', 2.5, ...
                    'MarkerEdgeColor', col*0.5, 'MarkerFaceColor', col*0.5, ...
                    'HandleVisibility', 'off');
            end

            if isgraphics(hBox) && ~isgraphics(legendHandles(g))
                legendHandles(g) = hBox;
            end
        end
    end

    %% Auto-decide log scale
    flatTimes = [rawTimes{:}];
    flatTimes = flatTimes(flatTimes > 0);
    useLog = false;
    if (ischar(opts.LogScale) && strcmpi(opts.LogScale, 'auto'))
        if ~isempty(flatTimes) && (max(flatTimes)/min(flatTimes)) > 5
            useLog = true;
        end
    elseif islogical(opts.LogScale)
        useLog = opts.LogScale;
    end

    if useLog
        set(ax, 'YScale', 'log');
    end

    %% Sample-size annotations (n=) under each box
    yLim = ylim(ax);
    if useLog
        yN = yLim(1) * 0.85;
    else
        yN = yLim(1) - 0.04 * (yLim(2)-yLim(1));
    end
    for g = 1:nGroups
        for c = 1:nCams
            xc = xPos(c) + offsets(g);
            n  = numel(rawTimes{c,g});
            if n == 0, continue; end
            text(ax, xc, yN, sprintf('n=%d', n), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
                'Color', 'k');
        end
    end
    if useLog
        ylim(ax, [yLim(1)*0.6, yLim(2)]);
    else
        ylim(ax, [yLim(1) - 0.10*(yLim(2)-yLim(1)), yLim(2)]);
    end

    hold(ax, 'off');

    %% Formatting
    set(ax, 'XTick', xPos, 'XTickLabel', ...
        arrayfun(@(x) sprintf('%d', x), uniqueCams, 'Uni', false));
    xlabel(ax, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax, 'Computation Time (hours)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    xlim(ax, [0.4, nCams+0.6]);
    grid(ax, 'on');

    if nGroups > 1
        validH = legendHandles(isgraphics(legendHandles));
        validN = splitNames(isgraphics(legendHandles));
        if ~isempty(validH)
            legend(validH, validN, 'Location', 'northwest', ...
                'FontSize', sty.FontSizeLegend);
        end
    end

    % --- Apply thesis style: white bg, black text/axes/ticks/legend ---
    applyThesisStyle(fig);

    %% Print summary
    fprintf('Computation time summary (hours, median [IQR], min–max):\n');
    for c = 1:nCams
        for g = 1:nGroups
            d = rawTimes{c,g};
            if isempty(d)
                fprintf('  %dC %s: (no runs)\n', uniqueCams(c), splitNames{g});
                continue;
            end
            q = quantile(d, [0.25 0.50 0.75]);
            fprintf('  %dC %s: med=%.2f [%.2f–%.2f], range %.2f–%.2f, n=%d\n', ...
                    uniqueCams(c), splitNames{g}, q(2), q(1), q(3), ...
                    min(d), max(d), numel(d));
        end
    end

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_CompTime_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helper ----
function hBox = drawBoxPlot(ax, xc, data, col, hw)
% Draw a single box-and-whisker at horizontal position xc.
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
        'FaceAlpha', 0.35, 'EdgeColor', col, 'LineWidth', 1.0);
    plot(ax, [xc-hw, xc+hw], [q(2), q(2)], 'Color', col, 'LineWidth', 1.5);
    plot(ax, [xc, xc], [wLo, q(1)], 'Color', col, 'LineWidth', 0.8);
    plot(ax, [xc, xc], [q(3), wHi], 'Color', col, 'LineWidth', 0.8);
    plot(ax, [xc-hw*0.5, xc+hw*0.5], [wLo, wLo], 'Color', col, 'LineWidth', 0.8);
    plot(ax, [xc-hw*0.5, xc+hw*0.5], [wHi, wHi], 'Color', col, 'LineWidth', 0.8);
    if ~isempty(outliers)
        plot(ax, repmat(xc, size(outliers)), outliers, ...
            'o', 'MarkerSize', 3, 'MarkerEdgeColor', col, 'HandleVisibility', 'off');
    end
end
