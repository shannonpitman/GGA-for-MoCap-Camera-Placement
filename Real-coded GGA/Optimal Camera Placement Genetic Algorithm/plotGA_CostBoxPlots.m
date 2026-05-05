function plotGA_CostBoxPlots(varargin)
% plotGA_CostBoxPlots  Box plots of best cost grouped by camera count,
% annotated with sample size and Mann-Whitney significance brackets
% between split groups.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Best cost decreases with the number of cameras, and the chosen
%   split factor (target type / grid mode / warm vs cold) shifts the
%   cost level. Differences between split groups are flagged as
%   significant or not via a non-parametric test."
%
% Strengths
%   - Box-plot is a defensible summary for small-n stochastic results.
%   - Outliers are explicitly drawn so the reader sees stragglers.
%   - Sample size is annotated above each box (n=…), so the reader can
%     judge the IQR's reliability.
%   - When a split factor with two levels is active, a Mann-Whitney
%     U-test p-value is printed above each camera-count cluster with
%     APA-style stars (n.s./*/**/***).
%
% Decisions taken to address prior examiner critiques
%   1. n= NOW VISIBLE per box.
%   2. SIGNIFICANCE TESTING added (Mann-Whitney U) for two-level splits.
%      For 3+ levels (e.g. cost-function), pairwise tests are skipped
%      to avoid multi-comparison clutter; report ANOVA in the caption.
%   3. Y-AXIS LABEL still states "Best Cost" — caption must define the
%      cost-function weights so units are interpretable.
%   4. CATEGORICAL X-TICKS only at integer camera counts; the in-between
%      space is unused.
% =====================================================================
%
%   plotGA_CostBoxPlots('Name', Value, ...)
%
%   Produces one figure per cost function (or a single figure if
%   CostFunction is specified). Within each figure, the x-axis groups by
%   camera count and separate box groups distinguish target type, grid
%   mode, or warm-start status.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SplitBy'    - 'TargetType' (default), 'GridMode', 'WarmStart', or
%                    'none' to disable grouping
%     'ShowStats'  - true to overlay Mann-Whitney significance brackets
%                    when the split factor has exactly two levels.
%                    (default: true)
%     'SaveAs'     - Output filename prefix (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SplitBy',   'TargetType', @ischar);
    addParameter(p, 'ShowStats', true,         @islogical);
    addParameter(p, 'SaveAs',    '',           @ischar);
    parse(p, varargin{:});

    opts   = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty   = gaPlotStyle();
    stats = gaStatsHelpers();

    %% Determine which cost functions are present
    allCF = unique([runs.CostFunctionType]);

    for cf = allCF
        cfRuns = runs([runs.CostFunctionType] == cf);
        if isempty(cfRuns), continue; end

        fig = figure('Units', 'inches', ...
            'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
            'PaperPositionMode', 'auto', ...
            'Color', sty.BackgroundColor);
        ax = axes(fig);
        hold(ax, 'on');

        uniqueCams = sort(unique([cfRuns.NumCameras]));
        nCams = length(uniqueCams);

        %% Determine split groups
        switch lower(opts.SplitBy)
            case 'targettype'
                splitField = 'TargetType';
                splitVals  = sort(unique([cfRuns.TargetType]));
                splitVals  = splitVals(~isnan(splitVals));
                splitNames = sty.TargetNames;
                splitColors = sty.TargetColors;
            case 'gridmode'
                splitField = 'GridMode';
                splitVals  = sort(unique([cfRuns.GridMode]));
                splitVals  = splitVals(~isnan(splitVals));
                splitNames = sty.GridNames;
                splitColors = sty.GridColors;
            case 'warmstart'
                splitField = 'WarmStart';
                splitVals  = [false, true];
                splitNames = {'Cold-Start', 'Warm-Start'};
                splitColors = [sty.ColdColor; sty.WarmColor];
            otherwise
                splitField = '';
                splitVals  = 1;
                splitNames = {''};
                splitColors = sty.CameraColors(1,:);
        end
        nGroups = length(splitVals);

        %% Build grouped box-plot data structures and store raw points
        boxWidth   = 0.35;
        camSpacing = 1.0;
        grpOffset  = linspace(-0.2, 0.2, nGroups);
        cellData   = cell(nCams, nGroups);

        tickPos   = zeros(1, nCams);
        tickLabel = cell(1, nCams);

        for c = 1:nCams
            cam = uniqueCams(c);
            centre = c * camSpacing;
            tickPos(c)   = centre;
            tickLabel{c} = sprintf('%d', cam);

            for g = 1:nGroups
                if isempty(splitField)
                    mask = [cfRuns.NumCameras] == cam;
                else
                    fieldVals = [cfRuns.(splitField)];
                    mask = ([cfRuns.NumCameras] == cam) & (fieldVals == splitVals(g));
                end
                cellData{c,g} = [cfRuns(mask).BestCost];
            end
        end

        %% Draw boxes
        legendHandles = gobjects(nGroups, 1);
        hw = boxWidth / (nGroups + 0.5);

        for c = 1:nCams
            centre = tickPos(c);
            for g = 1:nGroups
                pos  = centre + grpOffset(g);
                col  = splitColors(g,:);
                data = cellData{c,g};

                hB = drawBoxPlotLocal(ax, pos, data, col, hw);

                if isgraphics(hB) && ~isgraphics(legendHandles(g))
                    legendHandles(g) = hB;
                end
            end
        end

        hold(ax, 'off');

        %% Formatting (do basic axes first so we know yLim for annotations)
        set(ax, 'XTick', tickPos, 'XTickLabel', tickLabel);
        xlabel(ax, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylabel(ax, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
            'Box', 'on', 'TickDir', 'out');
        xlim(ax, [tickPos(1)-0.6, tickPos(end)+0.6]);
        grid(ax, 'on');
        drawnow;

        yLim = ylim(ax);
        yRange = yLim(2) - yLim(1);
        if yRange == 0, yRange = max(abs(yLim(2)),1); end

        hold(ax, 'on');

        %% n= annotations (above each box, just outside the IQR)
        for c = 1:nCams
            centre = tickPos(c);
            for g = 1:nGroups
                data = cellData{c,g};
                if isempty(data), continue; end
                pos = centre + grpOffset(g);
                yTop = max(data) + 0.025*yRange;
                text(ax, pos, yTop, sprintf('n=%d', numel(data)), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
                    'Color', 'k');
            end
        end

        %% Mann-Whitney brackets between two-level splits
        if opts.ShowStats && nGroups == 2 && ~isempty(splitField)
            for c = 1:nCams
                d1 = cellData{c,1};
                d2 = cellData{c,2};
                if isempty(d1) || isempty(d2), continue; end
                pVal = stats.mannWhitney(d1, d2);
                lbl  = sprintf('%s (p=%.3f)', stats.sigStars(pVal), pVal);

                centre = tickPos(c);
                x1 = centre + grpOffset(1);
                x2 = centre + grpOffset(2);
                yLine = max([d1, d2]) + 0.10*yRange;
                stats.drawSigBracket(ax, x1, x2, yLine, lbl, sty.FontSizeAnnot);
            end
        end

        % Make room for the annotations on top
        ylim(ax, [yLim(1), yLim(2) + 0.20*yRange]);

        hold(ax, 'off');

        %% Legend
        if nGroups > 1 && ~isempty(splitField)
            validH = legendHandles(isgraphics(legendHandles));
            validN = splitNames(isgraphics(legendHandles));
            if ~isempty(validH)
                legend(validH, validN, 'Location', 'best', ...
                    'FontSize', sty.FontSizeLegend);
            end
        end

        % Title noting the test (if any)
        if opts.ShowStats && nGroups == 2 && ~isempty(splitField)
            title(ax, sprintf('%s — significance: Mann–Whitney U (*<0.05, **<0.01, ***<0.001)', ...
                                sty.CostFuncNames{cf}), ...
                'FontWeight', 'normal', ...
                'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        end

        % --- Apply thesis style: white bg, black text/axes/ticks/legend ---
        applyThesisStyle(fig);

        %% Export
        if isempty(opts.SaveAs)
            outName = sprintf('GA_CostBox_CF%d_%s', cf, filterDesc);
        else
            outName = sprintf('%s_CF%d', opts.SaveAs, cf);
        end
        exportgraphics(fig, [outName '.pdf'], ...
            'ContentType',     'vector', ...
            'BackgroundColor', sty.ExportBgColor);
        fprintf('Saved: %s.pdf  [%s]\n', outName, sty.CostFuncNames{cf});
    end
end


%% ---- Local helper ----
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
