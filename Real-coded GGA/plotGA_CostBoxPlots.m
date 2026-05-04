function plotGA_CostBoxPlots(varargin)
% plotGA_CostBoxPlots  Box plots of best cost grouped by camera count.
%   plotGA_CostBoxPlots('Name', Value, ...)
%
%   Produces one figure per cost function (or a single figure if
%   CostFunction is specified). Within each figure, the x-axis groups by
%   camera count and separate box groups distinguish target type (UAV/UGV).
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SplitBy'    - 'TargetType' (default), 'GridMode', 'WarmStart', or
%                    'none' to disable grouping
%     'SaveAs'     - Output filename prefix (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SplitBy', 'TargetType', @ischar);
    addParameter(p, 'SaveAs',  '',           @ischar);
    parse(p, varargin{:});

    opts   = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty = gaPlotStyle();

    %% Determine which cost functions are present
    allCF = unique([runs.CostFunctionType]);

    for cf = allCF
        cfRuns = runs([runs.CostFunctionType] == cf);
        if isempty(cfRuns), continue; end

        fig = figure('Units', 'inches', ...
            'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
            'PaperPositionMode', 'auto', 'Color', 'w');
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

        %% Build grouped box plot data
        allData   = [];
        groupPos  = [];
        groupCol  = [];
        tickPos   = [];
        tickLabel = {};

        boxWidth  = 0.35;
        camSpacing = 1.0;
        grpOffset = linspace(-0.2, 0.2, nGroups);

        for c = 1:nCams
            cam = uniqueCams(c);
            centre = c * camSpacing;
            tickPos(end+1)   = centre; %#ok<AGROW>
            tickLabel{end+1} = sprintf('%d', cam); %#ok<AGROW>

            for g = 1:nGroups
                if isempty(splitField)
                    mask = [cfRuns.NumCameras] == cam;
                else
                    fieldVals = [cfRuns.(splitField)];
                    mask = ([cfRuns.NumCameras] == cam) & (fieldVals == splitVals(g));
                end

                costs = [cfRuns(mask).BestCost];
                nPts  = length(costs);
                pos   = centre + grpOffset(g);

                allData  = [allData;  costs(:)];          %#ok<AGROW>
                groupPos = [groupPos; repmat(pos, nPts, 1)]; %#ok<AGROW>
                groupCol = [groupCol; repmat(g,   nPts, 1)]; %#ok<AGROW>
            end
        end

        %% Draw box plots manually for colour control
        legendHandles = gobjects(nGroups, 1);

        for g = 1:nGroups
            gMask = groupCol == g;
            positions = unique(groupPos(gMask));

            for pos = positions'
                pMask = gMask & (groupPos == pos);
                data  = allData(pMask);
                if isempty(data), continue; end

                q = quantile(data, [0.25, 0.50, 0.75]);
                iqr_val = q(3) - q(1);
                wLo = max(min(data), q(1) - 1.5*iqr_val);
                wHi = min(max(data), q(3) + 1.5*iqr_val);
                outliers = data(data < wLo | data > wHi);

                hw = boxWidth / (nGroups + 0.5);
                col = splitColors(g,:);

                % Box
                hBox = fill(ax, ...
                    [pos-hw, pos+hw, pos+hw, pos-hw], ...
                    [q(1), q(1), q(3), q(3)], ...
                    col, 'FaceAlpha', 0.35, 'EdgeColor', col, ...
                    'LineWidth', 1.0);
                % Median line
                plot(ax, [pos-hw, pos+hw], [q(2), q(2)], ...
                    'Color', col, 'LineWidth', 1.5);
                % Whiskers
                plot(ax, [pos, pos], [wLo, q(1)], 'Color', col, 'LineWidth', 0.8);
                plot(ax, [pos, pos], [q(3), wHi], 'Color', col, 'LineWidth', 0.8);
                plot(ax, [pos-hw*0.5, pos+hw*0.5], [wLo, wLo], 'Color', col, 'LineWidth', 0.8);
                plot(ax, [pos-hw*0.5, pos+hw*0.5], [wHi, wHi], 'Color', col, 'LineWidth', 0.8);
                % Outliers
                if ~isempty(outliers)
                    plot(ax, repmat(pos, size(outliers)), outliers, ...
                        'o', 'MarkerSize', 3, 'MarkerEdgeColor', col);
                end

                % Stash first box handle per group for legend
                if isgraphics(legendHandles(g), 'matlab.graphics.primitive.Patch')
                    % already have one
                else
                    legendHandles(g) = hBox;
                end
            end
        end

        hold(ax, 'off');

        %% Formatting
        set(ax, 'XTick', tickPos, 'XTickLabel', tickLabel);
        xlabel(ax, 'Number of Cameras', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylabel(ax, 'Best Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
            'Box', 'on', 'TickDir', 'out');
        xlim(ax, [tickPos(1)-0.6, tickPos(end)+0.6]);
        grid(ax, 'on');
        set(ax, 'GridAlpha', 0.15);

        % Legend
        if nGroups > 1 && ~isempty(splitField)
            validH = legendHandles(isgraphics(legendHandles));
            validN = splitNames(isgraphics(legendHandles));
            if ~isempty(validH)
                legend(validH, validN, 'Location', 'best', ...
                    'FontSize', sty.FontSizeLegend);
            end
        end

        %% Export
        if isempty(opts.SaveAs)
            outName = sprintf('GA_CostBox_CF%d_%s', cf, filterDesc);
        else
            outName = sprintf('%s_CF%d', opts.SaveAs, cf);
        end
        exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
        fprintf('Saved: %s.pdf  [%s]\n', outName, sty.CostFuncNames{cf});
    end
end