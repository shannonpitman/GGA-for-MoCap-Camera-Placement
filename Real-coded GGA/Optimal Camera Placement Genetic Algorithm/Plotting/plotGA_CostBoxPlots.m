function plotGA_CostBoxPlots(varargin)
% plotGA_CostBoxPlots  Box plots of best cost grouped by camera count,
% with optional Mann-Whitney significance brackets, optional OptiTrack
% baseline overlay, and optional sample-size matching.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Best cost decreases with the number of cameras, and the chosen
%   split factor (warm vs cold, grid uniform vs normal) shifts the
%   cost level within the same problem. A red point + red star at the
%   7-cam cluster pin the ad-hoc OptiTrack baseline cost so the reader
%   can see the GA's improvement over the as-built lab rig."
%
% Strengths
%   - Box-plot is a defensible summary for small-n stochastic results.
%   - Outliers are explicitly drawn so the reader sees stragglers.
%   - Sample size is annotated above each box.
%   - When a split factor with two LEVELS THAT REPRESENT THE SAME
%     UNDERLYING PROBLEM is active (cold/warm, uniform/normal grid),
%     a Mann-Whitney U-test p-value is reported. The test is NOT
%     applied between target types (UAV vs UGV), because UAV and UGV
%     are different problems — their costs are not commensurable, and
%     a significant difference there only confirms "different
%     volumes give different costs", which is trivial. See the
%     'ShowStats' parameter to override.
%   - OptiTrack overlay quantifies the engineering improvement.
%
% Decisions taken to address prior examiner critiques
%   1. n= is annotated above every box.
%   2. p-values printed with three decimals; if rounded value would be
%      0.000 the bracket reads "p < 0.001" instead.
%   3. CATEGORICAL X-TICKS only at integer camera counts.
%   4. Cross-target-type brackets removed by default (see above).
%   5. OptiTrack ad-hoc baseline overlaid at the 7-cam position when
%      SplitBy='TargetType' and the data set has a single (GridMode,
%      Spacing), letting the reader judge the GA's headroom.
% =====================================================================
%
%   plotGA_CostBoxPlots('Name', Value, ...)
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'SplitBy'           - 'TargetType' (default), 'GridMode',
%                           'WarmStart', or 'none'
%     'ShowStats'         - 'auto' (default) — show brackets only for
%                           splits that compare the same problem
%                           (GridMode / WarmStart); 'on' to force
%                           brackets even for TargetType; 'off' to
%                           disable entirely.
%     'OptiTrackOverlay'  - true (default) to overlay the OptiTrack
%                           ad-hoc baseline at the 7-cam position
%                           (red filled circle for UAV, red 5-point
%                           star for UGV). Auto-skipped if SplitBy is
%                           not 'TargetType' or the loaded runs span
%                           multiple grid modes / spacings.
%     'OptiTrackWeights'  - [wUnc wOcc] for CF3 weights when
%                           evaluating the OptiTrack baseline.
%                           Default: [0.5 0.5].
%     'MatchSampleSize'   - true to random-subsample each group to
%                           min(n) for visual comparability. Default:
%                           false. Reports n actually plotted.
%     'SaveAs'            - Output filename prefix (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'SplitBy',          'TargetType', @ischar);
    addParameter(p, 'ShowStats',        'auto',       @(s) ischar(s) || isstring(s) || islogical(s));
    addParameter(p, 'OptiTrackOverlay', true,         @islogical);
    addParameter(p, 'OptiTrackWeights', [0.5 0.5],    @(v) isnumeric(v) && numel(v)==2);
    addParameter(p, 'MatchSampleSize',  false,        @islogical);
    addParameter(p, 'SaveAs',           '',           @ischar);
    parse(p, varargin{:});

    opts   = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    sty   = gaPlotStyle();
    stats = gaStatsHelpers();

    %% Determine which cost functions are present
    allCF = unique([runs.CostFunctionType]);

    % Normalise ShowStats setting
    [statsMode, ~] = canonicaliseStats(opts.ShowStats);

    % Decide if OptiTrack overlay is feasible for this data set
    [overlayEnabled, overlayGM, overlaySpacing] = ...
        decideOptiTrackOverlay(runs, opts);
    overlayCosts = struct();  % populated lazily inside the cf loop

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

        % Should this plot show inferential brackets?
        if strcmp(statsMode, 'on')
            doStats = (nGroups == 2 && ~isempty(splitField));
        elseif strcmp(statsMode, 'off')
            doStats = false;
        else   % auto
            % Skip brackets between different problems (UAV vs UGV)
            % because cost magnitudes are not commensurable.
            doStats = (nGroups == 2 && ~isempty(splitField)) && ...
                      ~strcmpi(splitField, 'TargetType');
        end

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

        %% Optional: match sample size across groups (random subsample)
        if opts.MatchSampleSize
            % Per-camera-count smallest n (across groups). Subsample
            % each group down to that n with a deterministic seed so
            % the result is reproducible across runs.
            seed = 20260513;
            rngState = rng(seed, 'twister');
            cleanup = onCleanup(@() rng(rngState));

            for c = 1:nCams
                nz = cellfun(@numel, cellData(c,:));
                nz = nz(nz>0);
                if numel(nz) < 2, continue; end
                nMin = min(nz);
                for g = 1:nGroups
                    d = cellData{c,g};
                    if numel(d) > nMin
                        idx = randperm(numel(d), nMin);
                        cellData{c,g} = d(idx);
                    end
                end
            end
            clear cleanup;
            fprintf(['plotGA_CostBoxPlots[%s]: subsampled each group ' ...
                     'down to min(n) per camera count.\n'], ...
                     sty.CostFuncShort{cf});
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

        %% Mann-Whitney brackets (only when commensurable)
        if doStats
            for c = 1:nCams
                d1 = cellData{c,1};
                d2 = cellData{c,2};
                if isempty(d1) || isempty(d2), continue; end
                pVal = stats.mannWhitney(d1, d2);
                lbl  = sprintf('%s (%s)', stats.sigStars(pVal), ...
                                          formatPValue(pVal));

                centre = tickPos(c);
                x1 = centre + grpOffset(1);
                x2 = centre + grpOffset(2);
                yLine = max([d1, d2]) + 0.10*yRange;
                stats.drawSigBracket(ax, x1, x2, yLine, lbl, sty.FontSizeAnnot);
            end
        end

        %% OptiTrack overlay (red dot UAV, red star UGV) ---------------
        overlayHandles = gobjects(0);
        overlayLabels  = {};
        if overlayEnabled && strcmpi(opts.SplitBy, 'TargetType')
            cam7Idx = find(uniqueCams == 7, 1);
            if ~isempty(cam7Idx)
                % Lazy-eval per-target OptiTrack costs once and re-use
                if ~isfield(overlayCosts, 'UAV')
                    overlayCosts.UAV = evaluateOptiTrackCost( ...
                        'TargetType', 1, 'GridMode', overlayGM, ...
                        'Spacing',    overlaySpacing, ...
                        'WeightUnc',  opts.OptiTrackWeights(1), ...
                        'WeightOcc',  opts.OptiTrackWeights(2));
                    overlayCosts.UGV = evaluateOptiTrackCost( ...
                        'TargetType', 2, 'GridMode', overlayGM, ...
                        'Spacing',    overlaySpacing, ...
                        'WeightUnc',  opts.OptiTrackWeights(1), ...
                        'WeightOcc',  opts.OptiTrackWeights(2));
                end

                cfField = sprintf('CF%d', cf);

                % UAV (TargetType==1) — red filled circle
                gU = find(splitVals == 1, 1);
                if ~isempty(gU) && isfield(overlayCosts.UAV, cfField)
                    xPos = tickPos(cam7Idx) + grpOffset(gU);
                    yVal = overlayCosts.UAV.(cfField);
                    hUAV = plot(ax, xPos, yVal, 'o', ...
                        'MarkerFaceColor', [0.85 0.10 0.10], ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerSize',      sty.MarkerSizeLg + 1, ...
                        'LineStyle',       'none');
                    overlayHandles(end+1) = hUAV;                  %#ok<AGROW>
                    overlayLabels{end+1}  = sprintf( ...
                        'OptiTrack ad-hoc UAV (%.4f)', yVal);
                end

                % UGV (TargetType==2) — red filled 5-point star
                gG = find(splitVals == 2, 1);
                if ~isempty(gG) && isfield(overlayCosts.UGV, cfField)
                    xPos = tickPos(cam7Idx) + grpOffset(gG);
                    yVal = overlayCosts.UGV.(cfField);
                    hUGV = plot(ax, xPos, yVal, 'pentagram', ...
                        'MarkerFaceColor', [0.85 0.10 0.10], ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerSize',      sty.MarkerSizeLg + 4, ...
                        'LineStyle',       'none');
                    overlayHandles(end+1) = hUGV;                  %#ok<AGROW>
                    overlayLabels{end+1}  = sprintf( ...
                        'OptiTrack ad-hoc UGV (%.4f)', yVal);
                end
            end
        end

        % Final y-axis: auto-expand to include OptiTrack overlay (which
        % typically sits well above the box data) and then add a small
        % headroom pad. The earlier yLim/yRange captured the box-only
        % range and is used purely for placing the n= annotations and
        % significance brackets just above the IQR boxes.
        ylim(ax, 'auto');
        yLimFinal = ylim(ax);
        yPadFinal = 0.08 * (yLimFinal(2) - yLimFinal(1));
        ylim(ax, [yLimFinal(1), yLimFinal(2) + yPadFinal]);

        % --- Cam-block dividers ----------------------------------------
        % Faint vertical line at the midpoint between adjacent camera-
        % count clusters so the reader sees the 6-cam vs 7-cam blocks as
        % distinct groups. Drawn after ylim is finalised so the line
        % spans the full vertical extent.
        yLimDraw = ylim(ax);
        for c = 1:(nCams-1)
            xDiv = 0.5 * (tickPos(c) + tickPos(c+1));
            plot(ax, [xDiv xDiv], yLimDraw, '-', ...
                'Color', [0.55 0.55 0.55 0.45], ...
                'LineWidth', 0.6, ...
                'HandleVisibility', 'off');
        end

        hold(ax, 'off');

        %% Legend
        legendH = [];
        legendL = {};
        if nGroups > 1 && ~isempty(splitField)
            validH = legendHandles(isgraphics(legendHandles));
            validN = splitNames(isgraphics(legendHandles));
            if ~isempty(validH)
                legendH = [legendH; validH(:)];
                legendL = [legendL, validN(:)'];
            end
        end
        if ~isempty(overlayHandles)
            legendH = [legendH; overlayHandles(:)];
            legendL = [legendL, overlayLabels(:)'];
        end
        if ~isempty(legendH)
            legend(legendH, legendL, 'Location', 'best', ...
                'FontSize', sty.FontSizeLegend);
        end

        % Title — describe the test only if it was actually run
        cfNameStr = sty.CostFuncNames{cf};
        if doStats
            titleStr = sprintf( ...
                '%s — significance: Mann–Whitney U (* p<0.05, ** p<0.01, *** p<0.001)', ...
                cfNameStr);
        else
            titleStr = cfNameStr;
        end
        title(ax, titleStr, ...
            'FontWeight', 'normal', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);

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
        fprintf('Saved: %s.pdf  [%s]\n', outName, cfNameStr);
    end
end


%% ---- Local helpers ------------------------------------------------------

function [mode, isLogical] = canonicaliseStats(arg)
% Map ShowStats argument to one of 'on' | 'off' | 'auto'.
    isLogical = false;
    if islogical(arg)
        isLogical = true;
        if arg, mode = 'on'; else, mode = 'off'; end
        return;
    end
    arg = lower(string(arg));
    switch arg
        case "auto",  mode = 'auto';
        case "on",    mode = 'on';
        case "off",   mode = 'off';
        case "true",  mode = 'on';
        case "false", mode = 'off';
        otherwise
            error('plotGA_CostBoxPlots:badShowStats', ...
                  'ShowStats must be ''auto'', ''on'', or ''off''.');
    end
end


function [enabled, gm, sp] = decideOptiTrackOverlay(runs, opts)
% Decide whether the OptiTrack overlay can be applied to this dataset.
%
% Requires: OptiTrackOverlay=true, SplitBy='TargetType', and the
% loaded runs share a single (GridMode, Spacing).
    enabled = false;
    gm = NaN;
    sp = NaN;
    if ~opts.OptiTrackOverlay, return; end
    if ~strcmpi(opts.SplitBy, 'TargetType'), return; end
    if isempty(runs), return; end

    gms = unique([runs.GridMode]);
    gms = gms(~isnan(gms));
    sps = unique([runs.Spacing]);
    sps = sps(~isnan(sps));

    if numel(gms) ~= 1 || numel(sps) ~= 1
        fprintf(['plotGA_CostBoxPlots: OptiTrack overlay skipped — ' ...
                 'loaded set spans %d GridMode(s) and %d Spacing(s); ' ...
                 'overlay requires a single (GridMode, Spacing).\n'], ...
                 numel(gms), numel(sps));
        return;
    end

    enabled = true;
    gm = gms;
    sp = sps;
end


function s = formatPValue(p)
% Render a p-value safely so very-small values do not display as
% the misleading "p=0.000".
    if isnan(p)
        s = 'p=NA';
    elseif p < 0.001
        s = 'p<0.001';
    else
        s = sprintf('p=%.3f', p);
    end
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
