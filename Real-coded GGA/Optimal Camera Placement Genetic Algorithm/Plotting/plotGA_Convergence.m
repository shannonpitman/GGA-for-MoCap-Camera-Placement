function plotGA_Convergence(varargin)
% plotGA_Convergence  Convergence curves with median + IQR ribbon, the
% winning run highlighted, robustness statistic in the legend, and a
% side-by-side Cold | Warm subplot layout (default).
%
% =====================================================================
% EXAMINER REVIEW (questions an external examiner would push on)
% =====================================================================
% What this plot claims to show
%   "Under the filtered conditions, the GA reliably converges to a
%   low-cost solution, with bounded run-to-run variance AND a stated
%   proportion of runs that landed near the per-start-strategy best.
%   The single best run's full convergence path is shown so the reader
%   can judge how quickly that winning trajectory dropped, and the
%   run's chronological index is annotated so the reader sees how many
%   independent restarts were needed before the winner appeared."
%
% Strengths
%   - Side-by-side Cold | Warm split avoids the false-premature-
%     convergence reading that arises when two distinct populations are
%     pooled (their distinct asymptotes look like stragglers).
%   - "Within X% of best" robustness is computed *per panel* so the
%     reader sees how stable each start strategy is on its own terms.
%   - The best-run index (1..N within its strategy) is annotated, so the
%     reader can judge how many GA restarts they would typically need.
%   - Median + IQR ribbon + faint individual traces show both central
%     tendency and spread, robust to skew and zero-bound.
%   - Top-10 average sub-curve gives evidence of population-level
%     improvement, not just elitist drift.
%
% Robustness against log inconsistency
%   - Runs whose MaxIt differs from the modal MaxIt of the filtered
%     batch are dropped with a console warning. This handles the case
%     where a stray pilot run is mixed in with the production batch.
% =====================================================================
%
%   plotGA_Convergence('Name', Value, ...)
%
%   Loads individual run .mat files (via RunFilename in the log) to
%   extract generation-by-generation cost histories. Default layout is
%   a 1x2 tiled figure: left = Cold-start runs, right = Warm-start
%   runs. Pass 'WarmStart', true/false to force a single panel, or set
%   'SplitWarmCold', 'off'.
%
%   Per panel, plots:
%     - Faint individual best-cost traces (one per run)
%     - Median best-cost curve + shaded Q25-Q75 IQR ribbon (default)
%       OR mean +/- 1 std envelope (CentralStat='mean')
%     - Median (or mean) top-10-average curve across runs
%     - Full trace of the winning run highlighted in deep red, with
%       an inline callout next to the line's end stating final cost
%       and the winning run's chronological index (e.g. "Run 7/10").
%     - Horizontal dashed line at the panel's best cost
%     - Robustness annotation in the legend: "X/N runs within Y% of best"
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'RunDir'         - Directory containing run .mat files (default: pwd)
%     'CentralStat'    - 'median' (default) or 'mean'
%     'Tolerance'      - Fractional tolerance for the "within X%"
%                        robustness count (default: 0.05 → 5%).
%     'ShowAvgCost'    - true to plot population-average curve (default: false)
%     'ShowTopTen'     - true to plot top-10-average curve (default: true)
%     'ShowIndividual' - true to show individual run traces (default: true)
%     'ShowBestLine'   - true to show overall best cost yline (default: true)
%     'ShowBestRun'    - true to overlay the winning run's full trace
%                        (default: true)
%     'AnnotateInit'   - true to annotate initial-generation cost
%                        (default: true)
%     'LogScale'       - true for semilogy (default: false)
%     'MaxIt'          - scalar MaxIt to filter to; if empty (default),
%                        the modal MaxIt across the filtered batch is
%                        used and runs at any other MaxIt are dropped
%                        with a warning.
%     'SplitWarmCold'  - 'auto' (default), 'on', or 'off'.
%                        'auto': split into Cold | Warm panels unless
%                        WarmStart was passed in the caller args.
%     'SaveAs'         - Output filename without extension (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;  % pass unmatched to loadGARuns
    addParameter(p, 'RunDir',         pwd,      @ischar);
    addParameter(p, 'CentralStat',    'median', @(s) ischar(s) || isstring(s));
    addParameter(p, 'Tolerance',      0.05,     @(x) isnumeric(x) && isscalar(x) && x>=0);
    addParameter(p, 'ShowAvgCost',    false,    @islogical);
    addParameter(p, 'ShowTopTen',     true,     @islogical);
    addParameter(p, 'ShowIndividual', true,     @islogical);
    addParameter(p, 'ShowBestLine',   true,     @islogical);
    addParameter(p, 'ShowBestRun',    true,     @islogical);
    addParameter(p, 'AnnotateInit',   true,     @islogical);
    addParameter(p, 'LogScale',       false,    @islogical);
    addParameter(p, 'MaxIt',          [],       @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SplitWarmCold',  'auto',   @(s) ischar(s) || isstring(s));
    addParameter(p, 'SaveAs',         '',       @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    centralStat = lower(string(opts.CentralStat));
    if ~ismember(centralStat, ["median", "mean"])
        error('plotGA_Convergence:badCentralStat', ...
            'CentralStat must be ''median'' or ''mean''.');
    end

    % Detect whether the caller has constrained WarmStart explicitly.
    explicitWarm = isfield(p.Unmatched, 'WarmStart');

    splitMode = lower(string(opts.SplitWarmCold));
    switch splitMode
        case "auto"
            doSplit = ~explicitWarm;
        case "on"
            doSplit = true;
        case "off"
            doSplit = false;
        otherwise
            error('plotGA_Convergence:badSplit', ...
                'SplitWarmCold must be ''auto'', ''on'', or ''off''.');
    end

    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    sty = gaPlotStyle();

    %% Load the combined run set
    [runs, filterDesc] = loadGARuns(loaderArgs{:});

    %% Optional: split by warm/cold
    if doSplit
        warmFlags = logical([runs.WarmStart]);
        coldRuns = runs(~warmFlags);
        warmRuns = runs( warmFlags);

        if isempty(coldRuns) || isempty(warmRuns)
            fprintf(['plotGA_Convergence: only %d cold / %d warm runs ' ...
                     '— falling back to single panel.\n'], ...
                     length(coldRuns), length(warmRuns));
            doSplit = false;
        end
    end

    %% =========================================================
    %  Render
    %  =========================================================
    if doSplit
        % --- 1x2 figure: Cold | Warm -------------------------------
        fig = figure('Units', 'inches', ...
            'Position', [1 1 sty.FigWidthDouble sty.FigHeightWide], ...
            'PaperPositionMode', 'auto', ...
            'Color', sty.BackgroundColor);
        tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', ...
                                    'TileSpacing', 'compact');

        ax1 = nexttile(tl, 1);
        meta1 = renderConvergencePanel(ax1, coldRuns, 'Cold-start', opts, sty);

        ax2 = nexttile(tl, 2);
        meta2 = renderConvergencePanel(ax2, warmRuns, 'Warm-start', opts, sty);

        % Link y-axes so the visual comparison is direct
        if ~isempty(meta1) && ~isempty(meta2)
            yloAll = min(meta1.ylo, meta2.ylo);
            yhiAll = max(meta1.yhi, meta2.yhi);
            ylim(ax1, [yloAll, yhiAll]);
            ylim(ax2, [yloAll, yhiAll]);
        end

        % Hide y-axis label on right panel — they share the scale
        if ~isempty(meta2)
            ylabel(ax2, '');
        end

        % Supertitle describing the filter condition
        title(tl, prettifyFilterDesc(filterDesc), ...
            'FontSize', sty.FontSizeTitle + 1, ...
            'FontName', sty.FontName, ...
            'FontWeight', 'bold');

        %% Summary print
        fprintf('\nConvergence summary (%s):\n', filterDesc);
        printPanelSummary('Cold-start', meta1, opts);
        printPanelSummary('Warm-start', meta2, opts);

    else
        % --- Single panel -------------------------------------------
        fig = figure('Units', 'inches', ...
            'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
            'PaperPositionMode', 'auto', ...
            'Color', sty.BackgroundColor);
        ax = axes(fig);
        if explicitWarm
            if p.Unmatched.WarmStart
                titleStr = 'Warm-start';
            else
                titleStr = 'Cold-start';
            end
        else
            titleStr = '';   % no split, no title
        end
        meta = renderConvergencePanel(ax, runs, titleStr, opts, sty);
        if ~isempty(meta)
            ylim(ax, [meta.ylo, meta.yhi]);
        end
        fprintf('\nConvergence summary (%s):\n', filterDesc);
        printPanelSummary(titleStr, meta, opts);
    end

    %% Lock thesis style (white bg, black text) and export
    applyThesisStyle(fig);

    if isempty(opts.SaveAs)
        outName = sprintf('GA_Convergence_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


% =========================================================================
% renderConvergencePanel — Plot a single panel on the supplied axes.
% =========================================================================
function meta = renderConvergencePanel(ax, runs, titleStr, opts, sty)
% Returns meta with fields: nValid, ylo, yhi, finalCentral, finalBest,
% nWithin, pctWithin, bestRunIdx, droppedMaxIt.
%
% Loads per-run convergence histories, filters out runs whose MaxIt
% deviates from the modal MaxIt (with a warning), and plots all curves
% on the given axes.

    meta = [];
    nIn  = length(runs);

    %% Discover modal MaxIt and load histories
    matPaths  = strings(nIn, 1);
    runMaxIt  = nan(nIn, 1);
    runHists  = cell(nIn, 1);
    runAvgs   = cell(nIn, 1);
    runTops   = cell(nIn, 1);
    runIdxs   = nan(nIn, 1);   % chronological index in the *input* set

    for i = 1:nIn
        if ~isfield(runs(i), 'RunFilename') || isempty(runs(i).RunFilename)
            warning('Run %d has no RunFilename — skipping.', i);
            continue;
        end
        if isfield(runs(i), 'NumCameras')
            matPath = resolveRunPath(runs(i).RunFilename, runs(i).NumCameras, opts.RunDir);
        else
            matPath = resolveRunPath(runs(i).RunFilename, [], opts.RunDir);
        end
        if ~isfile(matPath)
            warning('File not found: %s — skipping.', matPath);
            continue;
        end
        d = load(matPath, 'saveData');
        if ~isfield(d, 'saveData')
            continue;
        end
        sd = d.saveData;
        runHists{i} = sd.ConvergenceHistory(:)';
        if isfield(sd, 'AvgCostHistory'),       runAvgs{i} = sd.AvgCostHistory(:)';       end
        if isfield(sd, 'TopTenAvgCostHistory'), runTops{i} = sd.TopTenAvgCostHistory(:)'; end
        runMaxIt(i) = length(runHists{i});
        runIdxs(i)  = i;
        matPaths(i) = matPath;
    end

    %% Pick MaxIt (user-supplied wins; else modal across batch)
    haveData = ~cellfun(@isempty, runHists);
    if ~any(haveData)
        warning('No usable runs in this panel.');
        return;
    end

    if ~isempty(opts.MaxIt)
        targetMaxIt = opts.MaxIt;
    else
        targetMaxIt = mode(runMaxIt(haveData));
    end

    keepMask = haveData & (runMaxIt == targetMaxIt);
    nDropped = sum(haveData & ~keepMask);
    if nDropped > 0
        fprintf(['  [%s] Dropped %d run(s) whose MaxIt != %d (the modal ' ...
                 'value). MaxIt values present: %s\n'], ...
                 titleStr, nDropped, targetMaxIt, ...
                 mat2str(unique(runMaxIt(haveData))'));
    end

    keepIdx = find(keepMask);
    nValid  = numel(keepIdx);
    if nValid == 0
        warning('Panel %s has 0 valid runs after MaxIt filter.', titleStr);
        return;
    end

    nGen = targetMaxIt;
    gens = 1:nGen;
    bestHists   = zeros(nValid, nGen);
    avgHists    = [];
    topTenHists = [];
    for j = 1:nValid
        k = keepIdx(j);
        bestHists(j, :) = runHists{k};
        if ~isempty(runAvgs{k})
            avgHists(end+1, :) = runAvgs{k};         %#ok<AGROW>
        end
        if ~isempty(runTops{k})
            topTenHists(end+1, :) = runTops{k};     %#ok<AGROW>
        end
    end
    finalBests = bestHists(:, end);

    %% Stats
    centralStat = lower(string(opts.CentralStat));
    if centralStat == "median"
        bestCentral = median(bestHists, 1);
        bestUpper   = quantile(bestHists, 0.75, 1);
        bestLower   = quantile(bestHists, 0.25, 1);
        envLabel    = 'Q25–Q75 IQR';
        centralLbl  = 'median';
    else
        bestCentral = mean(bestHists, 1);
        bestStd     = std(bestHists, 0, 1);
        bestUpper   = bestCentral + bestStd;
        bestLower   = max(bestCentral - bestStd, 0);
        envLabel    = '\pm1 std';
        centralLbl  = 'mean';
    end

    overallBest    = min(finalBests);
    overallBestRow = find(finalBests == overallBest, 1);
    overallBestRunIdx = keepIdx(overallBestRow);

    tol = opts.Tolerance;
    if overallBest > 0
        thresh = overallBest * (1 + tol);
    else
        thresh = overallBest + tol;
    end
    nWithin   = sum(finalBests <= thresh);
    pctWithin = 100 * nWithin / nValid;

    %% Plot
    hold(ax, 'on');
    plotFn = @plot;
    if opts.LogScale
        plotFn = @semilogy;
    end

    % Individual traces
    if opts.ShowIndividual
        for i = 1:nValid
            plotFn(ax, gens, bestHists(i,:), ...
                'Color', [0.5 0.5 0.5 0.25], ...
                'LineWidth', sty.LineWidthThin, ...
                'HandleVisibility', 'off');
        end
    end

    % Envelope ribbon
    fillX = [gens, fliplr(gens)];
    fillY = [bestUpper, fliplr(bestLower)];
    if opts.LogScale
        fillY = max(fillY, 1e-12);
    end
    fill(ax, fillX, fillY, sty.CostFuncColors(3,:), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % Central tendency line
    hCentral = plotFn(ax, gens, bestCentral, '-', ...
        'Color', sty.CostFuncColors(3,:), 'LineWidth', sty.LineWidth + 0.4);
    legendHandles = hCentral;
    legendLabels  = {sprintf('Best cost %s (%s, n = %d)', ...
                              centralLbl, envLabel, nValid)};

    % Envelope boundary lines
    plotFn(ax, gens, bestUpper, '-', ...
        'Color', [sty.CostFuncColors(3,:) 0.4], ...
        'LineWidth', sty.LineWidthThin, 'HandleVisibility', 'off');
    plotFn(ax, gens, bestLower, '-', ...
        'Color', [sty.CostFuncColors(3,:) 0.4], ...
        'LineWidth', sty.LineWidthThin, 'HandleVisibility', 'off');

    % Top-10 avg
    if opts.ShowTopTen && ~isempty(topTenHists)
        if centralStat == "median"
            topTenCentral = median(topTenHists, 1);
        else
            topTenCentral = mean(topTenHists, 1);
        end
        hTopTen = plotFn(ax, gens, topTenCentral, '--', ...
            'Color', sty.CostFuncColors(1,:), 'LineWidth', sty.LineWidth);
        legendHandles(end+1) = hTopTen;
        legendLabels{end+1}  = sprintf('Top-10 avg %s (n = %d)', ...
                                        centralLbl, size(topTenHists,1));
    end

    % Population average (optional)
    if opts.ShowAvgCost && ~isempty(avgHists)
        if centralStat == "median"
            avgCentral = median(avgHists, 1);
        else
            avgCentral = mean(avgHists, 1);
        end
        hAvg = plotFn(ax, gens, avgCentral, ':', ...
            'Color', sty.CostFuncColors(2,:), 'LineWidth', sty.LineWidth);
        legendHandles(end+1) = hAvg;
        legendLabels{end+1}  = sprintf('Population avg %s (n = %d)', ...
                                        centralLbl, size(avgHists,1));
    end

    % --- Panel-strategy label used in legend captions ---------------
    % Anchor every "best" / "robustness" string to the strategy this
    % panel is showing, so the reader can see at a glance that the
    % 5% reference is the cold-start best vs the warm-start best.
    if ~isempty(titleStr)
        strategyLbl = lower(titleStr);   % e.g. 'cold-start'
        bestRefLbl  = sprintf('%s best', strategyLbl);
    else
        strategyLbl = '';
        bestRefLbl  = 'overall best';
    end

    % Overall best yline
    if opts.ShowBestLine
        hBest = yline(ax, overallBest, '-.', ...
            'Color', [0.7 0.15 0.15], 'LineWidth', sty.LineWidth, ...
            'Alpha', 0.8);
        legendHandles(end+1) = hBest;
        legendLabels{end+1}  = sprintf('%s = %.4f', ...
                                        capitalise(bestRefLbl), overallBest);
    end

    % Best-run trace (red)
    if opts.ShowBestRun
        bestRunCol = [0.70 0.15 0.15];
        bestTrace  = bestHists(overallBestRow, :);
        hBestRun = plotFn(ax, gens, bestTrace, '-', ...
            'Color', bestRunCol, 'LineWidth', sty.LineWidth + 0.6);
        legendHandles(end+1) = hBestRun;
        legendLabels{end+1}  = sprintf('Best run trace: run %d/%d (final = %.4f)', ...
                                        overallBestRunIdx, nValid, bestTrace(end));

        % Endpoint marker
        plot(ax, nGen, bestTrace(end), 'o', ...
            'MarkerFaceColor', bestRunCol, 'MarkerEdgeColor', 'k', ...
            'MarkerSize', sty.MarkerSizeLg, 'HandleVisibility', 'off');

        % Inline callout *next to the line end*. Extend the x-axis so
        % the box sits to the right of the endpoint with a short
        % leader, rather than overlapping the trace.
        calloutText = sprintf(' Run %d/%d\n Final = %.4f', ...
                              overallBestRunIdx, nValid, bestTrace(end));
        % Place to the right of the endpoint; alignment LEFT
        xCallout = nGen + 0.015 * nGen;
        yCallout = bestTrace(end);
        text(ax, xCallout, yCallout, calloutText, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'middle', ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'BackgroundColor', 'w', 'EdgeColor', bestRunCol, ...
            'Color', 'k', 'Margin', 4);
    end

    % Robustness as a no-line legend entry — name the strategy so
    % the comparison reference is unambiguous (cold-vs-cold, warm-vs-warm).
    hRob = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none');
    legendHandles(end+1) = hRob;
    legendLabels{end+1}  = sprintf('Within %.0f%% of %s: %d/%d (%.0f%%)', ...
                                    100*tol, bestRefLbl, ...
                                    nWithin, nValid, pctWithin);

    % Initial-population marker
    if opts.AnnotateInit
        if centralStat == "median"
            initCentral = median(bestHists(:,1));
        else
            initCentral = mean(bestHists(:,1));
        end
        plot(ax, 1, initCentral, 'kv', ...
            'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', 'k', ...
            'MarkerSize', sty.MarkerSize, 'HandleVisibility', 'off');
        text(ax, 1 + 0.02*nGen, initCentral, ...
            sprintf('init %s = %.3f', centralLbl, initCentral), ...
            'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
            'VerticalAlignment', 'bottom', 'Color', 'k');
    end

    if opts.LogScale
        set(ax, 'YScale', 'log');
    end

    hold(ax, 'off');

    %% Formatting
    xlabel(ax, 'Generation', 'FontSize', sty.FontSizeAxis, ...
                              'FontName', sty.FontName);
    ylabel(ax, 'Cost', 'FontSize', sty.FontSizeAxis, ...
                       'FontName', sty.FontName);
    if ~isempty(titleStr)
        title(ax, titleStr, 'FontSize', sty.FontSizeTitle, ...
                            'FontName', sty.FontName, ...
                            'FontWeight', 'normal');
    end
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
            'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');

    yLo = max(0, min(bestLower) * 0.9);
    yHi = max(bestUpper) * 1.05;
    if opts.LogScale
        yLo = max(1e-6, yLo);
    end
    if ~isfinite(yLo) || ~isfinite(yHi) || yHi <= yLo
        yLo = 0; yHi = 1;
    end

    % Extend x-axis to give the callout room to breathe
    xlim(ax, [1, nGen * 1.18]);

    legend(ax, legendHandles, legendLabels, ...
        'Location', 'northeast', 'FontSize', sty.FontSizeLegend);

    %% Pack metadata
    meta.nValid       = nValid;
    meta.ylo          = yLo;
    meta.yhi          = yHi;
    meta.finalCentral = bestCentral(end);
    meta.finalIQRLo   = bestLower(end);
    meta.finalIQRHi   = bestUpper(end);
    meta.finalBest    = overallBest;
    meta.bestRunIdx   = overallBestRunIdx;
    meta.nWithin      = nWithin;
    meta.pctWithin    = pctWithin;
    meta.nDropped     = nDropped;
    meta.maxIt        = targetMaxIt;
end


% =========================================================================
function s = capitalise(s)
% Lightweight title-case for legend strings.
    if isempty(s), return; end
    s(1) = upper(s(1));
end


% =========================================================================
function pretty = prettifyFilterDesc(desc)
% Convert e.g. "7C_CF3_TT1_GM1_sp100cm" into a human-readable supertitle.
    tokens = strsplit(desc, '_');
    parts  = strings(1, 0);
    cfNames = {'Resolution Uncertainty', 'Dynamic Occlusion', 'Combined'};
    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform grid', 'Normal grid'};
    for k = 1:numel(tokens)
        t = tokens{k};
        if endsWith(t, 'C') && all(isstrprop(t(1:end-1), 'digit'))
            n = str2double(t(1:end-1));
            parts(end+1) = sprintf('%d cameras', n);                         %#ok<AGROW>
        elseif startsWith(t, 'CF')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(cfNames)
                parts(end+1) = sprintf('CF%d (%s)', n, cfNames{n});          %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'TT')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(ttNames)
                parts(end+1) = ttNames{n};                                  %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'GM')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(gmNames)
                parts(end+1) = gmNames{n};                                  %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'sp') && endsWith(t, 'cm')
            num = str2double(t(3:end-2));
            parts(end+1) = sprintf('%.2f m spacing', num/100);              %#ok<AGROW>
        else
            parts(end+1) = string(t);                                       %#ok<AGROW>
        end
    end
    pretty = char(strjoin(parts, '  ·  '));
end


% =========================================================================
function printPanelSummary(label, meta, opts)
    if isempty(meta)
        fprintf('  [%s] no data\n', label);
        return;
    end
    fprintf('  [%s]\n', label);
    fprintf('    Runs:                   %d (dropped %d for MaxIt mismatch)\n', ...
                                          meta.nValid, meta.nDropped);
    fprintf('    Generations:            %d\n', meta.maxIt);
    fprintf('    Final central cost:     %.6f\n', meta.finalCentral);
    fprintf('    Final envelope:         [%.6f, %.6f]\n', ...
                                          meta.finalIQRLo, meta.finalIQRHi);
    fprintf('    Panel best:             %.6f (run %d/%d)\n', ...
                                          meta.finalBest, ...
                                          meta.bestRunIdx, meta.nValid);
    fprintf('    Robustness:             %d/%d (%.1f%%) within %.0f%% of best\n', ...
                                          meta.nWithin, meta.nValid, ...
                                          meta.pctWithin, 100*opts.Tolerance);
end
