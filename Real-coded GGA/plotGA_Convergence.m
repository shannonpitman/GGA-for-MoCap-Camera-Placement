function plotGA_Convergence(varargin)
% plotGA_Convergence  Convergence curves with mean +/- std envelope.
%   plotGA_Convergence('Name', Value, ...)
%
%   Loads individual run .mat files (via RunFilename in the log) to extract
%   generation-by-generation cost histories. Plots:
%     - Faint individual best-cost traces (one per run)
%     - Mean best-cost curve with shaded +/- 1 std envelope
%     - Mean top-10-average curve (averaged across runs)
%     - Horizontal dashed line at the overall best cost found
%
%   REQUIRED: all filtered runs must share the same MaxIt so the histories
%   align. Use loadGARuns filters to ensure a homogeneous set.
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'RunDir'       - Directory containing individual run .mat files
%                      (default: pwd)
%     'ShowAvgCost'  - true to also plot population-average envelope
%                      (default: false)
%     'ShowTopTen'   - true to also plot top-10-average curve
%                      (default: true)
%     'ShowIndividual' - true to show individual run traces (default: true)
%     'ShowBestLine' - true to show overall best cost line (default: true)
%     'LogScale'     - true for semilogy (default: false)
%     'SaveAs'       - Output filename without extension (default: auto)
%     (all loadGARuns parameters are also accepted)

    p = inputParser;
    p.KeepUnmatched = true;  % pass unmatched to loadGARuns
    addParameter(p, 'RunDir',         pwd,   @ischar);
    addParameter(p, 'ShowAvgCost',    false, @islogical);
    addParameter(p, 'ShowTopTen',     true,  @islogical);
    addParameter(p, 'ShowIndividual', true,  @islogical);
    addParameter(p, 'ShowBestLine',   true,  @islogical);
    addParameter(p, 'LogScale',       false, @islogical);
    addParameter(p, 'SaveAs',         '',    @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    %% Load filtered runs
    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    nRuns = length(runs);
    sty = gaPlotStyle();

    %% Load convergence histories from individual .mat files
    bestHists   = [];
    avgHists    = [];
    topTenHists = [];
    finalBests  = [];

    for i = 1:nRuns
        if ~isfield(runs(i), 'RunFilename') || isempty(runs(i).RunFilename)
            warning('Run %d has no RunFilename — skipping.', i);
            continue;
        end

        matPath = fullfile(opts.RunDir, runs(i).RunFilename);
        if ~isfile(matPath)
            warning('File not found: %s — skipping.', matPath);
            continue;
        end

        data = load(matPath, 'saveData');
        if ~isfield(data, 'saveData')
            warning('No saveData in %s — skipping.', matPath);
            continue;
        end

        sd = data.saveData;

        % Convergence histories (row vectors, one per generation)
        bestHists(end+1, :) = sd.ConvergenceHistory(:)';       %#ok<AGROW>
        finalBests(end+1)   = sd.ConvergenceHistory(end);       %#ok<AGROW>

        if isfield(sd, 'AvgCostHistory')
            avgHists(end+1, :) = sd.AvgCostHistory(:)';         %#ok<AGROW>
        end
        if isfield(sd, 'TopTenAvgCostHistory')
            topTenHists(end+1, :) = sd.TopTenAvgCostHistory(:)'; %#ok<AGROW>
        end
    end

    nValid = size(bestHists, 1);
    if nValid == 0
        error('plotGA_Convergence:noData', ...
            'No valid convergence histories found. Check RunDir and RunFilename.');
    end
    fprintf('Loaded convergence data from %d / %d runs.\n', nValid, nRuns);

    nGen = size(bestHists, 2);
    gens = 1:nGen;

    %% Compute statistics
    bestMean  = mean(bestHists, 1);
    bestStd   = std(bestHists, 0, 1);
    bestUpper = bestMean + bestStd;
    bestLower = max(bestMean - bestStd, 0);  % clamp above zero

    overallBest    = min(finalBests);
    overallBestIdx = find(finalBests == overallBest, 1);

    %% Create figure
    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', 'Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    plotFn = @plot;  % default linear scale
    if opts.LogScale
        plotFn = @semilogy;
    end

    % --- Individual run traces (faint) ---
    if opts.ShowIndividual
        for i = 1:nValid
            plotFn(ax, gens, bestHists(i,:), ...
                'Color', [0.5 0.5 0.5 0.25], ...   % light grey, low alpha
                'LineWidth', sty.LineWidthThin);
        end
    end

    % --- Std envelope (shaded) ---
    fillX = [gens, fliplr(gens)];
    fillY = [bestUpper, fliplr(bestLower)];
    if opts.LogScale
        fillY = max(fillY, 1e-12);
    end
    hEnv = fill(ax, fillX, fillY, sty.CostFuncColors(3,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % --- Mean best-cost line ---
    hMean = plotFn(ax, gens, bestMean, '-', ...
        'Color', sty.CostFuncColors(3,:), 'LineWidth', sty.LineWidth + 0.4);

    legendHandles = [hMean];
    legendLabels  = {sprintf('Best cost mean (n = %d)', nValid)};

    % --- Envelope boundary lines (thin, same colour) ---
    plotFn(ax, gens, bestUpper, '-', ...
        'Color', [sty.CostFuncColors(3,:) 0.4], ...
        'LineWidth', sty.LineWidthThin, 'HandleVisibility', 'off');
    plotFn(ax, gens, bestLower, '-', ...
        'Color', [sty.CostFuncColors(3,:) 0.4], ...
        'LineWidth', sty.LineWidthThin, 'HandleVisibility', 'off');

    % --- Top-10 average (mean across all runs' top-10 averages) ---
    if opts.ShowTopTen && ~isempty(topTenHists)
        topTenMean = mean(topTenHists, 1);
        hTopTen = plotFn(ax, gens, topTenMean, '--', ...
            'Color', sty.CostFuncColors(1,:), 'LineWidth', sty.LineWidth);
        legendHandles(end+1) = hTopTen;
        legendLabels{end+1}  = sprintf('Top-10 average mean (n = %d)', size(topTenHists,1));
    end

    % --- Population average (optional) ---
    if opts.ShowAvgCost && ~isempty(avgHists)
        avgMean = mean(avgHists, 1);
        hAvg = plotFn(ax, gens, avgMean, ':', ...
            'Color', sty.CostFuncColors(2,:), 'LineWidth', sty.LineWidth);
        legendHandles(end+1) = hAvg;
        legendLabels{end+1}  = sprintf('Population average mean (n = %d)', size(avgHists,1));
    end

    % --- Overall best cost line ---
    if opts.ShowBestLine
        hBest = yline(ax, overallBest, '-.', ...
            'Color', [0.7 0.15 0.15], 'LineWidth', sty.LineWidth, ...
            'Alpha', 0.8);
        legendHandles(end+1) = hBest;
        legendLabels{end+1}  = sprintf('Overall best = %.4f', overallBest);
    end

    if opts.LogScale
        set(ax, 'YScale', 'log');
    end

    hold(ax, 'off');

    %% Formatting
    xlabel(ax, 'Generation', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax, 'Cost', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.15);

    % Set y-axis limits with some breathing room
    yLo = max(0, min(bestLower) * 0.9);
    yHi = max(bestUpper) * 1.05;
    if opts.LogScale
        yLo = max(1e-6, yLo);
    end
    ylim(ax, [yLo, yHi]);

    legend(legendHandles, legendLabels, ...
        'Location', 'northeast', 'FontSize', sty.FontSizeLegend);

    %% Print summary
    fprintf('\nConvergence summary (%s):\n', filterDesc);
    fprintf('  Runs:           %d\n', nValid);
    fprintf('  Generations:    %d\n', nGen);
    fprintf('  Mean final cost:  %.6f +/- %.6f\n', bestMean(end), bestStd(end));
    fprintf('  Overall best:     %.6f (run %d)\n', overallBest, overallBestIdx);
    if opts.ShowTopTen && ~isempty(topTenHists)
        fprintf('  Top-10 avg final: %.6f\n', topTenMean(end));
    end

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_Convergence_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', outName);
end