function plotGA_PopulationDiversity(varargin)
% plotGA_PopulationDiversity  Population diversity vs generation.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "The GA actually explores the search space — population diversity
%   declines over generations as selection pressure pulls individuals
%   towards good basins, justifying the choice of GA over a single
%   hill-climber. A run that *does not* lose diversity is one that
%   has not yet converged; a run that loses diversity *too quickly*
%   has prematurely converged."
%
% Why this plot belongs in the methods section
%   plotGA_Convergence shows what the GA produced. THIS plot shows what
%   the GA *did*. Together they let an examiner judge whether the
%   stochastic search actually visited multiple basins (i.e. whether a
%   GA was even the right choice over a deterministic hill-climber).
%
% Diversity metric
%   The current GA save format records AvgCostHistory, ConvergenceHistory
%   (best-cost), and TopTenAvgCostHistory only — no per-locus gene
%   diversity. So this function uses the COST-SPREAD PROXY as its
%   primary metric:
%
%       diversity(g) = (avg_cost(g) - best_cost(g)) / max(avg_cost(g), eps)
%
%   This is bounded in [0, 1+] and decays as the population concentrates
%   around the best individual. It is a proxy for genotype diversity,
%   not a direct measure — but it is monotonically related under the
%   reasonable assumption that lower-cost individuals are genetically
%   close to the elite. The plot is annotated to make the proxy
%   explicit.
%
%   If you later add a precomputed PopulationDiversityHistory or
%   PopulationStdHistory field to saveData, this function will
%   automatically prefer that over the proxy.
%   Lookup priority:
%     1. saveData.PopulationDiversityHistory  (preferred; precomputed)
%     2. saveData.PopulationStdHistory        (per-gen mean per-locus std)
%     3. proxy from AvgCostHistory + ConvergenceHistory  (current path)
%
% Strengths
%   - Median + IQR ribbon across runs (robust, consistent with
%     plotGA_Convergence).
%   - Sample size n is shown in the legend.
%   - Optional overlay of mean best-cost curve from plotGA_Convergence
%     so the reader sees diversity collapse aligned with cost descent.
% =====================================================================
%
%   plotGA_PopulationDiversity('Name', Value, ...)
%
%   Name-Value Parameters (passed through to loadGARuns, plus):
%     'RunDir'         - Directory containing run .mat files (default: pwd)
%     'OverlayCost'    - true to overlay normalised median-best-cost
%                        on a secondary y-axis (default: true)
%     'LogScale'       - true for semilogy on diversity axis (default: false)
%     'SaveAs'         - Output filename (default: auto)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'RunDir',      pwd,   @ischar);
    addParameter(p, 'OverlayCost', true,  @islogical);
    addParameter(p, 'LogScale',    false, @islogical);
    addParameter(p, 'MaxIt',       [],    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SaveAs',      '',    @ischar);
    parse(p, varargin{:});

    opts = p.Results;
    loaderArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
    loaderArgs = loaderArgs(:)';

    [runs, filterDesc] = loadGARuns(loaderArgs{:});
    nRuns = length(runs);
    sty = gaPlotStyle();

    %% Pass 1 — load every usable run, regardless of MaxIt -------------
    runDiv   = cell(nRuns, 1);
    runBest  = cell(nRuns, 1);
    runMaxIt = nan(nRuns, 1);
    proxyTmp = false(nRuns, 1);

    for i = 1:nRuns
        if ~isfield(runs(i), 'RunFilename') || isempty(runs(i).RunFilename)
            continue;
        end
        if isfield(runs(i), 'NumCameras')
            matPath = resolveRunPath(runs(i).RunFilename, runs(i).NumCameras, opts.RunDir);
        else
            matPath = resolveRunPath(runs(i).RunFilename, [], opts.RunDir);
        end
        if ~isfile(matPath), continue; end

        data = load(matPath, 'saveData');
        if ~isfield(data, 'saveData'), continue; end
        sd = data.saveData;

        thisDiv = [];
        if isfield(sd, 'PopulationDiversityHistory') && ~isempty(sd.PopulationDiversityHistory)
            thisDiv = sd.PopulationDiversityHistory(:)';
        elseif isfield(sd, 'PopulationStdHistory') && ~isempty(sd.PopulationStdHistory)
            thisDiv = sd.PopulationStdHistory(:)';
        elseif isfield(sd, 'AvgCostHistory') && isfield(sd, 'ConvergenceHistory') && ...
               numel(sd.AvgCostHistory) == numel(sd.ConvergenceHistory)
            avgC = sd.AvgCostHistory(:)';
            bC   = sd.ConvergenceHistory(:)';
            thisDiv = (avgC - bC) ./ max(avgC, eps);
            proxyTmp(i) = true;
        end
        if isempty(thisDiv), continue; end

        runDiv{i}   = thisDiv;
        runBest{i}  = sd.ConvergenceHistory(:)';
        runMaxIt(i) = numel(thisDiv);
    end

    %% Pass 2 — pick the modal MaxIt and drop the rest -----------------
    haveData = ~cellfun(@isempty, runDiv);
    if ~any(haveData)
        warning(['plotGA_PopulationDiversity: no diversity history found in ' ...
                 'any run (need PopulationDiversityHistory, PopulationStdHistory, ' ...
                 'or AvgCostHistory + ConvergenceHistory). Skipping.']);
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
        fprintf(['plotGA_PopulationDiversity: dropped %d run(s) whose MaxIt ' ...
                 '!= %d (the modal value). MaxIt values present: %s\n'], ...
                 nDropped, targetMaxIt, mat2str(unique(runMaxIt(haveData))'));
    end

    keepIdx = find(keepMask);
    nValid  = numel(keepIdx);
    if nValid == 0
        warning('plotGA_PopulationDiversity: no runs left after MaxIt filter.');
        return;
    end

    divHists  = zeros(nValid, targetMaxIt);
    bestHists = zeros(nValid, targetMaxIt);
    proxyUsed = false(1, nValid);
    for j = 1:nValid
        k = keepIdx(j);
        divHists(j,:)  = runDiv{k};
        bestHists(j,:) = runBest{k};
        proxyUsed(j)   = proxyTmp(k);
    end

    nGen = size(divHists, 2);
    gens = 1:nGen;

    %% Median + IQR
    divMed   = median(divHists, 1);
    divQ25   = quantile(divHists, 0.25, 1);
    divQ75   = quantile(divHists, 0.75, 1);

    %% Plot
    fig = figure('Units', 'inches', ...
        'Position', [1 1 sty.FigWidthFull sty.FigHeight], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);
    ax = axes(fig);
    hold(ax, 'on');

    % Faint individual traces
    for i = 1:nValid
        plot(ax, gens, divHists(i,:), ...
            'Color', [0.5 0.5 0.5 0.25], 'LineWidth', sty.LineWidthThin);
    end

    % IQR ribbon
    fillX = [gens, fliplr(gens)];
    fillY = [divQ75, fliplr(divQ25)];
    fill(ax, fillX, fillY, sty.CostFuncColors(2,:), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Median line
    hMed = plot(ax, gens, divMed, '-', ...
        'Color', sty.CostFuncColors(2,:), 'LineWidth', sty.LineWidth + 0.4);

    legendHandles = hMed;
    if any(proxyUsed)
        divLabel = sprintf('Diversity proxy median (Q25–Q75, n = %d)', nValid);
    else
        divLabel = sprintf('Population diversity median (Q25–Q75, n = %d)', nValid);
    end
    legendLabels = {divLabel};

    % Optional cost overlay on right axis
    if opts.OverlayCost && ~isempty(bestHists)
        bestMed = median(bestHists, 1);
        normCost = bestMed / max(bestMed);
        yyaxis(ax, 'right');
        hCost = plot(ax, gens, normCost, '--', ...
            'Color', sty.CostFuncColors(3,:), 'LineWidth', sty.LineWidth);
        ylabel(ax, 'Best cost (median, normalised)', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylim(ax, [0 1.05]);
        yyaxis(ax, 'left');

        legendHandles(end+1) = hCost;
        legendLabels{end+1}  = 'Best cost median (normalised)';
    end

    if opts.LogScale
        set(ax, 'YScale', 'log');
    end

    hold(ax, 'off');

    %% Formatting
    xlabel(ax, 'Generation', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    ylabel(ax, 'Population diversity (a.u.)', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');

    legend(legendHandles, legendLabels, ...
        'Location', 'northeast', 'FontSize', sty.FontSizeLegend);

    if any(proxyUsed)
        title(ax, ...
          'Diversity computed as cost-spread proxy: (avg − best) / avg', ...
          'FontWeight', 'normal', 'FontSize', sty.FontSizeAxis, ...
          'FontName', sty.FontName);
    end

    applyThesisStyle(fig);

    %% Print summary
    fprintf('\nDiversity summary (%s):\n', filterDesc);
    fprintf('  Runs:               %d\n', nValid);
    fprintf('  Generations:        %d\n', nGen);
    fprintf('  Initial diversity median: %.4g\n', divMed(1));
    fprintf('  Final diversity median:   %.4g\n', divMed(end));
    if any(proxyUsed)
        fprintf('  NOTE: %d of %d runs used the cost-spread proxy.\n', ...
                sum(proxyUsed), nValid);
    end

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('GA_Diversity_%s', filterDesc);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end
