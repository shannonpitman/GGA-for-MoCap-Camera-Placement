function results = batchRun7Cameras(varargin)
% BATCHRUN7CAMERAS  Focused 7-camera rerun on the corrected cost code.
%
% Runs the Tier-0 sanity check + Tier-1 thesis matrix from the post-bugfix
% rerun plan, restricted to 7 cameras / UAV / Uniform grid:
%
%   for CF in {1, 2, 3}
%       run NumColdRepeats cold-start GAs
%       run NumWarmRepeats warm-start GAs, each seeded from the BEST
%       chromosome found across every completed run for this CF so far
%       (cold + previously-finished warm). If a warm run produces a
%       worse cost than the running best, the next warm run still seeds
%       from the running best — never from the most-recent run.
%
% Defaults: 5 cold + 5 warm per cost function = 30 runs total. Each
% finished run is written through saveResults so it appears in the master
% runLog with the corrected cost values.
%
% USAGE
%   batchRun7Cameras()                              % full Tier-1 (CF1,CF2,CF3)
%   batchRun7Cameras('CostFunctions', 3)            % CF3 only
%   batchRun7Cameras('NumColdRepeats', 3, ...
%                    'NumWarmRepeats', 3)
%   batchRun7Cameras('Tier0Only', true)             % just one quick sanity run
%   batchRun7Cameras('DryRun', true)                % print plan, no GA execution
%
% Returns a struct array with one entry per completed run.

    addProjectPaths();

    %% Inputs
    p = inputParser;
    addParameter(p, 'CostFunctions',   [1 2 3], @isnumeric);
    addParameter(p, 'NumColdRepeats',  5,       @isnumeric);
    addParameter(p, 'NumWarmRepeats',  5,       @isnumeric);
    addParameter(p, 'Spacing',         1.0,     @isnumeric);
    addParameter(p, 'TargetType',      1,       @isnumeric); % 1=UAV
    addParameter(p, 'GridMode',        1,       @isnumeric); % 1=Uniform
    addParameter(p, 'MaxGenerations',  100,     @isnumeric);
    % PopulationSize [] => auto-scale as numCams * params-per-camera * 10
    % (= problem.nVar * 10). Pass a numeric value to force a fixed population.
    addParameter(p, 'PopulationSize',  [],      @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'Volume',          [-4 4; -4 4; 0 4], @isnumeric);
    addParameter(p, 'CamLowerBounds',  [-5 -4.5 0  -pi -pi/2 -pi], @isnumeric);
    addParameter(p, 'CamUpperBounds',  [ 5  4.5 4.8 pi  pi/2  pi], @isnumeric);
    addParameter(p, 'Tier0Only',       false, @islogical);
    addParameter(p, 'DryRun',          false, @islogical);
    addParameter(p, 'SuppressPlots',   true,  @islogical);
    parse(p, varargin{:});
    cfg = p.Results;

    numCams = 7;

    % Population scales with chromosome length unless a fixed PopulationSize
    % was supplied: nPop = numCams * numParams * 10 (7 cams -> 420).
    numParams = numel(cfg.CamLowerBounds);
    if isempty(cfg.PopulationSize)
        popSize = numCams * numParams * 10;
    else
        popSize = cfg.PopulationSize;
    end

    if cfg.Tier0Only
        cfg.CostFunctions  = 3;
        cfg.NumColdRepeats = 1;
        cfg.NumWarmRepeats = 0;
        cfg.MaxGenerations = max(round(cfg.MaxGenerations / 2), 50);
        fprintf('[Tier 0] sanity run: 7C / CF3 / %d generations / single seed\n', ...
            cfg.MaxGenerations);
    end

    %% Pre-flight summary
    nTotal = numel(cfg.CostFunctions) * (cfg.NumColdRepeats + cfg.NumWarmRepeats);
    fprintf('\n  7-CAMERA BATCH (%d total runs)\n', nTotal);
    fprintf('  Cost functions:   %s\n', mat2str(cfg.CostFunctions));
    fprintf('  Cold per CF:      %d\n', cfg.NumColdRepeats);
    fprintf('  Warm per CF:      %d (each seeded from running best, cold+warm)\n', ...
        cfg.NumWarmRepeats);
    fprintf('  Target / grid:    UAV / Uniform / spacing %.2fm\n', cfg.Spacing);
    fprintf('  GA budget:        %d generations x population %d%s\n\n', ...
        cfg.MaxGenerations, popSize, ...
        ternary(isempty(cfg.PopulationSize), ' (auto: 7 x 6 x 10)', ' (fixed)'));
    if cfg.DryRun
        fprintf('  DRY RUN — no GA executed.\n\n');
        results = struct([]);
        return;
    end

    %% Set up the (mostly fixed) per-run inputs once
    volume   = cfg.Volume;
    if cfg.TargetType == 2
        volume(3,:) = [0 0.5];
    end
    spacing  = cfg.Spacing;

    %% Suppress figures inside RunGA / visualizeCameraCoverage
    if cfg.SuppressPlots
        prevVis = get(0, 'DefaultFigureVisible');
        set(0, 'DefaultFigureVisible', 'off');
        cleanup = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));
    end

    results   = struct([]);
    runIdxAll = 0;
    batchTic  = tic;

    %% Outer loop: cost functions
    for cfType = cfg.CostFunctions(:)'
        fprintf('=== CF%d : cold-start phase ===\n', cfType);

        % Track the running best chromosome and cost for THIS CF.
        % This is what every warm run will be seeded from, so a warm run
        % that lands a worse cost never poisons the next warm run.
        bestSoFarCost  = inf;
        bestSoFarChrom = [];
        bestSoFarLabel = '';

        for r = 1:cfg.NumColdRepeats
            runIdxAll = runIdxAll + 1;
            fprintf('  cold rep %d/%d  [%d/%d total]\n', ...
                r, cfg.NumColdRepeats, runIdxAll, nTotal);

            specs = makeSpecs(numCams, cfType, cfg, volume, spacing);
            specs.warmStart       = false;
            specs.warmChromosomes = [];

            problem = setupProblem(numCams, cfType, ...
                cfg.CamUpperBounds, cfg.CamLowerBounds);
            params  = setupGAparams(cfg.MaxGenerations, popSize);

            tic;
            out = RunGA(problem, params, specs);
            elapsed = toc;

            covStats = visualizeCameraCoverage(out, specs);
            saveResults(out, specs, params, elapsed, cfType, false, covStats);

            if out.bestsol.Cost < bestSoFarCost
                bestSoFarCost  = out.bestsol.Cost;
                bestSoFarChrom = out.bestsol.Chromosome;
                bestSoFarLabel = sprintf('cold rep %d', r);
            end

            results(end+1).CF          = cfType;          %#ok<AGROW>
            results(end).WarmStart     = false;
            results(end).Repeat        = r;
            results(end).BestCost      = out.bestsol.Cost;
            results(end).Elapsed       = elapsed;
            results(end).MinCoverage   = covStats.minCoverage;
            results(end).ZeroPct       = covStats.zeroCamerasPercent;

            fprintf('    cost = %.6f | min cov = %d | 0-cov = %.1f%% | %.1fs\n', ...
                out.bestsol.Cost, covStats.minCoverage, ...
                covStats.zeroCamerasPercent, elapsed);
            close all;
        end

        if cfg.NumWarmRepeats > 0
            fprintf('=== CF%d : warm-start phase (initial seed = %s, cost = %.6f) ===\n', ...
                cfType, bestSoFarLabel, bestSoFarCost);

            for r = 1:cfg.NumWarmRepeats
                runIdxAll = runIdxAll + 1;

                % Always seed from the running best; print where it came from.
                seed     = bestSoFarChrom;
                seedFrom = bestSoFarLabel;
                seedCost = bestSoFarCost;
                fprintf('  warm rep %d/%d  [%d/%d total]   seed: %s (%.6f)\n', ...
                    r, cfg.NumWarmRepeats, runIdxAll, nTotal, seedFrom, seedCost);

                specs = makeSpecs(numCams, cfType, cfg, volume, spacing);

                problem  = setupProblem(numCams, cfType, ...
                    cfg.CamUpperBounds, cfg.CamLowerBounds);
                perturb  = Mutate(seed, 1, 0.5);
                perturb  = max(perturb, problem.VarMin);
                perturb  = min(perturb, problem.VarMax);
                specs.warmStart       = true;
                specs.warmChromosomes = [seed; perturb];

                params = setupGAparams(cfg.MaxGenerations, popSize);

                tic;
                out = RunGA(problem, params, specs);
                elapsed = toc;

                covStats = visualizeCameraCoverage(out, specs);
                saveResults(out, specs, params, elapsed, cfType, true, covStats);

                % Update the running best ONLY if this warm run beat it.
                if out.bestsol.Cost < bestSoFarCost
                    bestSoFarCost  = out.bestsol.Cost;
                    bestSoFarChrom = out.bestsol.Chromosome;
                    bestSoFarLabel = sprintf('warm rep %d', r);
                    improvedStr = sprintf('NEW BEST (was %.6f)', seedCost);
                else
                    improvedStr = 'no improvement';
                end

                results(end+1).CF        = cfType;       %#ok<AGROW>
                results(end).WarmStart   = true;
                results(end).Repeat      = r;
                results(end).BestCost    = out.bestsol.Cost;
                results(end).Elapsed     = elapsed;
                results(end).MinCoverage = covStats.minCoverage;
                results(end).ZeroPct     = covStats.zeroCamerasPercent;

                fprintf('    cost = %.6f | min cov = %d | 0-cov = %.1f%% | %.1fs | %s\n', ...
                    out.bestsol.Cost, covStats.minCoverage, ...
                    covStats.zeroCamerasPercent, elapsed, improvedStr);
                close all;
            end

            fprintf('=== CF%d : final running best = %s, cost = %.6f ===\n', ...
                cfType, bestSoFarLabel, bestSoFarCost);
        end
    end

    %% Final summary
    totalSec = toc(batchTic);
    fprintf('\n  BATCH COMPLETE — %.1f min wall time\n', totalSec/60);
    printSummary(results);
end


% =====================================================================
function specs = makeSpecs(numCams, cfType, cfg, volume, spacing) %#ok<INUSL>
    specs = setupHardwareSpecs(numCams);
    specs.WeightUncertainty = 0.5;
    specs.WeightOcclusion   = 0.5;
    specs.TargetType        = cfg.TargetType;
    specs.TargetMode        = cfg.GridMode;
    specs.Target            = generateTargetSpace(volume, cfg.GridMode, spacing);
    specs.NumPoints         = size(specs.Target, 1);
    specs.spacing           = spacing;
    specs.SectionCentres    = generateSectionCentres(numCams, volume);
    specs                   = setupCostParams(specs);
end


% =====================================================================
function printSummary(results)
    if isempty(results), return; end
    cfs = unique([results.CF]);
    fprintf('\n  %-4s %-6s %-8s %-9s %-9s %-9s\n', ...
        'CF', 'Type', 'Best', 'Mean', 'Std', '#runs');
    fprintf('  %s\n', repmat('-', 1, 50));
    for cfType = cfs(:)'
        for warm = [false true]
            sel = ([results.CF] == cfType) & ([results.WarmStart] == warm);
            if ~any(sel), continue; end
            costs = [results(sel).BestCost];
            fprintf('  %-4d %-6s %-8.5f %-9.5f %-9.5f %-9d\n', ...
                cfType, ternary(warm, 'warm', 'cold'), ...
                min(costs), mean(costs), std(costs), numel(costs));
        end
    end
end


function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
