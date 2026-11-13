function R = abWarmStart(varargin)
%ABWARMSTART  A/B study of the warm-start seeding scheme.
%
%   R = abWarmStart()                         % full study (sim machine)
%   R = abWarmStart('Cameras', 7)             % 7C only, half the time
%   R = abWarmStart('Quick', true)            % smoke test, minutes
%   R = abWarmStart('Resume', 'Results/Tuning/ABWarmStart/<stamp>/abWarmStart.mat')
%
%   Design (per camera count; UAV, uniform grid, CF3, lab preset):
%     Step 0  'Step0Runs' fresh cold runs, saved normally (saveResults), so
%             they enter the run log and the seed archive together with
%             their final populations. Gives the pool > 20 candidates.
%     Main    2 x 2 arms, every arm seeded from the SAME frozen archive
%             snapshot taken after step 0 (arms never feed each other):
%               method  'best_mutant' (best + mutated copies, old scheme)
%                       'top_n'       (best distinct archived solutions)
%               seeds   'SeedCounts' (default [10 21]; 21 = 5% of 420)
%             'Repeats' paired repeats: repeat r of every arm uses the same
%             RNG seed for the GA, so the guided-random part of the initial
%             population and the GA's random draws are shared across arms.
%   Main-study runs are NOT added to the archive or the master run log.
%
%   Output: Results/Tuning/ABWarmStart/<stamp>/abWarmStart.mat (saved after
%   every run, resumable) and the report from reportABWarmStart.

    root = addProjectPaths();
    p = inputParser;
    addParameter(p, 'Preset',     'optitrack_lab', @ischar);
    addParameter(p, 'Cameras',    [7 6],  @isnumeric);
    addParameter(p, 'CostFunc',   3,      @isnumeric);
    addParameter(p, 'TargetType', 1,      @isnumeric);
    addParameter(p, 'GridMode',   1,      @isnumeric);
    addParameter(p, 'Spacing',    1.0,    @isnumeric);
    addParameter(p, 'Repeats',    15,     @isnumeric);
    addParameter(p, 'SeedCounts', [10 21], @isnumeric);
    addParameter(p, 'Methods',    {'best_mutant', 'top_n'}, @iscell);
    addParameter(p, 'MinDist',    0.25,   @isnumeric);
    addParameter(p, 'Step0Runs',  5,      @isnumeric);
    addParameter(p, 'Quick',      false,  @islogical);
    addParameter(p, 'Resume',     '',     @ischar);
    parse(p, varargin{:});
    o = p.Results;

    cfgArgs = {};
    if o.Quick
        % Read-only smoke test: no step 0 (nothing written to the log or
        % archive), tiny GA, two repeats.
        o.Cameras = o.Cameras(end);  o.Repeats = 2;  o.Step0Runs = 0;
        o.SeedCounts = [2 4];
        cfgArgs = {'MaxGenerations', 3, 'PopulationSize', 24};
    end
    cfg = runConfig(o.Preset, cfgArgs{:});

    %% Output / resume
    if ~isempty(o.Resume)
        S = load(o.Resume, 'R');
        R = S.R;  outFile = o.Resume;
        fprintf('Resuming %s: %d of %d runs done\n', outFile, sum(~isnan([R.runs.BestCost])), numel(R.runs));
    else
        stamp = string(datetime('now'), 'yyyyMMdd_HHmmss');
        if o.Quick, stamp = stamp + "_quick"; end
        outDir = fullfile(root, 'Results', 'Tuning', 'ABWarmStart', char(stamp));
        mkdir(outDir);
        outFile = fullfile(outDir, 'abWarmStart.mat');
        R.opts = o;  R.cfg = cfg;  R.started = datetime('now');
        R.step0 = struct('NumCams', {}, 'RunFile', {}, 'BestCost', {});
        R.snapshot = struct('NumCams', {}, 'Key', {}, 'Pool', {});
        R.runs = buildRuns(o);
        R.step0Done = false;
        save(outFile, 'R');
    end
    o = R.opts;  cfg = R.cfg;

    %% Step 0: fresh cold runs into the archive, then freeze the snapshot
    if ~R.step0Done
        if ~o.Quick
            backfillSeedArchive('Preset', o.Preset);   % idempotent
        end
        for nc = o.Cameras
            for k = 1:o.Step0Runs
                fprintf('\n[step 0] %dC cold run %d/%d\n', nc, k, o.Step0Runs);
                [specs, problem, params] = buildRunSpecs(cfg, nc, o.CostFunc, o.TargetType, o.GridMode, o.Spacing);
                t = tic;
                out = RunGA(problem, params, specs);
                el = toc(t);
                cov = visualizeCameraCoverage(out, specs);
                close all;
                f = saveResults(out, specs, params, el, o.CostFunc, false, cov);
                R.step0(end+1) = struct('NumCams', nc, 'RunFile', f, 'BestCost', out.bestsol.Cost);
                save(outFile, 'R');
            end
        end
        for nc = o.Cameras
            key = conditionKey(nc, o.CostFunc, o.TargetType, o.GridMode, o.Spacing);
            pool = loadSeedArchive(o.Preset, key);
            if isempty(pool)
                error('abWarmStart:EmptyArchive', ...
                    'No archive for %s. Run backfillSeedArchive first.', key);
            end
            R.snapshot(end+1) = struct('NumCams', nc, 'Key', key, 'Pool', pool);
            fprintf('Snapshot %s: %d archived solutions, best %.4f\n', key, numel(pool), pool(1).Cost);
        end
        R.step0Done = true;
        save(outFile, 'R');
    end

    %% Main study
    if ~o.Quick
        set(0, 'DefaultFigureVisible', 'off');
    end
    todo = find(isnan([R.runs.BestCost]));
    for n = 1:numel(todo)
        i = todo(n);
        r = R.runs(i);
        fprintf('\n[%d/%d] %dC  %-11s  %2d seeds  rep %d\n', ...
            numel(R.runs) - numel(todo) + n, numel(R.runs), r.NumCams, r.Method, r.NumSeeds, r.Repeat);

        [specs, problem, params] = buildRunSpecs(cfg, r.NumCams, o.CostFunc, o.TargetType, o.GridMode, o.Spacing);
        pool = R.snapshot([R.snapshot.NumCams] == r.NumCams).Pool;

        % Mutants depend on the repeat only, not on the arm order
        rng(r.RngSeed, 'twister');
        [seeds, info] = selectWarmSeeds(pool, r.NumSeeds, r.Method, o.MinDist, problem);
        specs.warmStart = true;
        specs.warmChromosomes = seeds;

        rng(r.RngSeed + 1e6, 'twister');     % same GA stream for every arm of this repeat
        t = tic;
        out = RunGA(problem, params, specs);
        el = toc(t);

        R.runs(i).BestCost      = out.bestsol.Cost;
        R.runs(i).BestChromosome = out.bestsol.Chromosome;
        R.runs(i).Convergence   = out.bestcost(:).';
        R.runs(i).InitialBest   = min(info.Costs);
        R.runs(i).SeedInfo      = info;
        R.runs(i).Elapsed       = el;
        fprintf('    best %.5f (best seed %.5f) in %.1f min\n', out.bestsol.Cost, min(info.Costs), el/60);
        save(outFile, 'R');
    end
    R.finished = datetime('now');
    save(outFile, 'R');

    reportABWarmStart(outFile);
end

function runs = buildRuns(o)
    runs = struct('NumCams', {}, 'Method', {}, 'NumSeeds', {}, 'Repeat', {}, 'RngSeed', {}, ...
        'BestCost', {}, 'BestChromosome', {}, 'Convergence', {}, 'InitialBest', {}, ...
        'SeedInfo', {}, 'Elapsed', {});
    % Interleave arms within a repeat so a stopped study is still balanced
    for nc = o.Cameras
        for rep = 1:o.Repeats
            for m = 1:numel(o.Methods)
                for s = o.SeedCounts
                    runs(end+1) = struct('NumCams', nc, 'Method', o.Methods{m}, 'NumSeeds', s, ...
                        'Repeat', rep, 'RngSeed', 1000*nc + rep, 'BestCost', NaN, ...
                        'BestChromosome', [], 'Convergence', [], 'InitialBest', NaN, ...
                        'SeedInfo', [], 'Elapsed', NaN); %#ok<AGROW>
                end
            end
        end
    end
end
