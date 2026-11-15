function batchRunGA(varargin)
% Systematic batch testing across all experimental conditions.
% Number of cameras
% Cost function type (1=ResUncert, 2=DynOccl, 3=Combined)
% Target type (1=UAV, 2=UGV)
% Grid discretisation mode (1=Uniform, 2=Normal)
% Grid spacing
% Warm-start (false first, then true, seeded from the condition's seed archive)

% Make sure every code subfolder is on the MATLAB path.
addProjectPaths();

%% Inputs
% Defaults come from runConfig. Pass 'Preset' to switch setup (e.g.
% batchRunGA('Preset','lowcost_tripod')); any other runConfig field (Volume,
% Hardware, Weights, MutationRate, ...) can also be overridden by name.
presetArg = find(strcmpi(varargin(1:2:end), 'Preset'), 1);
if isempty(presetArg)
    preset = 'optitrack_lab';
else
    preset = varargin{2*presetArg};
end
base = runConfig(preset);

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'Preset', preset, @ischar);

% Experimental design
addParameter(p, 'CameraRange',   base.CameraRange,   @isnumeric);
addParameter(p, 'CostFunctions', base.CostFunctions, @isnumeric);
addParameter(p, 'TargetTypes',   base.TargetTypes,   @isnumeric); % 1=UAV, 2=UGV
addParameter(p, 'GridModes',     base.GridModes,     @isnumeric); % 1=Uniform, 2=Normal
addParameter(p, 'Spacings',      base.Spacings,      @isnumeric); % metres (x-y on the grid)
addParameter(p, 'UGV_MaxHeight', base.UGV_MaxHeight, @isnumeric);
% UGV z-axis spacing is decoupled from x-y (default 3 z-layers on a 0.5 m slab).
addParameter(p, 'UGV_ZSpacing',  base.UGV_ZSpacing,  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NumRepeats',    base.NumRepeats,    @isnumeric);
addParameter(p, 'SkipWarmStart', base.SkipWarmStart, @islogical);
% WarmOnly: schedule only the warm runs (seeded from the existing archive),
% e.g. to re-do warm runs with a new seeding scheme without repeating colds.
addParameter(p, 'WarmOnly', false, @islogical);

% GA parameters. PopulationSize [] => numCams * params-per-camera * PopulationScale.
addParameter(p, 'MaxGenerations', base.MaxGenerations, @isnumeric);
addParameter(p, 'PopulationSize', base.PopulationSize, @(x) isempty(x) || isnumeric(x));

% Workspace volume & camera search box
addParameter(p, 'Volume',         base.Volume,         @isnumeric);
addParameter(p, 'CamLowerBounds', base.CamLowerBounds, @isnumeric);
addParameter(p, 'CamUpperBounds', base.CamUpperBounds, @isnumeric);

% Execution control
addParameter(p, 'DryRun', false, @islogical);
addParameter(p, 'ResumeFrom', 1,  @isnumeric);
addParameter(p, 'ResumeLog', '', @ischar);
addParameter(p, 'OutputDir', '', @ischar);
addParameter(p, 'SuppressPlots',  true, @islogical);

parse(p, varargin{:});

% Full run configuration = preset + every override. Unmatched names must be
% runConfig fields (runConfig errors on anything else).
unmatched = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)].';
cfg = runConfig(preset, unmatched{:});
for f = fieldnames(p.Results).'
    cfg.(f{1}) = p.Results.(f{1});
end

%% Tests
% Cold-start runs first, then warm-start runs (which reference cold results).
% Phase-ordered so all CF1/CF2 runs execute before any CF3 run: the CF3
% normalisation table is derived mid-run from the completed CF1/CF2 results.
schedule = buildSchedule(cfg);
schedule = orderByPhase(schedule);
totalRuns = length(schedule);

fprintf('\n');
fprintf('  BATCH GA TESTING - %d total runs\n', totalRuns);
fprintf('Cameras: %s\n', mat2str(cfg.CameraRange));
fprintf('Cost functions: %s\n', mat2str(cfg.CostFunctions));
fprintf('Target types: %s\n', mat2str(cfg.TargetTypes));
fprintf('Grid modes: %s\n', mat2str(cfg.GridModes));
fprintf('Spacings (x-y): %s m\n', mat2str(cfg.Spacings));
if any(cfg.TargetTypes == 2)
    fprintf('UGV slab: z in [0, %.3g] m, z-spacing fixed at %.3g m (%d layers)\n', ...
        cfg.UGV_MaxHeight, cfg.UGV_ZSpacing, ...
        numel(0:cfg.UGV_ZSpacing:cfg.UGV_MaxHeight));
end
fprintf('Repeats: %d\n', cfg.NumRepeats);
fprintf('Warm-start: %s\n', string(~cfg.SkipWarmStart));
fprintf('GA generations: %d\n', cfg.MaxGenerations);
if isempty(cfg.PopulationSize)
    fprintf('Population size: auto (numCams x %d x 10)\n', numel(cfg.CamLowerBounds));
else
    fprintf('Population size: %d (fixed)\n', cfg.PopulationSize);
end

% Estimate runtime
coldRuns = sum(~[schedule.WarmStart]);
warmRuns = totalRuns - coldRuns;
fprintf('Cold-start runs: %d\n', coldRuns);
fprintf('Warm-start runs: %d\n', warmRuns);

% Print schedule preview
printSchedule(schedule);

if cfg.DryRun
    fprintf('\n DRY RUN - no optimisations executed \n');
    fprintf('To run, call: batchRunGA() without DryRun flag.\n\n');
    return;
end

%% Ouput directory
if ~isempty(cfg.OutputDir)
    if ~isfolder(cfg.OutputDir)
        mkdir(cfg.OutputDir);
    end
    originalDir = pwd;
    cd(cfg.OutputDir);
    cleanupObj = onCleanup(@() cd(originalDir));
end

%% Batch log
% Build or reload schedule
if ~isempty(cfg.ResumeLog)
    % Resume from saved batch log
    if ~isfile(cfg.ResumeLog)
        error('Resume log not found: %s', cfg.ResumeLog);
    end
    loaded = load(cfg.ResumeLog, 'schedule', 'batchConfig', 'batchTimestamp');
    schedule = loaded.schedule;
    batchTimestamp = loaded.batchTimestamp;
    batchConfig = loaded.batchConfig;
    totalRuns = length(schedule);
    
    % Find first non-done run
    doneFlags = strcmp({schedule.Status}, 'done');
    startIdx = find(~doneFlags, 1);
    if isempty(startIdx)
        fprintf('All %d runs already completed!\n', totalRuns);
        return;
    end
    
    nDone = sum(doneFlags);
    fprintf('\n  RESUMING from saved log: %s\n', cfg.ResumeLog);
    fprintf('  %d/%d runs already completed — starting from run %d\n\n', ...
        nDone, totalRuns, startIdx);
    
    batchLogFile = cfg.ResumeLog;  % keep writing to same file
else
    schedule = buildSchedule(cfg);
    schedule = orderByPhase(schedule);
    totalRuns = length(schedule);

    % ... existing initialisation code (print summary, dry run, etc.) ...
    
    for i = 1:totalRuns
        schedule(i).Status     = 'pending';
        schedule(i).BestCost   = NaN;
        schedule(i).ElapsedTime = NaN;
        schedule(i).BestChromosome = [];
        schedule(i).RunFilename = '';
        schedule(i).ErrorMsg   = '';
    end
    
    batchTimestamp = datetime('now');
    batchConfig = cfg;
    batchDateStr = string(batchTimestamp, 'yyyyMMdd_HHmmss');

    % Batch logs live alongside the master run log under Results/Logs/.
    logsDir = fullfile(addProjectPaths(), 'Results', 'Logs');
    if ~isfolder(logsDir), mkdir(logsDir); end
    batchLogFile = fullfile(logsDir, sprintf('BatchLog_%s.mat', batchDateStr));
    startIdx = max(1, cfg.ResumeFrom);
end

%% Suppress figures

if cfg.SuppressPlots
    originalVisibility = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    restorePlots = onCleanup(@() set(0, 'DefaultFigureVisible', originalVisibility));
end

%% Execution Station
batchTic = tic;

% Two-pass normalisation: CF3 runs read a normTable derived from the CF1/CF2
% runs that precede them (phase-ordered above). Built lazily the first time a
% CF3 run is reached this session, from every completed CF1/CF2 run so far.
% This also covers resume: on resume the CF1/CF2 runs are already 'done'.
normRebuilt = false;
frozenPools = containers.Map();   % ConditionKey -> archive snapshot (FreezeWarmPool)
haveCF12 = any(ismember([schedule.CostFunc], [1 2]));

for runIdx = startIdx:totalRuns
    s = schedule(runIdx);

    % Refresh CF3 normalisation from completed CF1/CF2 runs before the first
    % CF3 run. utopia = min single-objective cost (a genuine lower bound at the
    % full production budget), so CF3 = (raw - utopia)/norm stays >= 0.
    if s.CostFunc == 3 && ~normRebuilt && haveCF12
        fprintf('  Deriving CF3 normalisation from completed CF1/CF2 runs...\n');
        try
            buildNormFromBatch(schedule, batchConfig);
        catch ME
            warning('batchRunGA:normBuild', ...
                'normTable build failed (%s); CF3 will use existing/default norms.', ...
                ME.message);
        end
        normRebuilt = true;
    end

    fprintf('Run %d/%d [%s] \n', runIdx, totalRuns, ...
        string(datetime('now'), 'HH:mm:ss'));
    fprintf('    %dC | CF%d | TT%d | GM%d | sp=%.2f | WS=%d | Rep%d\n', ...
        s.NumCams, s.CostFunc, s.TargetType, s.GridMode, ...
        s.Spacing, s.WarmStart, s.RepeatIdx);
    
    schedule(runIdx).Status = 'running';
    
    try
        numCams = s.NumCams;
        costFunctionType = s.CostFunc;
        warmStartUsed = s.WarmStart;

        % Target space, hardware, cost params, bounds and GA params
        [specs, problem, params] = buildRunSpecs(cfg, numCams, costFunctionType, ...
            s.TargetType, s.GridMode, s.Spacing);

        % Warm-start chromosomes from the condition's seed archive, which
        % every completed run (cold and warm, this batch and earlier ones)
        % has added to. See WarmStart/selectWarmSeeds and runConfig.
        if warmStartUsed
            if cfg.FreezeWarmPool
                if ~isKey(frozenPools, s.ConditionKey)
                    frozenPools(s.ConditionKey) = loadSeedArchive(cfg.Preset, s.ConditionKey);
                end
                [seeds, seedInfo] = warmSeedsForRun(cfg, s.ConditionKey, problem, params, ...
                    frozenPools(s.ConditionKey));
            else
                [seeds, seedInfo] = warmSeedsForRun(cfg, s.ConditionKey, problem, params);
            end
            if ~isempty(seeds)
                specs.warmStart = true;
                specs.warmChromosomes = seeds;
                specs.WarmSeedInfo = seedInfo;
                schedule(runIdx).WarmSeedChromosome = seeds;
            end
        end
        
        %% Run GA 
        tic;
        out = RunGA(problem, params, specs);
        elapsedTime = toc;
        
        %% Post-process 
        coverageStats = visualizeCameraCoverage(out, specs);
        
        if ~cfg.SuppressPlots
            plotResults(out, specs, params, elapsedTime);
        end
        
        % Use existing saveResults to maintain log compatibility
        runFile = saveResults(out, specs, params, elapsedTime, ...
            costFunctionType, warmStartUsed, coverageStats);
        
        %% Record in batch schedule 
        schedule(runIdx).Status = 'done';
        schedule(runIdx).BestCost = out.bestsol.Cost;
        schedule(runIdx).ElapsedTime = elapsedTime;
        schedule(runIdx).BestChromosome = out.bestsol.Chromosome;
        schedule(runIdx).CoverageStats = coverageStats;
        
        schedule(runIdx).RunFilename = runFile;
        
        fprintf('DONE: Cost=%.6f | Time=%.1fs (%.1f min)\n\n', ...
            out.bestsol.Cost, elapsedTime, elapsedTime/60);
        
    catch ME
        schedule(runIdx).Status = 'failed';
        schedule(runIdx).ErrorMsg = ME.message;
        fprintf('    FAILED: %s\n\n', ME.message);
    end
    
    % Save progress after every run
    save(batchLogFile, 'schedule', 'batchConfig', 'batchTimestamp');
    
    % Close all figures to prevent memory overoad
    close all;
    
    % Progress summary
    elapsed = toc(batchTic);
    completed = sum(strcmp({schedule.Status}, 'done'));
    remaining = totalRuns - runIdx;
    if completed > 0
        avgTime = elapsed / (runIdx - startIdx + 1);
        etaSeconds = avgTime * remaining;
        fprintf('    Progress: %d/%d done | ETA: %.1f hours\n\n', ...
            completed, totalRuns, etaSeconds/3600);
    end
end

%% FINAL SUMMARY
totalElapsed = toc(batchTic);
printBatchSummary(schedule, totalElapsed, batchLogFile);

end


%% Scheduler

function schedule = buildSchedule(cfg)
% Full factorial batch. WarmOnly keeps only the warm runs.
    schedule = buildFullSchedule(cfg);
    if isfield(cfg, 'WarmOnly') && cfg.WarmOnly
        if cfg.SkipWarmStart
            error('batchRunGA:WarmOnlyNoWarm', 'WarmOnly with SkipWarmStart leaves nothing to run.');
        end
        schedule = schedule([schedule.WarmStart]);
        for i = 1:numel(schedule), schedule(i).RunIndex = i; end
    end
end

function schedule = buildFullSchedule(cfg)

    schedule = struct([]);
    idx = 0;
    
    % Phase 1: All cold-start runs
    for nc = cfg.CameraRange
        for cf = cfg.CostFunctions
            for tt = cfg.TargetTypes
                for gm = cfg.GridModes
                    for sp = cfg.Spacings
                        for rep = 1:cfg.NumRepeats
                            idx = idx + 1;
                            schedule(idx).RunIndex    = idx;
                            schedule(idx).NumCams     = nc;
                            schedule(idx).CostFunc    = cf;
                            schedule(idx).TargetType  = tt;
                            schedule(idx).GridMode    = gm;
                            schedule(idx).Spacing     = sp;
                            schedule(idx).WarmStart   = false;
                            schedule(idx).RepeatIdx   = rep;
                            schedule(idx).WarmSeedChromosome = [];
                            schedule(idx).ConditionKey = sprintf( ...
                                '%dC_CF%d_TT%d_GM%d_sp%.2f', ...
                                nc, cf, tt, gm, sp);
                        end
                    end
                end
            end
        end
    end
    
    % Phase 2: Warm-start runs (NumRepeats per condition, seeded at run time
    % from the seed archive, which by then holds this batch's cold runs)
    if ~cfg.SkipWarmStart
        uniqueConditions = unique({schedule.ConditionKey});
        for c = 1:length(uniqueConditions)
            condKey = uniqueConditions{c};
            % Parse condition parameters from key
            ref = schedule(strcmp({schedule.ConditionKey}, condKey));
            ref = ref(1); % grab first for parameters
            
            for rep = 1:cfg.NumRepeats
                idx = idx + 1;
                schedule(idx).RunIndex    = idx;
                schedule(idx).NumCams     = ref.NumCams;
                schedule(idx).CostFunc    = ref.CostFunc;
                schedule(idx).TargetType  = ref.TargetType;
                schedule(idx).GridMode    = ref.GridMode;
                schedule(idx).Spacing     = ref.Spacing;
                schedule(idx).WarmStart   = true;
                schedule(idx).RepeatIdx   = rep;
                schedule(idx).WarmSeedChromosome = []; % filled at runtime
                schedule(idx).ConditionKey = condKey;
            end
        end
    end
end


function printSchedule(schedule)
% Print a compact preview of the test schedule.
    
    coldRuns = schedule(~[schedule.WarmStart]);
    warmRuns = schedule([schedule.WarmStart]);
    
    fprintf('--- COLD-START RUNS ---\n');
    fprintf('%-5s %-6s %-4s %-4s %-4s %-7s %-4s\n', ...
        'Run', 'Cams', 'CF', 'TT', 'GM', 'Space', 'Rep');
    fprintf('%s\n', repmat('-', 1, 40));
    
    % Show first few and last few if many
    nCold = length(coldRuns);
    if nCold <= 20
        showIdx = 1:nCold;
    else
        showIdx = [1:10, nCold-4:nCold];
    end
    
    for i = showIdx
        s = coldRuns(i);
        fprintf('%-5d %-6d %-4d %-4d %-4d %-7.2f %-4d\n', ...
            s.RunIndex, s.NumCams, s.CostFunc, s.TargetType, ...
            s.GridMode, s.Spacing, s.RepeatIdx);
        if i == 10 && nCold > 20
            fprintf('  ... (%d rows omitted) ...\n', nCold - 15);
        end
    end
    
    if ~isempty(warmRuns)
        fprintf('\n--- WARM-START RUNS ---\n');
        nWarm = length(warmRuns);
        fprintf('  %d warm-start runs (seeded from best cold run per condition)\n', nWarm);
    end
    
    fprintf('\n');
end


function printBatchSummary(schedule, totalElapsed, logFile)
% Print final batch summary with per-condition statistics.

    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  BATCH COMPLETE\n');
    fprintf('=========================================================\n');
    fprintf('  Total wall time: %.1f hours (%.1f min)\n', ...
        totalElapsed/3600, totalElapsed/60);
    
    nDone   = sum(strcmp({schedule.Status}, 'done'));
    nFailed = sum(strcmp({schedule.Status}, 'failed'));
    nPending = sum(strcmp({schedule.Status}, 'pending'));
    
    fprintf('  Completed: %d | Failed: %d | Pending: %d\n', ...
        nDone, nFailed, nPending);
    fprintf('  Results log: %s\n', logFile);
    
    % Per-condition summary table
    if nDone > 0
        doneIdx = strcmp({schedule.Status}, 'done');
        doneRuns = schedule(doneIdx);
        
        conditions = unique({doneRuns.ConditionKey});
        
        fprintf('\n  %-25s  %-8s %-8s %-8s %-10s %-6s\n', ...
            'Condition', 'Cold', 'Warm', 'Best', 'Avg Time', 'Runs');
        fprintf('  %s\n', repmat('-', 1, 75));
        
        for c = 1:length(conditions)
            condKey = conditions{c};
            condRuns = doneRuns(strcmp({doneRuns.ConditionKey}, condKey));
            
            coldCosts = [condRuns(~[condRuns.WarmStart]).BestCost];
            warmCosts = [condRuns([condRuns.WarmStart]).BestCost];
            allCosts  = [condRuns.BestCost];
            allTimes  = [condRuns.ElapsedTime];
            
            coldBest = min(coldCosts);
            warmBest = min(warmCosts);
            overallBest = min(allCosts);
            avgTime = mean(allTimes);
            
            coldStr = sprintf('%.4f', coldBest);
            if isempty(warmCosts)
                warmStr = '  -   ';
            else
                warmStr = sprintf('%.4f', warmBest);
            end
            
            fprintf('  %-25s  %-8s %-8s %-8.4f %-10.1fs %-6d\n', ...
                condKey, coldStr, warmStr, overallBest, avgTime, length(condRuns));
        end
    end
    
    % Report any failures
    if nFailed > 0
        fprintf('\n  FAILED RUNS:\n');
        failedRuns = schedule(strcmp({schedule.Status}, 'failed'));
        for i = 1:length(failedRuns)
            fprintf('    Run %d (%s): %s\n', ...
                failedRuns(i).RunIndex, failedRuns(i).ConditionKey, ...
                failedRuns(i).ErrorMsg);
        end
    end
    
    fprintf('\n=========================================================\n\n');
end


function schedule = orderByPhase(schedule)
% Execute all CF1/CF2 runs before any CF3 run so the CF3 normalisation can be
% derived from completed single-objective runs. Stable within each phase, so
% the existing cold-before-warm ordering is preserved. RunIndex is renumbered
% to match the new execution order.
    if isempty(schedule)
        return;
    end
    phase = 1 + ([schedule.CostFunc] == 3);   % 1 = CF1/CF2, 2 = CF3
    [~, ord] = sort(phase);   % MATLAB sort is stable: ties keep input order
    schedule = schedule(ord);
    for i = 1:numel(schedule)
        schedule(i).RunIndex = i;
    end
end