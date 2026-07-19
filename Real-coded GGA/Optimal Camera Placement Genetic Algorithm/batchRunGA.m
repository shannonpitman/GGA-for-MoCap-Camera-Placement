function batchRunGA(varargin)
% Systematic batch testing across all experimental conditions.
% Number of cameras
% Cost function type (1=ResUncert, 2=DynOccl, 3=Combined)
% Target type (1=UAV, 2=UGV)
% Grid discretisation mode (1=Uniform, 2=Normal)
% Grid spacing
% Warm-start (false first, then true seeded from best cold run)

% Make sure every code subfolder is on the MATLAB path.
addProjectPaths();

%% Inputs
p = inputParser;

% Experimental design
addParameter(p, 'CameraRange', [6 7 8], @isnumeric);
addParameter(p, 'CostFunctions', [1 2 3], @isnumeric);
addParameter(p, 'TargetTypes', [1 2], @isnumeric); % 1=UAV, 2=UGV
addParameter(p, 'GridModes', [1 2], @isnumeric); % 1=Uniform, 2=Normal
addParameter(p, 'Spacings', [1.0], @isnumeric); % metres (x-y on the grid)
addParameter(p, 'UGV_MaxHeight', 0.5, @isnumeric);
% UGV z-axis spacing is decoupled from x-y. Default 0.25 m gives 3 z-layers
% on a 0.5 m slab (z = [0, 0.25, 0.5]). x-y spacing is governed by Spacings
% so UGV x-y can match UAV without forcing 33x33 in-plane resolution.
addParameter(p, 'UGV_ZSpacing', 0.25, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NumRepeats', 5, @isnumeric);
addParameter(p, 'SkipWarmStart', false, @islogical);

% GA parameters
% PopulationSize [] => auto-scale per run as numCams * params-per-camera * 10
% (= problem.nVar * 10). Pass a numeric value to force a fixed population.
addParameter(p, 'MaxGenerations', 100,  @isnumeric);
addParameter(p, 'PopulationSize', [],   @(x) isempty(x) || isnumeric(x));

% Workspace volume & mounting constraints
addParameter(p, 'Volume', [-4 4; -4 4; 0 4], @isnumeric);
addParameter(p, 'CamLowerBounds', [-5 -4.5 0  -pi -pi/2 -pi], @isnumeric);
addParameter(p, 'CamUpperBounds', [ 5  4.5 4.8 pi  pi/2  pi], @isnumeric);

% Execution control
addParameter(p, 'DryRun', false, @islogical);
addParameter(p, 'ResumeFrom', 1,  @isnumeric);
addParameter(p, 'ResumeLog', '', @ischar);
addParameter(p, 'OutputDir', '', @ischar);
addParameter(p, 'SuppressPlots',  true, @islogical);

parse(p, varargin{:});
cfg = p.Results;

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
    logsDir = fullfile(fileparts(mfilename('fullpath')), 'Results', 'Logs');
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
        targetType = s.TargetType;
        targetMode = s.GridMode;
        spacing = s.Spacing;
        warmStartUsed = s.WarmStart;
        volume = cfg.Volume;
        if targetType == 2
            % UGV floor-slab volume
            volume(3,:) = [0, cfg.UGV_MaxHeight];
            % Anisotropic spacing for UGV: x-y honour the swept Spacing,
            % z is held fixed (default 3 layers on a 0.5 m slab). This was
            % previously clamped to UGV_MaxHeight/2 isotropic, which forced
            % the in-plane grid to 0.25 m and made UGV runs ~7x slower than
            % UAV. See generateTargetSpace for the vector-spacing contract.
            zSp = min(cfg.UGV_ZSpacing, cfg.UGV_MaxHeight); % safety clamp
            targetSpacing = [spacing, spacing, zSp];   % vector, for grid build
        else
            targetSpacing = spacing;                   % scalar, isotropic
        end
        
        specs = setupHardwareSpecs(numCams);
        
        % Warm-start chromosomes
        if warmStartUsed && ~isempty(s.WarmSeedChromosome)
            seedChrom = s.WarmSeedChromosome;
            problem_temp = setupProblem(numCams, costFunctionType, ...
                cfg.CamUpperBounds, cfg.CamLowerBounds);
            perturbed = Mutate(seedChrom, 1, 0.5);
            perturbed = max(perturbed, problem_temp.VarMin);
            perturbed = min(perturbed, problem_temp.VarMax);
            specs.warmStart = true;
            specs.warmChromosomes = [seedChrom; perturbed];
        else
            specs.warmStart = false;
            specs.warmChromosomes = [];
        end
        
        % Cost function weights
        specs.WeightUncertainty = 0.5;
        specs.WeightOcclusion   = 0.5;
        
        % Target space
        specs.TargetType = targetType;
        specs.TargetMode = targetMode;
        specs.Target = generateTargetSpace(volume, targetMode, targetSpacing);
        specs.NumPoints = size(specs.Target, 1);
        specs.spacing = spacing;                 % scalar x-y; kept for log
        if targetType == 2
            specs.spacingZ = zSp;                % UGV z spacing (info only)
        end
        
        % Section centres
        specs.SectionCentres = generateSectionCentres(numCams, volume);
        
        % Cost parameters
        specs = setupCostParams(specs);
        
        % Problem definition
        problem = setupProblem(numCams, costFunctionType, ...
            cfg.CamUpperBounds, cfg.CamLowerBounds);
        
        % GA parameters. Population scales with chromosome length unless a
        % fixed PopulationSize was supplied: nPop = numCams * numParams * 10.
        numParams = numel(cfg.CamLowerBounds);
        if isempty(cfg.PopulationSize)
            popSize = numCams * numParams * 10;
        else
            popSize = cfg.PopulationSize;
        end
        params = setupGAparams(cfg.MaxGenerations, popSize);
        
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
        saveResults(out, specs, params, elapsedTime, ...
            costFunctionType, warmStartUsed, coverageStats);
        
        %% Record in batch schedule 
        schedule(runIdx).Status = 'done';
        schedule(runIdx).BestCost = out.bestsol.Cost;
        schedule(runIdx).ElapsedTime = elapsedTime;
        schedule(runIdx).BestChromosome = out.bestsol.Chromosome;
        schedule(runIdx).CoverageStats = coverageStats;
        
        latestMat = dir(sprintf('%dCams_Run_*.mat', numCams));
        if ~isempty(latestMat)
            [~, newest] = max([latestMat.datenum]);
            schedule(runIdx).RunFilename = latestMat(newest).name;
        end
        
        fprintf('DONE: Cost=%.6f | Time=%.1fs (%.1f min)\n\n', ...
            out.bestsol.Cost, elapsedTime, elapsedTime/60);
        
        %% Propagate best chromosome to warm-start runs
        % If this is a cold-start run, find any warm-start runs that
        % share the same condition and update their seed chromosome
        if ~warmStartUsed
            schedule = propagateWarmSeed(schedule, runIdx);
        end
        
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
% Full factorial batch

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
    
    % Phase 2: Warm-start runs (one per condition, seeded from best cold)
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


function schedule = propagateWarmSeed(schedule, completedIdx)
% After a cold-start run completes, find warm-start runs with the same
% condition key. If the completed run's cost is better than any previously
% stored seed, update the seed chromosome.

    condKey = schedule(completedIdx).ConditionKey;
    newCost = schedule(completedIdx).BestCost;
    newChrom = schedule(completedIdx).BestChromosome;
    
    for i = 1:length(schedule)
        if schedule(i).WarmStart && ...
           strcmp(schedule(i).ConditionKey, condKey) && ...
           strcmp(schedule(i).Status, 'pending')
            
            % Only update if this chromosome is better than current seed
            if isempty(schedule(i).WarmSeedChromosome)
                schedule(i).WarmSeedChromosome = newChrom;
            else
                % Check if we already have a cost stored; compare
                % For simplicity, always take the latest completed cold run's best
                % (the last cold run to finish will have seen whether warm-seeding
                %  from the best of N repeats matters)
                coldDone = find(~[schedule.WarmStart] & ...
                    strcmp({schedule.ConditionKey}, condKey) & ...
                    strcmp({schedule.Status}, 'done'));
                if ~isempty(coldDone)
                    coldCosts = [schedule(coldDone).BestCost];
                    [~, bestColdIdx] = min(coldCosts);
                    schedule(i).WarmSeedChromosome = ...
                        schedule(coldDone(bestColdIdx)).BestChromosome;
                end
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