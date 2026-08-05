function results = surveyGACoverage(varargin)
% SURVEYGACOVERAGE  Coverage stats for every 7-cam CF3 GA run, by scenario.
%
%   results = surveyGACoverage()
%   results = surveyGACoverage('NumCameras', 7, 'CostFunction', 3)
%
%   Loops over all matching GA runs in the log, computes per-target camera
%   coverage for each (using perTargetCoverage), groups them by scenario
%   (TargetType / GridMode / Spacing) and prints a table per scenario
%   sorted by AVERAGE coverage (highest first). The OptiTrack ad-hoc rig is
%   appended to each scenario for reference.
%
%   This exists because plotHeatmap_GAvsOptiTrack only ever shows the
%   LOWEST-COST run (min BestCost), which is not necessarily the run with
%   the best coverage. Use this to see which run actually maximises 2+ /
%   average coverage.
%
%   Name-Value Parameters (all forwarded to loadGARuns):
%     'NumCameras'   - default 7
%     'CostFunction' - default 3 (combined). Pass [] for all cost functions.
%     'LogFile'      - default Results/Logs/GGA_RunsLog.mat
%
%   OUTPUT
%     results - struct array, one row per run, with fields:
%               RunFilename, TargetType, GridMode, Spacing, BestCost,
%               Avg, ZeroPct, OnePct, TwoPlusPct
%
%   See also: perTargetCoverage, plotHeatmap_GAvsOptiTrack, loadGARuns

    addProjectPaths();

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'NumCameras',   7);
    addParameter(p, 'CostFunction', 3);
    parse(p, varargin{:});
    numCams = p.Results.NumCameras;
    cf      = p.Results.CostFunction;

    % Pull matching runs from the log.
    loadArgs = {'NumCameras', numCams};
    if ~isempty(cf)
        loadArgs = [loadArgs, {'CostFunction', cf}];
    end
    runLog = loadGARuns(loadArgs{:});

    % Backfill grouping fields if older rows lack them.
    for fn = {'TargetType', 'GridMode', 'Spacing'}
        if ~isfield(runLog, fn{1}), [runLog.(fn{1})] = deal(NaN); end
    end

    nRuns = numel(runLog);
    results = struct('RunFilename', {}, 'TargetType', {}, 'GridMode', {}, ...
        'Spacing', {}, 'BestCost', {}, 'Avg', {}, 'ZeroPct', {}, ...
        'OnePct', {}, 'TwoPlusPct', {});

    fprintf('Computing coverage for %d run(s)...\n', nRuns);
    for i = 1:nRuns
        r = runLog(i);
        try
            matPath = resolveRunPath(r.RunFilename, r.NumCameras);
            L  = load(matPath, 'saveData');
            sd = L.saveData;
            [~, st] = perTargetCoverage(sd.BestSolution.Chromosome, sd.Specifications);
        catch ME
            fprintf('  Skipped %s: %s\n', r.RunFilename, ME.message);
            continue;
        end
        results(end+1) = struct( ...
            'RunFilename', r.RunFilename, ...
            'TargetType',  r.TargetType, ...
            'GridMode',    r.GridMode, ...
            'Spacing',     r.Spacing, ...
            'BestCost',    r.BestCost, ...
            'Avg',         st.avg, ...
            'ZeroPct',     st.zeroPct, ...
            'OnePct',      st.onePct, ...
            'TwoPlusPct',  st.twoPlusPct); %#ok<AGROW>
    end

    if isempty(results)
        warning('No runs produced coverage stats.');
        return;
    end

    % Build scenario keys and print one sorted table per scenario.
    keys = arrayfun(@(x) sprintf('TT%d_GM%d_sp%.2f', ...
        x.TargetType, x.GridMode, x.Spacing), results, 'UniformOutput', false);
    [uKeys, ~, grp] = unique(keys);

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    for g = 1:numel(uKeys)
        sub = results(grp == g);
        ex  = sub(1);
        if ex.TargetType >= 1 && ex.TargetType <= 2, ttStr = ttNames{ex.TargetType}; else, ttStr = '?'; end
        if ex.GridMode   >= 1 && ex.GridMode   <= 2, gmStr = gmNames{ex.GridMode};   else, gmStr = '?'; end

        % Sort runs by average coverage, descending.
        [~, ord] = sort([sub.Avg], 'descend');
        sub = sub(ord);

        fprintf('\n==========================================================\n');
        fprintf('  %s | %s grid | sp = %.2f m | CF%s | %d runs\n', ...
            ttStr, gmStr, ex.Spacing, num2str(cf), numel(sub));
        fprintf('==========================================================\n');
        fprintf('  %-32s %8s %7s %7s %7s %7s\n', ...
            'RunFilename', 'Cost', 'Avg', '0-cam', '1-cam', '2+cam');
        fprintf('  %s\n', repmat('-', 1, 78));
        for j = 1:numel(sub)
            fprintf('  %-32s %8.4f %7.2f %6.1f%% %6.1f%% %6.1f%%\n', ...
                sub(j).RunFilename, sub(j).BestCost, sub(j).Avg, ...
                sub(j).ZeroPct, sub(j).OnePct, sub(j).TwoPlusPct);
        end

        % OptiTrack reference for this scenario.
        try
            optiChrom = buildOptiTrackChromosome();
            matPath = resolveRunPath(sub(1).RunFilename, numCams);
            sd = load(matPath, 'saveData').saveData;     % borrow specs.Target
            [~, ost] = perTargetCoverage(optiChrom, sd.Specifications);
            fprintf('  %s\n', repmat('-', 1, 78));
            fprintf('  %-32s %8s %7.2f %6.1f%% %6.1f%% %6.1f%%\n', ...
                'OptiTrack ad-hoc (reference)', '--', ost.avg, ...
                ost.zeroPct, ost.onePct, ost.twoPlusPct);
        catch ME
            fprintf('  (OptiTrack reference unavailable: %s)\n', ME.message);
        end
    end
    fprintf('\n');
end
