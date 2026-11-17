function summary = backfillSeedArchive(varargin)
%BACKFILLSEEDARCHIVE  Fill the seed archive from runs already on disk.
%
%   summary = backfillSeedArchive()                       % lab preset, every run
%   summary = backfillSeedArchive('Preset', 'openmv_2m')
%   summary = backfillSeedArchive('Since', datetime(2026,10,5))
%
%   Adds each logged run's best chromosome (and its FinalTopK, for runs saved
%   after the archive existed) to Results/SeedArchive/<preset>/<condition>.mat.
%   Only rotation-vector runs are taken (older Euler-gene runs were scored by
%   superseded cost code). Already-archived chromosomes are skipped, so this
%   is safe to re-run.
%
%   'Since' defaults to 2026-10-05, the start of the re-run with the fixed
%   occlusion/mount/rotvec code; earlier runs used other cost surfaces.

    addProjectPaths();
    p = inputParser;
    addParameter(p, 'Preset', 'optitrack_lab', @ischar);
    addParameter(p, 'Since',  datetime(2026, 10, 5), @(x) isdatetime(x));
    parse(p, varargin{:});
    o = p.Results;

    runLog = loadGARuns('Preset', o.Preset);
    added = containers.Map();
    nSkip = 0;
    for i = 1:numel(runLog)
        r = runLog(i);
        if r.Timestamp < o.Since, nSkip = nSkip + 1; continue; end
        f = resolveRunPath(r.RunFilename, r.NumCameras);
        if ~isfile(f), nSkip = nSkip + 1; continue; end
        S = load(f, 'saveData');
        d = S.saveData;
        if ~isfield(d.Specifications, 'Parameterisation') || ...
                ~strcmp(d.Specifications.Parameterisation, 'rotvec')
            nSkip = nSkip + 1; continue;
        end
        key = conditionKey(r.NumCameras, r.CostFunctionType, r.TargetType, r.GridMode, r.Spacing);
        C = d.BestSolution.Chromosome;  J = d.BestCost;  src = {'best'};
        if isfield(d, 'FinalTopK')
            C = [C; d.FinalTopK.Chromosomes(2:end, :)]; %#ok<AGROW>
            J = [J; d.FinalTopK.Costs(2:end)];          %#ok<AGROW>
            src = [src; repmat({'finalpop'}, numel(J) - 1, 1)]; %#ok<AGROW>
        end
        n = updateSeedArchive(o.Preset, key, C, J, src, r.RunFilename, r.WarmStart);
        if isKey(added, key), added(key) = added(key) + n; else, added(key) = n; end
    end

    keys_ = added.keys;
    summary = table(keys_(:), cell2mat(added.values(keys_)).', ...
        'VariableNames', {'Condition', 'Added'});
    for k = 1:height(summary)
        summary.PoolSize(k) = numel(loadSeedArchive(o.Preset, summary.Condition{k}));
        pool = loadSeedArchive(o.Preset, summary.Condition{k});
        summary.BestCost(k) = pool(1).Cost;
    end
    fprintf('Backfill (%s): %d runs skipped (older / missing / Euler genes)\n', o.Preset, nSkip);
    disp(summary);
end
