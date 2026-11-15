function T = compareGenerations(varargin)
%COMPAREGENERATIONS  Generation-100 vs final (150) cost, and where the gain is.
%
%   T = compareGenerations()                          % lab preset, CF3, runs since 2026-11-15
%   T = compareGenerations('Preset', 'openmv_lab')
%   T = compareGenerations('Since', datetime(2026,11,15), 'WarmStart', true)
%
%   Uses runs saved with a Snapshot (saveResults, SnapshotGeneration). For
%   each condition (cameras x target x grid):
%     - median CF3 cost at the snapshot generation and at the end, the
%       median % improvement and a paired Wilcoxon signed-rank p-value
%     - how the improvement splits between the weighted resolution term
%       (Junc) and the weighted occlusion term (Jocc): median dJunc, dJocc
%       and the share of the total drop due to each
%     - for comparison, the median of the October warm runs (best + 1
%       mutant, 100 generations) of the same condition, if any
%   Writes Results/Analysis/compareGenerations_<preset>_<date>.txt and
%   returns the table.

    root = addProjectPaths();
    p = inputParser;
    addParameter(p, 'Preset',    'optitrack_lab', @ischar);
    addParameter(p, 'CostFunc',  3, @isnumeric);
    addParameter(p, 'Since',     datetime(2026, 11, 15), @isdatetime);
    addParameter(p, 'WarmStart', [], @(x) isempty(x) || islogical(x));
    parse(p, varargin{:});
    o = p.Results;

    args = {'Preset', o.Preset, 'CostFunction', o.CostFunc};
    if ~isempty(o.WarmStart), args = [args, {'WarmStart', o.WarmStart}]; end
    L = loadGARuns(args{:});
    L = L([L.Timestamp] >= o.Since);

    rows = struct('Key', {}, 'Warm', {}, 'J100', {}, 'J150', {}, 'Junc100', {}, 'Junc150', {}, ...
                  'Jocc100', {}, 'Jocc150', {});
    for i = 1:numel(L)
        S = load(resolveRunPath(L(i).RunFilename, L(i).NumCameras), 'saveData');
        d = S.saveData;
        if ~isfield(d, 'Snapshot') || ~isfield(d, 'Breakdown'), continue; end
        rows(end+1) = struct('Key', conditionKey(L(i).NumCameras, L(i).CostFunctionType, ...
                L(i).TargetType, L(i).GridMode, L(i).Spacing), 'Warm', L(i).WarmStart, ...
            'J100', d.Breakdown.Snapshot.J, 'J150', d.Breakdown.Final.J, ...
            'Junc100', d.Breakdown.Snapshot.Junc, 'Junc150', d.Breakdown.Final.Junc, ...
            'Jocc100', d.Breakdown.Snapshot.Jocc, 'Jocc150', d.Breakdown.Final.Jocc); %#ok<AGROW>
    end
    if isempty(rows)
        error('compareGenerations:NoRuns', 'No runs with a generation snapshot since %s.', string(o.Since));
    end

    % October warm runs (old seeding, 100 generations) for reference
    try
        Lold = loadGARuns('Preset', o.Preset, 'CostFunction', o.CostFunc, 'WarmStart', true);
        Lold = Lold([Lold.Timestamp] < datetime(2026, 11, 1) & [Lold.Timestamp] >= datetime(2026, 10, 5));
    catch
        Lold = [];
    end

    keys = unique({rows.Key});
    T = table();
    for k = 1:numel(keys)
        for w = unique([rows.Warm])
            r = rows(strcmp({rows.Key}, keys{k}) & [rows.Warm] == w);
            if isempty(r), continue; end
            j100 = [r.J100];  j150 = [r.J150];
            dJ = j150 - j100;  dU = [r.Junc150] - [r.Junc100];  dO = [r.Jocc150] - [r.Jocc100];
            pW = NaN;
            if numel(r) >= 2 && any(dJ ~= 0), pW = signrank(j100, j150); end
            shareUnc = sum(dU) / sum(dJ);
            oct = NaN;
            if ~isempty(Lold) && w
                tok = sscanf(keys{k}, '%dC_CF%d_TT%d_GM%d_sp%f');
                m = [Lold.NumCameras] == tok(1) & [Lold.TargetType] == tok(3) & [Lold.GridMode] == tok(4);
                if any(m), oct = median([Lold(m).BestCost]); end
            end
            T = [T; table(string(keys{k}), w, numel(r), median(j100), median(j150), ...
                100*median(-dJ ./ j100), pW, median(dU), median(dO), shareUnc, oct, ...
                'VariableNames', {'Condition', 'Warm', 'N', 'MedianJ_100', 'MedianJ_150', ...
                'MedianGainPct', 'p_Wilcoxon', 'Median_dJunc', 'Median_dJocc', ...
                'ShareOfGain_Unc', 'Oct_Warm_Median_100'})]; %#ok<AGROW>
        end
    end

    outDir = fullfile(root, 'Results', 'Analysis');
    if ~isfolder(outDir), mkdir(outDir); end
    f = fullfile(outDir, sprintf('compareGenerations_%s_%s.txt', o.Preset, string(datetime('now'), 'yyyyMMdd')));
    fid = fopen(f, 'w');
    txt = evalc('disp(T)');
    fprintf(fid, 'Generation %s vs final, %s, CF%d, runs since %s\n', 'snapshot', o.Preset, o.CostFunc, string(o.Since));
    fprintf(fid, 'ShareOfGain_Unc: fraction of the total CF3 drop from the resolution term (rest = occlusion).\n\n%s', txt);
    fclose(fid);
    disp(T);
    fprintf('Written to %s\n', f);
end
