function runFullRerun(varargin)
%RUNFULLRERUN  The complete re-run in one call: CF1/CF2, normalisation, CF3.
%
%   runFullRerun()                      % 320 runs: CF1/CF2 x5, CF3 x10
%   runFullRerun('CF12Repeats', 10)     % 480 runs: everything x10
%   runFullRerun('DryRun', true)        % print both schedules, run nothing
%
%   Phase 1  batchRunGA with CostFunctions [1 2] (resolution-only and
%            occlusion-only runs; cold + warm, CF12Repeats each).
%   Phase 2  buildNormFromBatch on the phase-1 batch log: rebuilds the CF3
%            utopia/nadir table (Results/normTable.mat for the lab preset)
%            from the new CF1/CF2 results. The old table is archived.
%   Phase 3  batchRunGA with CostFunctions 3 (combined cost; cold + warm,
%            CF3Repeats each), which reads the new table.
%
%   All other settings (cameras [6 7], UAV+UGV, uniform+normal grids,
%   spacing, mutation, mount model, hardware) come from runConfig(Preset).
%
%   Resuming after an interruption
%     Phase 1 or 3 stopped part-way:
%       runFullRerun('StartPhase', 1, 'ResumeLog', 'Results/Logs/BatchLog_<date>.mat')
%       runFullRerun('StartPhase', 3, 'ResumeLog', 'Results/Logs/BatchLog_<date>.mat')
%     Phase 1 finished but phase 2/3 did not start:
%       runFullRerun('StartPhase', 2, 'CF12Log', 'Results/Logs/BatchLog_<date>.mat')

    root = addProjectPaths();

    p = inputParser;
    addParameter(p, 'Preset',      'optitrack_lab', @ischar);
    addParameter(p, 'CF12Repeats', 5,     @isnumeric);
    addParameter(p, 'CF3Repeats',  10,    @isnumeric);
    addParameter(p, 'StartPhase',  1,     @(x) any(x == [1 2 3]));
    addParameter(p, 'ResumeLog',   '',    @ischar);
    addParameter(p, 'CF12Log',     '',    @ischar);
    addParameter(p, 'DryRun',      false, @islogical);
    parse(p, varargin{:});
    o = p.Results;

    t0 = tic;
    cf12Log = o.CF12Log;

    %% Phase 1: resolution-only and occlusion-only runs
    if o.StartPhase <= 1
        banner('PHASE 1: CF1 + CF2', o.DryRun);
        args = {'Preset', o.Preset, 'CostFunctions', [1 2], 'NumRepeats', o.CF12Repeats, 'DryRun', o.DryRun};
        if ~isempty(o.ResumeLog), args = [args, {'ResumeLog', o.ResumeLog}]; end
        before = newestLog(root);
        batchRunGA(args{:});
        if ~o.DryRun
            cf12Log = newestLog(root);
            if ~isempty(o.ResumeLog), cf12Log = o.ResumeLog; end
            assert(~isempty(cf12Log) && ~strcmp(cf12Log, before) || ~isempty(o.ResumeLog), ...
                'runFullRerun:NoLog', 'Phase 1 did not write a batch log.');
        end
    end

    %% Phase 2: rebuild utopia/nadir normalisation from the new CF1/CF2 runs
    if o.StartPhase <= 2 && ~o.DryRun
        banner('PHASE 2: normalisation table', false);
        if isempty(cf12Log)
            error('runFullRerun:NoCF12Log', 'Pass ''CF12Log'' (the phase-1 batch log) to start at phase 2.');
        end
        S = load(cf12Log, 'schedule', 'batchConfig');
        nDone = sum(strcmp({S.schedule.Status}, 'done'));
        fprintf('  Using %s (%d of %d runs done)\n', cf12Log, nDone, numel(S.schedule));
        buildNormFromBatch(S.schedule, S.batchConfig);
    end

    %% Phase 3: combined-cost runs
    banner('PHASE 3: CF3', o.DryRun);
    args = {'Preset', o.Preset, 'CostFunctions', 3, 'NumRepeats', o.CF3Repeats, 'DryRun', o.DryRun};
    if o.StartPhase == 3 && ~isempty(o.ResumeLog), args = [args, {'ResumeLog', o.ResumeLog}]; end
    batchRunGA(args{:});

    fprintf('\n  runFullRerun finished in %.1f h\n', toc(t0)/3600);
end

function f = newestLog(root)
    d = dir(fullfile(root, 'Results', 'Logs', 'BatchLog_*.mat'));
    if isempty(d)
        f = '';
    else
        [~, i] = max([d.datenum]);
        f = fullfile(d(i).folder, d(i).name);
    end
end

function banner(txt, dry)
    if dry, txt = [txt ' (dry run)']; end
    fprintf('\n==================================================\n  %s\n==================================================\n', txt);
end
