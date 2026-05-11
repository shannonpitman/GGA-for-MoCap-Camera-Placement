function [patchedPath, summary] = patchBatchLogForUGVRerun(batchLogPath, varargin)
%PATCHBATCHLOGFORUGVRERUN  Prep a saved BatchLog for a UGV re-run.
%
%   [patchedPath, summary] = patchBatchLogForUGVRerun(batchLogPath, ...)
%
%   Loads a saved BatchLog .mat (the per-batch schedule file written under
%   Results/Logs/BatchLog_<ts>.mat), resets every UGV (TargetType==2) row
%   back to 'pending', clears stale per-row state, and optionally updates
%   each UGV row's Spacing field to a new x-y value. UAV rows are left
%   untouched so a subsequent batchRunGA('ResumeLog', patchedPath) will
%   skip them.
%
%   The patched schedule is saved to a NEW file by default
%   (BatchLog_<ts>_UGVrerun_<newts>.mat) so the original schedule is
%   preserved. Pass 'InPlace', true to overwrite the original instead.
%
%   PARAMETERS:
%       NewUGVSpacing   scalar new x-y spacing (m). If empty (default),
%                       UGV rows keep their original Spacing. Set this to
%                       the x-y value chosen from spacingSensitivity_UGV.
%       InPlace         if true, overwrite the source file. Default false.
%       DryRun          if true, print the plan but don't write. Default false.
%
%   USAGE:
%       % Inspect what would change:
%       [~, s] = patchBatchLogForUGVRerun('Results/Logs/BatchLog_20260510_xxx.mat', ...
%                                         'NewUGVSpacing', 0.5, 'DryRun', true);
%
%       % Commit (writes a new file alongside the original):
%       patchedPath = patchBatchLogForUGVRerun( ...
%           'Results/Logs/BatchLog_20260510_xxx.mat', 'NewUGVSpacing', 0.5);
%
%       % Resume from the patched schedule:
%       batchRunGA('ResumeLog', patchedPath, 'UGV_ZSpacing', 0.25);
%
%   See also: archiveUGVRuns, batchRunGA.

    %% Parse
    p = inputParser;
    addRequired(p, 'batchLogPath', @ischar);
    addParameter(p, 'NewUGVSpacing', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'InPlace',       false, @islogical);
    addParameter(p, 'DryRun',        false, @islogical);
    parse(p, batchLogPath, varargin{:});
    args = p.Results;

    if ~isfile(batchLogPath)
        error('patchBatchLogForUGVRerun:NoFile', 'Batch log not found: %s', batchLogPath);
    end

    %% Load
    loaded = load(batchLogPath, 'schedule', 'batchConfig', 'batchTimestamp');
    schedule = loaded.schedule;

    isUGV = arrayfun(@(s) isequal(s.TargetType, 2), schedule);
    ugvIdx = find(isUGV);

    summary.TotalRows     = numel(schedule);
    summary.UGVRows       = numel(ugvIdx);
    summary.UGVDoneBefore = sum(strcmp({schedule(isUGV).Status}, 'done'));
    summary.UGVPending    = sum(strcmp({schedule(isUGV).Status}, 'pending'));
    summary.UGVOther      = numel(ugvIdx) - summary.UGVDoneBefore - summary.UGVPending;
    summary.NewSpacing    = args.NewUGVSpacing;

    fprintf('\n  patchBatchLogForUGVRerun  --------------------------------\n');
    fprintf('  Source        : %s\n', batchLogPath);
    fprintf('  Total rows    : %d\n', summary.TotalRows);
    fprintf('  UGV rows      : %d  (done before: %d | pending: %d | other: %d)\n', ...
        summary.UGVRows, summary.UGVDoneBefore, summary.UGVPending, summary.UGVOther);
    if isempty(args.NewUGVSpacing)
        oldSpacings = unique([schedule(isUGV).Spacing]);
        fprintf('  Spacing       : keep existing UGV Spacing(s) = %s m\n', ...
            mat2str(oldSpacings));
    else
        fprintf('  Spacing       : rewriting UGV rows to %.3g m x-y\n', args.NewUGVSpacing);
    end
    fprintf('  Dry run       : %s\n', mat2str(args.DryRun));
    fprintf('  In-place      : %s\n\n', mat2str(args.InPlace));

    if isempty(ugvIdx)
        fprintf('  No UGV rows in schedule - nothing to do.\n\n');
        patchedPath = batchLogPath;
        return;
    end

    %% Patch the UGV rows
    for k = 1:numel(ugvIdx)
        i = ugvIdx(k);
        schedule(i).Status            = 'pending';
        schedule(i).BestCost          = NaN;
        schedule(i).ElapsedTime       = NaN;
        if isfield(schedule, 'BestChromosome')
            schedule(i).BestChromosome = [];
        end
        if isfield(schedule, 'RunFilename')
            schedule(i).RunFilename   = '';
        end
        if isfield(schedule, 'ErrorMsg')
            schedule(i).ErrorMsg      = '';
        end
        % Clear the cached warm seed so it gets re-derived from the new
        % cold runs at runtime.
        if isfield(schedule, 'WarmSeedChromosome')
            schedule(i).WarmSeedChromosome = [];
        end
        if ~isempty(args.NewUGVSpacing)
            schedule(i).Spacing = args.NewUGVSpacing;
        end
    end

    %% Decide output path
    if args.InPlace
        patchedPath = batchLogPath;
    else
        [pth, base, ext] = fileparts(batchLogPath);
        newTs = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        patchedPath = fullfile(pth, sprintf('%s_UGVrerun_%s%s', base, newTs, ext));
    end

    %% Preview / write
    fprintf('  Patched UGV rows: %d  (all reset to pending', numel(ugvIdx));
    if ~isempty(args.NewUGVSpacing)
        fprintf(', Spacing -> %.3g m', args.NewUGVSpacing);
    end
    fprintf(')\n');
    fprintf('  Output        : %s\n', patchedPath);

    if args.DryRun
        fprintf('  DRY RUN - no file written.\n\n');
        return;
    end

    % Save with the same variable names batchRunGA expects.
    batchConfig    = loaded.batchConfig; %#ok<NASGU>
    batchTimestamp = loaded.batchTimestamp; %#ok<NASGU>
    save(patchedPath, 'schedule', 'batchConfig', 'batchTimestamp');
    fprintf('  Wrote patched batch log.\n');
    fprintf('  Resume with:  batchRunGA(''ResumeLog'', ''%s'', ''UGV_ZSpacing'', 0.25);\n\n', ...
        patchedPath);
end
