function summary = archiveUGVRuns(varargin)
%ARCHIVEUGVRUNS  Move UGV-flagged GA runs out of Results/ and prune them
%                from the master log so they aren't treated as current.
%
%   summary = archiveUGVRuns()                 % dry run, interactive prompt
%   summary = archiveUGVRuns('DryRun', false)  % actually perform the moves
%   summary = archiveUGVRuns('ArchiveName', 'Results_UGV_fine_grid_archive')
%
%   Why: the legacy UGV runs used an isotropic 0.25 m grid (33x33x3 = 3267
%   target points) because batchRunGA clamped UGV spacing to UGV_MaxHeight/2.
%   With the new anisotropic spacing (coarse x-y, fixed z), older UGV best
%   costs are not comparable - move them aside so loadBestCF3Config and
%   downstream analysis only see the new runs.
%
%   What it does (in order):
%     1. Loads Results/Logs/GGA_RunsLog.mat
%     2. Writes a timestamped backup of that .mat
%     3. Lists every entry with TargetType == 2 (UGV)
%     4. Moves the matching .mat AND sister .txt to the archive folder,
%        preserving the Results/<N>Cams/ subfolder structure
%     5. Removes those entries from runLog and saves the pruned log
%
%   The function is dry-run-by-default. Inspect the summary, then re-run
%   with 'DryRun', false to commit.
%
%   summary fields:
%       NumUGVRuns            count of UGV entries found
%       Moved                 cell array {srcPath, dstPath} of files moved
%       Missing               cell array of expected files not found on disk
%       LogBackupPath         path to the pre-prune log backup
%       PrunedLogPath         path to the pruned log (same name overwritten)
%       Performed             true if files were actually moved
%
%   See also: saveResults, loadBestCF3Config.

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'DryRun',      true,                       @islogical);
    addParameter(p, 'ArchiveName', 'Results_UGV_fine_grid_archive', @ischar);
    addParameter(p, 'ProjectRoot', '',                         @ischar);
    parse(p, varargin{:});
    cfg = p.Results;

    if isempty(cfg.ProjectRoot)
        cfg.ProjectRoot = fileparts(fileparts(mfilename('fullpath')));   % .../<proj>
    end

    resultsDir = fullfile(cfg.ProjectRoot, 'Results');
    logFile    = fullfile(resultsDir, 'Logs', 'GGA_RunsLog.mat');
    archiveDir = fullfile(cfg.ProjectRoot, cfg.ArchiveName);

    if ~isfile(logFile)
        error('archiveUGVRuns:NoLog', 'Master log not found: %s', logFile);
    end

    %% Load + identify
    S = load(logFile, 'runLog');
    runLog = S.runLog;
    if ~isfield(runLog, 'TargetType'), [runLog.TargetType] = deal(NaN); end

    ugvMask = arrayfun(@(r) isequal(r.TargetType, 2), runLog);
    ugvIdx  = find(ugvMask);

    fprintf('\n  archiveUGVRuns  ----------------------------------------\n');
    fprintf('  Project root  : %s\n', cfg.ProjectRoot);
    fprintf('  Master log    : %s (entries: %d)\n', logFile, numel(runLog));
    fprintf('  Archive dest  : %s\n', archiveDir);
    fprintf('  UGV runs found: %d\n', numel(ugvIdx));
    fprintf('  Dry run       : %s\n\n', mat2str(cfg.DryRun));

    summary.NumUGVRuns    = numel(ugvIdx);
    summary.Moved         = {};
    summary.Missing       = {};
    summary.LogBackupPath = '';
    summary.PrunedLogPath = logFile;
    summary.Performed     = false;

    if isempty(ugvIdx)
        fprintf('  Nothing to do.\n');
        return;
    end

    %% Build per-entry move plan
    moves = repmat(struct('srcMat','','dstMat','','srcTxt','','dstTxt','', ...
                          'numCams',NaN,'logIdx',NaN,'cost',NaN,'ts',''), ...
                   numel(ugvIdx), 1);

    for k = 1:numel(ugvIdx)
        idx = ugvIdx(k);
        entry = runLog(idx);
        numCams = entry.NumCameras;
        baseMat = entry.RunFilename;
        baseTxt = regexprep(baseMat, '\.mat$', '.txt');

        srcDir = fullfile(resultsDir, sprintf('%dCams', numCams));
        dstDir = fullfile(archiveDir, sprintf('%dCams', numCams));

        moves(k).srcMat  = fullfile(srcDir, baseMat);
        moves(k).dstMat  = fullfile(dstDir, baseMat);
        moves(k).srcTxt  = fullfile(srcDir, baseTxt);
        moves(k).dstTxt  = fullfile(dstDir, baseTxt);
        moves(k).numCams = numCams;
        moves(k).logIdx  = idx;
        moves(k).cost    = entry.BestCost;
        moves(k).ts      = char(string(entry.Timestamp, 'yyyy-MM-dd HH:mm:ss'));
    end

    %% Preview
    fprintf('  Plan (idx | cams | cost | timestamp | file):\n');
    for k = 1:numel(moves)
        fprintf('   %4d | %2d | %.6f | %s | %s\n', ...
            moves(k).logIdx, moves(k).numCams, moves(k).cost, moves(k).ts, ...
            relativePath(moves(k).srcMat, cfg.ProjectRoot));
    end
    fprintf('\n');

    %% Execute or stop
    if cfg.DryRun
        fprintf('  DRY RUN - no files moved, log not modified.\n');
        fprintf('  Re-run with ("DryRun", false) to commit.\n\n');
        return;
    end

    %% Back up log
    ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    backupPath = fullfile(resultsDir, 'Logs', sprintf('GGA_RunsLog_pre_UGVarchive_%s.mat', ts));
    copyfile(logFile, backupPath);
    summary.LogBackupPath = backupPath;
    fprintf('  Log backup written to: %s\n', backupPath);

    %% Move files
    for k = 1:numel(moves)
        if ~isfolder(fileparts(moves(k).dstMat))
            mkdir(fileparts(moves(k).dstMat));
        end
        moved = false;
        if isfile(moves(k).srcMat)
            movefile(moves(k).srcMat, moves(k).dstMat);
            summary.Moved{end+1} = {moves(k).srcMat, moves(k).dstMat}; %#ok<AGROW>
            moved = true;
        else
            summary.Missing{end+1} = moves(k).srcMat; %#ok<AGROW>
        end
        if isfile(moves(k).srcTxt)
            movefile(moves(k).srcTxt, moves(k).dstTxt);
            summary.Moved{end+1} = {moves(k).srcTxt, moves(k).dstTxt}; %#ok<AGROW>
        elseif ~moved
            % Only flag txt as missing if the mat was also missing -
            % otherwise the txt was probably already archived in a prior run.
            summary.Missing{end+1} = moves(k).srcTxt; %#ok<AGROW>
        end
    end

    %% Prune log
    runLog(ugvMask) = [];
    save(logFile, 'runLog');
    summary.PrunedLogPath = logFile;
    summary.Performed     = true;

    fprintf('  Files moved   : %d\n', numel(summary.Moved));
    fprintf('  Files missing : %d\n', numel(summary.Missing));
    fprintf('  Log entries pruned: %d  (remaining: %d)\n', numel(ugvIdx), numel(runLog));
    fprintf('  Done.\n\n');
end


function rel = relativePath(absPath, root)
    if startsWith(absPath, root)
        rel = absPath(numel(root)+2:end);
    else
        rel = absPath;
    end
end
