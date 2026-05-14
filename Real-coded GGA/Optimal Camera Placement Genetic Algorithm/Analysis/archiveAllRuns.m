function summary = archiveAllRuns(varargin)
%ARCHIVEALLRUNS  Archive every existing GA result so a fresh batch can be
%                run from a clean slate.
%
%   summary = archiveAllRuns()                 % dry run, prints the plan
%   summary = archiveAllRuns('DryRun', false)  % actually perform the moves
%   summary = archiveAllRuns('ArchiveTag', 'pre_FOVfix')
%
%   Why: the camera FOV bug (PixelSize not passed to CentralCamera) means
%   every existing GA run was optimised against an unphysically wide FOV
%   (~99°/122°). After the fix, the cost function changes magnitude and
%   the GA-best chromosomes become apples-to-oranges versus future runs.
%   This helper moves the ENTIRE result tree out of the way under
%   Archive/<timestamped tag>/ so the next batch starts fresh.
%
%   What it moves (in order):
%     1. Results/Logs/GGA_RunsLog.mat       -> Archive/<tag>/Logs/
%     2. Results/Logs/<other>.mat backups   -> Archive/<tag>/Logs/
%     3. Results/<N>Cams/* (.mat + .txt)    -> Archive/<tag>/<N>Cams/
%     4. figures/* (.pdf etc.)              -> Archive/<tag>/figures/
%
%   The function is dry-run-by-default. Inspect the summary, then re-run
%   with 'DryRun', false to commit. After committing, the source folders
%   exist but are empty; the next batchRunGA invocation will recreate the
%   master log and per-cam-count subfolders.
%
%   summary fields:
%     ArchiveDir          target folder
%     PlannedFiles        cell array {srcPath, dstPath}
%     Performed           true if the moves were committed
%     LogBackupPath       absolute path to the moved master log
%
%   See also: archiveUGVRuns, saveResults, batchRunGA.

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'DryRun',      true,            @islogical);
    addParameter(p, 'ArchiveTag',  'pre_FOVfix',    @ischar);
    addParameter(p, 'ProjectRoot', '',              @ischar);
    parse(p, varargin{:});
    cfg = p.Results;

    if isempty(cfg.ProjectRoot)
        cfg.ProjectRoot = fileparts(fileparts(mfilename('fullpath')));
    end

    resultsDir = fullfile(cfg.ProjectRoot, 'Results');
    figsDir    = fullfile(cfg.ProjectRoot, 'figures');

    ts         = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    archiveDir = fullfile(cfg.ProjectRoot, 'Archive', ...
                          sprintf('%s_%s', ts, cfg.ArchiveTag));

    fprintf('\n  archiveAllRuns -------------------------------------------\n');
    fprintf('  Project root  : %s\n', cfg.ProjectRoot);
    fprintf('  Archive dest  : %s\n', archiveDir);
    fprintf('  Dry run       : %s\n\n', mat2str(cfg.DryRun));

    %% Build the move plan
    moves = {};

    % 1. Master log + any sibling backups in Results/Logs/
    logsSrc = fullfile(resultsDir, 'Logs');
    if isfolder(logsSrc)
        logFiles = dir(fullfile(logsSrc, '*.mat'));
        for k = 1:numel(logFiles)
            src = fullfile(logsSrc, logFiles(k).name);
            dst = fullfile(archiveDir, 'Logs', logFiles(k).name);
            moves(end+1, :) = {src, dst}; %#ok<AGROW>
        end
    end

    % 2. Per-cam-count subfolders under Results/
    if isfolder(resultsDir)
        camDirs = dir(fullfile(resultsDir, '*Cams'));
        camDirs = camDirs([camDirs.isdir]);
        for d = 1:numel(camDirs)
            srcDir = fullfile(resultsDir, camDirs(d).name);
            files  = dir(fullfile(srcDir, '*'));
            files  = files(~[files.isdir]);
            for k = 1:numel(files)
                src = fullfile(srcDir, files(k).name);
                dst = fullfile(archiveDir, camDirs(d).name, files(k).name);
                moves(end+1, :) = {src, dst}; %#ok<AGROW>
            end
        end
    end

    % 3. Generated figures
    if isfolder(figsDir)
        figFiles = dir(fullfile(figsDir, '*'));
        figFiles = figFiles(~[figFiles.isdir]);
        for k = 1:numel(figFiles)
            src = fullfile(figsDir, figFiles(k).name);
            dst = fullfile(archiveDir, 'figures', figFiles(k).name);
            moves(end+1, :) = {src, dst}; %#ok<AGROW>
        end
    end

    %% Preview
    summary.ArchiveDir    = archiveDir;
    summary.PlannedFiles  = moves;
    summary.Performed     = false;
    summary.LogBackupPath = '';

    fprintf('  Files to move : %d\n', size(moves, 1));
    if size(moves, 1) == 0
        fprintf('  Nothing to archive.\n');
        return;
    end

    % Group by destination folder for a tidy preview
    [dstFolders, ~, idx] = unique(cellfun(@fileparts, moves(:,2), 'UniformOutput', false));
    for f = 1:numel(dstFolders)
        n = sum(idx == f);
        fprintf('   %4d -> %s\n', n, relativePath(dstFolders{f}, cfg.ProjectRoot));
    end
    fprintf('\n');

    %% Execute or stop
    if cfg.DryRun
        fprintf('  DRY RUN - no files moved.\n');
        fprintf('  Re-run with archiveAllRuns(''DryRun'', false) to commit.\n\n');
        return;
    end

    %% Move files
    failures = {};
    for k = 1:size(moves, 1)
        src = moves{k, 1};
        dst = moves{k, 2};
        dstFolder = fileparts(dst);
        if ~isfolder(dstFolder), mkdir(dstFolder); end
        try
            movefile(src, dst);
        catch ME
            failures(end+1, :) = {src, ME.message}; %#ok<AGROW>
        end
        if endsWith(src, 'GGA_RunsLog.mat') && contains(src, fullfile('Results','Logs'))
            summary.LogBackupPath = dst;
        end
    end

    summary.Performed = true;

    fprintf('  Moved      : %d\n', size(moves, 1) - size(failures, 1));
    fprintf('  Failures   : %d\n', size(failures, 1));
    if ~isempty(failures)
        fprintf('  --- failed moves ---\n');
        for k = 1:size(failures, 1)
            fprintf('    %s : %s\n', failures{k,1}, failures{k,2});
        end
    end
    fprintf('  Archive at : %s\n', archiveDir);
    fprintf('  Done. Source folders are empty and ready for a fresh batchRunGA.\n\n');
end


function rel = relativePath(absPath, root)
    if startsWith(absPath, root)
        rel = absPath(numel(root)+2:end);
    else
        rel = absPath;
    end
end
