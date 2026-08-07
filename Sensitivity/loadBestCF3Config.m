function [chromosome, bestRun] = loadBestCF3Config(numCams, targetType, gridMode)
%LOADBESTCF3CONFIG  Load the best logged CF3 chromosome for a given case.
%
%   [chromosome, bestRun] = loadBestCF3Config(numCams) returns the chromosome
%   for the lowest-cost CF3 run with the given camera count, defaulting to
%   targetType=1 (UAV) and gridMode=1 (uniform).
%
%   loadBestCF3Config(numCams, targetType, gridMode) restricts the search
%   to runs matching those metadata fields.
%
%   This mirrors the filtering logic in analyseConfiguration.m so the
%   sensitivity scripts and the analyser pull the same "current best" run.

    if nargin < 2 || isempty(targetType), targetType = 1; end
    if nargin < 3 || isempty(gridMode),   gridMode   = 1; end

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    logFile = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat');
    if ~isfile(logFile)
        error('loadBestCF3Config:LogMissing', ...
            'Master log not found at %s. Run a CF3 batch first.', logFile);
    end

    S = load(logFile, 'runLog');
    runLog = S.runLog;

    % Backwards compatibility for legacy entries.
    if ~isfield(runLog, 'TargetType'), [runLog.TargetType] = deal(NaN); end
    if ~isfield(runLog, 'GridMode'),   [runLog.GridMode]   = deal(NaN); end

    mask = ([runLog.NumCameras]      == numCams) & ...
           ([runLog.CostFunctionType] == 3)      & ...
           ([runLog.TargetType]       == targetType) & ...
           ([runLog.GridMode]         == gridMode);
    matches = runLog(mask);
    if isempty(matches)
        error('loadBestCF3Config:NoMatch', ...
            'No CF3 entries for %d cams, TargetType=%d, GridMode=%d.', ...
            numCams, targetType, gridMode);
    end

    [~, bestIdx] = min([matches.BestCost]);
    bestRun = matches(bestIdx);

    runPath = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
    if ~isfile(runPath)
        error('loadBestCF3Config:FileMissing', ...
            'Run file not resolvable: %s', runPath);
    end
    sd = load(runPath, 'saveData');
    chromosome = sd.saveData.BestSolution.Chromosome(:)';
end
