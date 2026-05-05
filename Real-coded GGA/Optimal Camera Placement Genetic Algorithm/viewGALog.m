function viewGALog(varargin)
% Display summary of all logged GA runs with optional filtering
% 
% Usage:
%   viewGALog()                          % View all runs
%   viewGALog('NumCameras', 7)           % Only 7-camera runs
%   viewGALog('WarmStart', true)         % Only warm-start runs

    % Parse inputs
    p = inputParser;
    addParameter(p, 'LogFile', 'GA_RunsLog.mat', @ischar);
    addParameter(p, 'NumCameras', [], @isnumeric);
    addParameter(p, 'CostFunction', [], @isnumeric);
    addParameter(p, 'WarmStart', [], @islogical);
    parse(p, varargin{:});
    
    logFile = p.Results.LogFile;
    
    if ~isfile(logFile)
        error('Log file not found: %s', logFile);
    end
    
    load(logFile, 'runLog');
    totalRuns = length(runLog);
    
    % Apply filters
    keepIdx = true(totalRuns, 1);
    
    if ~isempty(p.Results.NumCameras)
        numCameras = [runLog.NumCameras];
        keepIdx = keepIdx & (numCameras == p.Results.NumCameras);
    end
    
    if ~isempty(p.Results.CostFunction)
        costFuncTypes = [runLog.CostFunctionType];
        keepIdx = keepIdx & (costFuncTypes == p.Results.CostFunction);
    end
    
    if ~isempty(p.Results.WarmStart)
        warmStarts = [runLog.WarmStart];
        keepIdx = keepIdx & (warmStarts == p.Results.WarmStart);
    end
    
    runLog = runLog(keepIdx);
    
    fprintf('\n=== GA Runs Log (%d of %d runs shown) ===\n\n', length(runLog), totalRuns);
    fprintf('%-4s %-20s %-8s %-10s %-12s %-10s %-6s\n', ...
        'Run', 'Timestamp', 'Cameras', 'CostFunc', 'Best Cost', 'Time (hr)', 'Warm?');
    fprintf('%s\n', repmat('-', 1, 80));
    
    costFuncNames = {'ResUnc', 'DynOcc', 'Combined'};
    
    for i = 1:length(runLog)
        warmStr = 'No';
        if runLog(i).WarmStart
            warmStr = 'Yes';
        end
        
        fprintf('%-4d %-20s %-8d %-10s %-12.4f %-10.1f %-6s\n', ...
            i, ...
            char(runLog(i).Timestamp), ...
            runLog(i).NumCameras, ...
            costFuncNames{runLog(i).CostFunctionType}, ...
            runLog(i).BestCost, ...
            runLog(i).ElapsedTime/3600, ...
            warmStr);
    end
    
    fprintf('\nTotal runs shown: %d\n', length(runLog));
end