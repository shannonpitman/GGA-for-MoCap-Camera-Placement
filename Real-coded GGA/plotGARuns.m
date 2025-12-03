function plotGARuns(varargin)
% Plot comparison of multiple GA runs from log file with filtering options
% 
% Usage:
%   plotGARuns()                          % Plot all runs
%   plotGARuns('NumCameras', 7)           % Only 7-camera runs
%   plotGARuns('CostFunction', 1)         % Only resolution uncertainty runs
%   plotGARuns('WarmStart', true)         % Only warm-start runs
%   plotGARuns('NumCameras', 7, 'WarmStart', false)  % Combined filters
%
% Available filters:
%   'NumCameras'     - Filter by number of cameras (scalar)
%   'CostFunction'   - Filter by cost function type (1, 2, or 3)
%   'WarmStart'      - Filter by warm-start usage (true/false)
%   'MatchParams'    - Only compare runs with identical GA parameters (true/false)
%   'LogFile'        - Specify log file (default: 'GA_RunsLog.mat')

    % Parse inputs
    p = inputParser;
    addParameter(p, 'LogFile', 'GA_RunsLog.mat', @ischar);
    addParameter(p, 'NumCameras', [], @isnumeric);
    addParameter(p, 'CostFunction', [], @isnumeric);
    addParameter(p, 'WarmStart', [], @islogical);
    addParameter(p, 'MatchParams', false, @islogical);
    parse(p, varargin{:});
    
    logFile = p.Results.LogFile;
    filterCams = p.Results.NumCameras;
    filterCostFunc = p.Results.CostFunction;
    filterWarmStart = p.Results.WarmStart;
    matchParams = p.Results.MatchParams;
    
    % Load log file
    if ~isfile(logFile)
        error('Log file not found: %s', logFile);
    end
    
    load(logFile, 'runLog');
    totalRuns = length(runLog);
    fprintf('Loaded %d total GA runs from log...\n', totalRuns);
    
    % Apply filters
    keepIdx = true(totalRuns, 1);
    
    if ~isempty(filterCams)
        numCameras = [runLog.NumCameras];
        keepIdx = keepIdx & (numCameras == filterCams);
        fprintf('Filtering: NumCameras = %d\n', filterCams);
    end
    
    if ~isempty(filterCostFunc)
        costFuncTypes = [runLog.CostFunctionType];
        keepIdx = keepIdx & (costFuncTypes == filterCostFunc);
        fprintf('Filtering: CostFunction = %d\n', filterCostFunc);
    end
    
    if ~isempty(filterWarmStart)
        warmStarts = [runLog.WarmStart];
        keepIdx = keepIdx & (warmStarts == filterWarmStart);
        fprintf('Filtering: WarmStart = %d\n', filterWarmStart);
    end
    
    % Filter by matching GA parameters
    if matchParams && sum(keepIdx) > 1
        % Get first kept run's parameters as reference
        refIdx = find(keepIdx, 1);
        refParams = runLog(refIdx).GAParams;
        
        for i = 1:totalRuns
            if keepIdx(i)
                % Compare all GA parameters
                if ~isequal(runLog(i).GAParams.MaxIt, refParams.MaxIt) || ...
                   ~isequal(runLog(i).GAParams.nPop, refParams.nPop) || ...
                   ~isequal(runLog(i).GAParams.beta, refParams.beta) || ...
                   ~isequal(runLog(i).GAParams.mu, refParams.mu)
                    keepIdx(i) = false;
                end
            end
        end
        fprintf('Filtering: Matching GA parameters only\n');
    end
    
    % Apply filter
    runLog = runLog(keepIdx);
    numRuns = length(runLog);
    
    if numRuns == 0
        error('No runs match the specified filters!');
    end
    
    fprintf('Plotting %d runs after filtering...\n\n', numRuns);
    
    % Extract data
    numCameras = [runLog.NumCameras];
    bestCosts = [runLog.BestCost];
    elapsedTimes = [runLog.ElapsedTime];
    timestamps = [runLog.Timestamp];
    costFuncTypes = [runLog.CostFunctionType];
    warmStarts = [runLog.WarmStart];
    
    % Cost function names for labels
    costFuncNames = {'Res. Uncertainty', 'Dynamic Occlusion', 'Combined'};
    
    % Create figure with subplots
    figure('Name', 'GA Runs Analysis', 'Position', [100, 100, 1400, 900]);
    
    %% Subplot 1: Best Cost vs Number of Cameras (colored by cost function)
    subplot(2,3,1);
    hold on;
    
    uniqueCostFunc = unique(costFuncTypes);
    colors = lines(3); % Get 3 distinct colors
    legendEntries = {};
    
    for cf = uniqueCostFunc
        idx = costFuncTypes == cf;
        scatter(numCameras(idx), bestCosts(idx), 100, colors(cf,:), ...
            'filled', 'MarkerFaceAlpha', 0.6);
        legendEntries{end+1} = costFuncNames{cf};
    end
    
    xlabel('Number of Cameras');
    ylabel('Best Cost');
    title('Best Cost vs Number of Cameras');
    grid on;
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'best');
    end
    
    %% Subplot 2: Best Cost comparison (Warm-start vs Cold-start)
    subplot(2,3,2);
    
    if any(warmStarts) && any(~warmStarts)
        % Both warm and cold starts exist
        warmIdx = warmStarts == true;
        coldIdx = warmStarts == false;
        
        hold on;
        scatter(numCameras(warmIdx), bestCosts(warmIdx), 100, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(numCameras(coldIdx), bestCosts(coldIdx), 100, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        
        xlabel('Number of Cameras');
        ylabel('Best Cost');
        title('Warm-Start vs Cold-Start Comparison');
        legend('Warm-Start', 'Cold-Start', 'Location', 'best');
        grid on;
    else
        % Only one type
        scatter(numCameras, bestCosts, 100, 'filled', 'MarkerFaceAlpha', 0.6);
        xlabel('Number of Cameras');
        ylabel('Best Cost');
        if all(warmStarts)
            title('Best Cost (Warm-Start Only)');
        else
            title('Best Cost (Cold-Start Only)');
        end
        grid on;
    end
    
    %% Subplot 3: Computation Time vs Number of Cameras
    subplot(2,3,3);
    scatter(numCameras, elapsedTimes, 100, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    uniqueCams = unique(numCameras);
    for cam = uniqueCams
        idx = numCameras == cam;
        meanTime = mean(elapsedTimes(idx));
        plot(cam, meanTime, 'rx', 'MarkerSize', 15, 'LineWidth', 2);
    end
    
    xlabel('Number of Cameras');
    ylabel('Computation Time (seconds)');
    title('Computation Time vs Number of Cameras');
    grid on;
    
    %% Subplot 4: Cost distribution by camera count
    subplot(2,3,4);
    
    if length(uniqueCams) > 1
        boxData = [];
        boxLabels = {};
        for cam = uniqueCams
            idx = numCameras == cam;
            boxData = [boxData; bestCosts(idx)'];
            boxLabels = [boxLabels; repmat({sprintf('%dC', cam)}, sum(idx), 1)];
        end
        boxplot(boxData, boxLabels);
        ylabel('Best Cost');
        xlabel('Number of Cameras');
        title('Cost Distribution by Camera Count');
        grid on;
    else
        histogram(bestCosts, 10, 'FaceColor', [0.3 0.6 0.9]);
        xlabel('Best Cost');
        ylabel('Frequency');
        title(sprintf('Cost Distribution (%d Cameras)', uniqueCams));
        grid on;
    end
    
    %% Subplot 5: Cost over Time (chronological)
    subplot(2,3,5);
    [sortedTimes, sortIdx] = sort(timestamps);
    plot(sortedTimes, bestCosts(sortIdx), '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('Run Timestamp');
    ylabel('Best Cost');
    title('Best Cost Over Time (Chronological)');
    grid on;
    xtickangle(45);
    
    %% Subplot 6: Statistics Table
    subplot(2,3,6);
    axis off;
    
    % Create statistics text
    statsText = sprintf('GA Runs Statistics\n');
    statsText = [statsText sprintf('==================\n\n')];
    statsText = [statsText sprintf('Total Filtered Runs: %d\n\n', numRuns)];
    
    % Group by camera count
    for cam = uniqueCams
        camIdx = numCameras == cam;
        numRuns_cam = sum(camIdx);
        
        statsText = [statsText sprintf('%d Cameras (%d runs):\n', cam, numRuns_cam)];
        
        % Stats by cost function type
        for cf = unique(costFuncTypes(camIdx))
            idx = camIdx & (costFuncTypes == cf);
            if sum(idx) > 0
                meanCost = mean(bestCosts(idx));
                stdCost = std(bestCosts(idx));
                minCost = min(bestCosts(idx));
                
                statsText = [statsText sprintf('  %s:\n', costFuncNames{cf})];
                statsText = [statsText sprintf('    Cost: %.2f±%.2f (min:%.2f)\n', ...
                    meanCost, stdCost, minCost)];
            end
        end
        
        % Warm-start vs cold-start comparison
        if any(warmStarts(camIdx)) && any(~warmStarts(camIdx))
            warmIdx = camIdx & warmStarts;
            coldIdx = camIdx & ~warmStarts;
            statsText = [statsText sprintf('  Warm: %.2f  Cold: %.2f\n', ...
                mean(bestCosts(warmIdx)), mean(bestCosts(coldIdx)))];
        end
        
        statsText = [statsText sprintf('\n')];
    end
    
    text(0.05, 0.95, statsText, 'FontSize', 9, 'FontName', 'FixedWidth', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none');
    
    %% Print summary to console
    fprintf('\n');
    fprintf('%s\n', statsText);
    
    % Save figure with descriptive name
    filterStr = '';
    if ~isempty(filterCams), filterStr = [filterStr sprintf('_%dCams', filterCams)]; end
    if ~isempty(filterCostFunc), filterStr = [filterStr sprintf('_CF%d', filterCostFunc)]; end
    if ~isempty(filterWarmStart), filterStr = [filterStr sprintf('_WS%d', filterWarmStart)]; end
    
    figFilename = sprintf('GA_RunsComparison%s.png', filterStr);
    saveas(gcf, figFilename);
    fprintf('Comparison plot saved to: %s\n', figFilename);
end