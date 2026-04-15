function plotGARuns(varargin)
    % Parse inputs
    p = inputParser;
    addParameter(p, 'LogFile', 'GA_RunsLog.mat', @ischar);
    addParameter(p, 'NumCameras', [], @isnumeric);
    addParameter(p, 'CostFunction', [], @isnumeric);
    addParameter(p, 'WarmStart', [], @islogical);
    addParameter(p, 'TargetType', [], @isnumeric);
    addParameter(p, 'GridMode', [], @isnumeric);
    addParameter(p, 'Spacing', [], @isnumeric);
    addParameter(p, 'MatchParams', false, @islogical);
    parse(p, varargin{:});
    
    logFile = p.Results.LogFile;
    filterCams = p.Results.NumCameras;
    filterCostFunc = p.Results.CostFunction;
    filterWarmStart = p.Results.WarmStart;
    filterTargetType = p.Results.TargetType;
    filterGridMode = p.Results.GridMode;
    filterSpacing = p.Results.Spacing;
    matchParams = p.Results.MatchParams;
    
    % Load log file
    if ~isfile(logFile)
        error('Log file not found: %s', logFile);
    end
    
    load(logFile, 'runLog');
    totalRuns = length(runLog);
    fprintf('Loaded %d total GA runs from log...\n', totalRuns);
    
    % --- Backwards compatibility: ensure new fields exist ---
    requiredFields = {'TargetType', 'GridMode', 'Spacing', 'NumTargetPoints'};
    for f = 1:length(requiredFields)
        if ~isfield(runLog, requiredFields{f})
            [runLog.(requiredFields{f})] = deal(NaN);
        end
    end
    
    % Apply filters
    keepIdx = true(totalRuns, 1);
    
    if ~isempty(filterCams)
        numCameras = [runLog.NumCameras];
        keepIdx = keepIdx & (numCameras == filterCams)';
        fprintf('Filtering: NumCameras = %d\n', filterCams);
    end
    
    if ~isempty(filterCostFunc)
        costFuncTypes = [runLog.CostFunctionType];
        keepIdx = keepIdx & (costFuncTypes == filterCostFunc)';
        fprintf('Filtering: CostFunction = %d\n', filterCostFunc);
    end
    
    if ~isempty(filterWarmStart)
        warmStarts = [runLog.WarmStart];
        keepIdx = keepIdx & (warmStarts == filterWarmStart)';
        fprintf('Filtering: WarmStart = %d\n', filterWarmStart);
    end
    
    if ~isempty(filterTargetType)
        targetTypes = [runLog.TargetType];
        keepIdx = keepIdx & (targetTypes == filterTargetType)';
        fprintf('Filtering: TargetType = %d\n', filterTargetType);
    end
    
    if ~isempty(filterGridMode)
        gridModes = [runLog.GridMode];
        keepIdx = keepIdx & (gridModes == filterGridMode)';
        fprintf('Filtering: GridMode = %d\n', filterGridMode);
    end
    
    if ~isempty(filterSpacing)
        spacings = [runLog.Spacing];
        keepIdx = keepIdx & (abs(spacings - filterSpacing) < 0.001)';
        fprintf('Filtering: Spacing = %.3f\n', filterSpacing);
    end
    
    % Filter by matching GA parameters
    if matchParams && sum(keepIdx) > 1
        refIdx = find(keepIdx, 1);
        refParams = runLog(refIdx).GAParams;
        
        for i = 1:totalRuns
            if keepIdx(i)
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
    targetTypes = [runLog.TargetType];
    gridModes = [runLog.GridMode];
    spacings = [runLog.Spacing];
    
    % Label maps
    costFuncNames = {'Res. Uncertainty', 'Dynamic Occlusion', 'Combined'};
    targetTypeNames = {'UAV', 'UGV'};
    gridModeNames = {'Uniform', 'Normal'};
    
    % Create figure with 3x3 subplots
    figure('Name', 'GA Runs Analysis', 'Position', [50, 50, 1600, 1000]);
    
    %% Subplot 1: Best Cost vs Number of Cameras (colored by cost function)

    subplot(3,3,1);
    hold on;
    
    uniqueCostFunc = unique(costFuncTypes);
    colors = lines(3);
    legendEntries = {};
    
    for cf = uniqueCostFunc
        idx = costFuncTypes == cf;
        scatter(numCameras(idx), bestCosts(idx), 100, colors(cf,:), ...
            'filled', 'MarkerFaceAlpha', 0.6);
        legendEntries{end+1} = costFuncNames{cf};
    end
    
    xlabel('Number of Cameras');
    ylabel('Best Cost');
    title('Best Cost vs Cameras (by Cost Function)');
    grid on;
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'best');
    end
    
    %% Subplot 2: Best Cost - Warm-start vs Cold-start
    subplot(3,3,2);
    
    if any(warmStarts) && any(~warmStarts)
        warmIdx = warmStarts == true;
        coldIdx = warmStarts == false;
        
        hold on;
        scatter(numCameras(warmIdx), bestCosts(warmIdx), 100, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(numCameras(coldIdx), bestCosts(coldIdx), 100, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        
        xlabel('Number of Cameras');
        ylabel('Best Cost');
        title('Warm-Start vs Cold-Start');
        legend('Warm-Start', 'Cold-Start', 'Location', 'best');
        grid on;
    else
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
    subplot(3,3,3);
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
    title('Computation Time vs Cameras');
    grid on;
    
    %% Subplot 4: Cost distribution by camera count
    subplot(3,3,4);
    
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
    
    %% Subplot 5: Cost by Target Type (UAV vs UGV)
    subplot(3,3,5);
    
    validTT = ~isnan(targetTypes);
    if any(validTT)
        uniqueTT = unique(targetTypes(validTT));
        
        if length(uniqueTT) > 1
            hold on;
            ttColors = [0.2 0.6 0.9; 0.9 0.5 0.2]; % blue=UAV, orange=UGV
            ttLegend = {};
            for t = 1:length(uniqueTT)
                tt = uniqueTT(t);
                idx = targetTypes == tt;
                scatter(numCameras(idx), bestCosts(idx), 100, ttColors(t,:), ...
                    'filled', 'MarkerFaceAlpha', 0.7);
                if tt >= 1 && tt <= 2
                    ttLegend{end+1} = targetTypeNames{tt};
                else
                    ttLegend{end+1} = sprintf('Type %d', tt);
                end
            end
            xlabel('Number of Cameras');
            ylabel('Best Cost');
            title('Cost by Target Type');
            legend(ttLegend, 'Location', 'best');
            grid on;
        else
            % Only one target type present
            scatter(numCameras(validTT), bestCosts(validTT), 100, 'filled', 'MarkerFaceAlpha', 0.6);
            xlabel('Number of Cameras');
            ylabel('Best Cost');
            if uniqueTT >= 1 && uniqueTT <= 2
                title(sprintf('Cost (%s Only)', targetTypeNames{uniqueTT}));
            else
                title('Cost by Target Type');
            end
            grid on;
        end
    else
        text(0.5, 0.5, 'No target type data', 'HorizontalAlignment', 'center');
        axis off;
        title('Cost by Target Type');
    end
    
    %% Subplot 6: Cost by Grid Mode (Uniform vs Normal)
    subplot(3,3,6);
    
    validGM = ~isnan(gridModes);
    if any(validGM)
        uniqueGM = unique(gridModes(validGM));
        
        if length(uniqueGM) > 1
            hold on;
            gmColors = [0.3 0.7 0.3; 0.7 0.3 0.7]; % green=Uniform, purple=Normal
            gmLegend = {};
            for g = 1:length(uniqueGM)
                gm = uniqueGM(g);
                idx = gridModes == gm;
                scatter(numCameras(idx), bestCosts(idx), 100, gmColors(g,:), ...
                    'filled', 'MarkerFaceAlpha', 0.7);
                if gm >= 1 && gm <= 2
                    gmLegend{end+1} = gridModeNames{gm};
                else
                    gmLegend{end+1} = sprintf('Mode %d', gm);
                end
            end
            xlabel('Number of Cameras');
            ylabel('Best Cost');
            title('Cost by Grid Mode');
            legend(gmLegend, 'Location', 'best');
            grid on;
        else
            scatter(numCameras(validGM), bestCosts(validGM), 100, 'filled', 'MarkerFaceAlpha', 0.6);
            xlabel('Number of Cameras');
            ylabel('Best Cost');
            if uniqueGM >= 1 && uniqueGM <= 2
                title(sprintf('Cost (%s Grid Only)', gridModeNames{uniqueGM}));
            else
                title('Cost by Grid Mode');
            end
            grid on;
        end
    else
        text(0.5, 0.5, 'No grid mode data', 'HorizontalAlignment', 'center');
        axis off;
        title('Cost by Grid Mode');
    end
    
    %%  Subplot 7: Grid Spacing Convergence
    subplot(3,3,7);
    
    validSP = ~isnan(spacings);
    if any(validSP)
        uniqueSP = unique(spacings(validSP));
        
        if length(uniqueSP) > 1
            % Group by camera count and plot cost vs spacing
            hold on;
            camColors = lines(length(uniqueCams));
            camLegend = {};
            
            for c = 1:length(uniqueCams)
                cam = uniqueCams(c);
                camMask = numCameras == cam & validSP;
                
                if any(camMask)
                    spVals = spacings(camMask);
                    costVals = bestCosts(camMask);
                    
                    % Compute median cost per spacing level
                    usp = unique(spVals);
                    medCosts = zeros(size(usp));
                    for s = 1:length(usp)
                        medCosts(s) = median(costVals(spVals == usp(s)));
                    end
                    
                    plot(usp, medCosts, '-o', 'Color', camColors(c,:), ...
                        'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', camColors(c,:));
                    
                    % Also scatter individual runs faintly
                    scatter(spVals, costVals, 30, camColors(c,:), ...
                        'filled', 'MarkerFaceAlpha', 0.3);
                    
                    camLegend{end+1} = sprintf('%dC', cam);
                end
            end
            
            xlabel('Grid Spacing (m)');
            ylabel('Best Cost');
            title('Cost vs Grid Spacing');
            legend(camLegend, 'Location', 'best');
            grid on;
        else
            % Only one spacing
            text(0.5, 0.5, sprintf('Single spacing: %.2fm', uniqueSP), ...
                'HorizontalAlignment', 'center');
            axis off;
            title('Grid Spacing Convergence');
        end
    else
        text(0.5, 0.5, 'No spacing data', 'HorizontalAlignment', 'center');
        axis off;
        title('Grid Spacing Convergence');
    end
    
    %% Subplot 8: Cost over Time 
    subplot(3,3,8);
    [sortedTimes, sortIdx] = sort(timestamps);
    plot(sortedTimes, bestCosts(sortIdx), '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('Run Timestamp');
    ylabel('Best Cost');
    title('Best Cost Over Time (Chronological)');
    grid on;
    xtickangle(45);
    
    %% Subplot 9: Statistics Table
    subplot(3,3,9);
    axis off;
    
    statsText = sprintf('GA Runs Statistics\n');
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
                
                statsText = [statsText sprintf('%s:\n', costFuncNames{cf})];
                statsText = [statsText sprintf('Cost: %.2f+/-%.2f (min:%.2f)\n', ...
                    meanCost, stdCost, minCost)];
            end
        end
        
        % Warm-start vs cold-start comparison
        if any(warmStarts(camIdx)) && any(~warmStarts(camIdx))
            warmIdx = camIdx & warmStarts;
            coldIdx = camIdx & ~warmStarts;
            statsText = [statsText sprintf('Warm: %.2f  Cold: %.2f\n', ...
                mean(bestCosts(warmIdx)), mean(bestCosts(coldIdx)))];
        end
        
        % Target type breakdown (if data exists)
        validTTcam = camIdx & ~isnan(targetTypes);
        if any(validTTcam)
            for tt = unique(targetTypes(validTTcam))
                ttIdx = camIdx & (targetTypes == tt);
                if sum(ttIdx) > 0 && tt >= 1 && tt <= 2
                    statsText = [statsText sprintf('  %s: %.2f+/-%.2f\n', ...
                        targetTypeNames{tt}, mean(bestCosts(ttIdx)), std(bestCosts(ttIdx)))];
                end
            end
        end
        
        statsText = [statsText sprintf('\n')];
    end
    
    text(0.05, 0.95, statsText, 'FontSize', 8, 'FontName', 'FixedWidth', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none');
    
    %% Print summary to console
    fprintf('\n');
    fprintf('%s\n', statsText);
    
    %% Save figure
    filterStr = '';
    if ~isempty(filterCams), filterStr = [filterStr sprintf('_%dCams', filterCams)]; end
    if ~isempty(filterCostFunc), filterStr = [filterStr sprintf('_CF%d', filterCostFunc)]; end
    if ~isempty(filterWarmStart), filterStr = [filterStr sprintf('_WS%d', filterWarmStart)]; end
    if ~isempty(filterTargetType), filterStr = [filterStr sprintf('_TT%d', filterTargetType)]; end
    if ~isempty(filterGridMode), filterStr = [filterStr sprintf('_GM%d', filterGridMode)]; end
    if ~isempty(filterSpacing), filterStr = [filterStr sprintf('_sp%.0fcm', filterSpacing*100)]; end
    
    figFilename = sprintf('GA_RunsComparison%s.png', filterStr);
    saveas(gcf, figFilename);
    fprintf('Comparison plot saved to: %s\n', figFilename);
end