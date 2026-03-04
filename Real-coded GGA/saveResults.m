function saveResults(out, specs, params, elapsedTime, costFunctionType, warmStartUsed, coverageStats)
% Outputs three files: data saved in a .mat, a .txt summary, and NB data
% collated for later comparison in a .mat log file

currentDateTime = datetime('now');
    dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
    
    %% Build save struct
    saveData.BestSolution = out.bestsol;
    saveData.BestCost = out.bestsol.Cost;
    saveData.CameraConfiguration = reshape(out.bestsol.Chromosome, 6, specs.Cams)';
    saveData.Specifications = specs;
    saveData.GAParams = params;
    saveData.ConvergenceHistory = out.bestcost;
    saveData.AvgCostHistory = out.avgcost;
    saveData.TopTenAvgCostHistory = out.topTenAvgCost;
    saveData.ElapsedTime = elapsedTime;
    saveData.Timestamp = currentDateTime;
    saveData.CoverageStats = coverageStats;

    % Per-camera breakdown
    for i = 1:specs.Cams
        chromStart = (i-1)*6 + 1;
        chromEnd = i*6;
        saveData.Cameras(i).Position = out.bestsol.Chromosome(chromStart:chromStart+2);
        saveData.Cameras(i).Orientation = out.bestsol.Chromosome(chromEnd-2:chromEnd);
        saveData.Cameras(i).OrientationDegrees = rad2deg(saveData.Cameras(i).Orientation);
    end

    %% MAT file
    matFilename = sprintf('%dCams_Run_%s.mat', specs.Cams, dateTimeStr);
    save(matFilename, 'saveData');
    fprintf('Results saved to: %s\n', matFilename);

    %% Text summary
    txtFilename = sprintf('%dCams_Run_%s.txt', specs.Cams, dateTimeStr);
    fid = fopen(txtFilename, 'w');

    fprintf(fid, 'Genetic Algorithm Camera Placement Results\n');
    fprintf(fid, '==========================================\n\n');
    fprintf(fid, 'Timestamp: %s\n', char(currentDateTime));
    fprintf(fid, 'Number of Cameras: %d\n', specs.Cams);
    fprintf(fid, 'Best Cost (Cost Function %d): %.6f\n', costFunctionType, saveData.BestCost);
    fprintf(fid, 'Cost Function Weights: %.2f Resolution Uncertainty and %.2f Dynamic Occlusion\n', ...
        specs.WeightUncertainty, specs.WeightOcclusion);
    fprintf(fid, 'Computation Time: %.2f min\n', elapsedTime / 60);
    fprintf(fid, 'Mutation Rate: %.2f\n', params.mu);
    fprintf(fid, 'Tournament Size: %d\n', params.Tournamentsize);
    fprintf(fid, 'Warm Starting Used: %s\n', string(warmStartUsed));
    fprintf(fid, 'Population Size: %d\n', params.nPop);
    fprintf(fid, 'Number of Generations: %d\n', params.MaxIt);

    fprintf(fid, '\nCamera Configurations:\n');
    fprintf(fid, '---------------------\n');
    for i = 1:specs.Cams
        fprintf(fid, '\nCamera %d:\n', i);
        fprintf(fid, '  Position (m):        [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Position);
        fprintf(fid, '  Orientation (rad):    [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Orientation);
        fprintf(fid, '  Orientation (deg):    [%.1f, %.1f, %.1f]\n', saveData.Cameras(i).OrientationDegrees);
    end

    fprintf(fid, '\nCamera Coverage Statistics:\n');
    fprintf(fid, '---------------------------\n');
    fprintf(fid, 'Total target points: %d\n', coverageStats.numPoints);
    fprintf(fid, 'Points with 0 cameras: %d (%.1f%%)\n', coverageStats.zeroCameras, coverageStats.zeroCamerasPercent);
    fprintf(fid, 'Points with 1 camera:  %d (%.1f%%)\n', coverageStats.oneCamera, coverageStats.oneCameraPercent);
    fprintf(fid, 'Points with 2+ cameras: %d (%.1f%%)\n', coverageStats.twoPlusCameras, coverageStats.twoPlusCamerasPercent);
    fprintf(fid, 'Average coverage: %.2f cameras per point\n', coverageStats.avgCoverage);
    fprintf(fid, 'Max coverage: %d | Min: %d | Median: %.2f\n', ...
        coverageStats.maxCoverage, coverageStats.minCoverage, coverageStats.medianCoverage);

    fclose(fid);
    fprintf('Summary saved to: %s\n', txtFilename);

    %% Append to master log
    masterLogFile = 'GGA_RunsLog.mat';

    newLogEntry.Timestamp = currentDateTime;
    newLogEntry.NumCameras = specs.Cams;
    newLogEntry.BestCost = saveData.BestCost;
    newLogEntry.ElapsedTime = elapsedTime;
    newLogEntry.CostFunctionType = costFunctionType;
    newLogEntry.WarmStart = warmStartUsed;
    newLogEntry.GAParams = params;
    newLogEntry.RunFilename = matFilename;

    if isfile(masterLogFile)
        load(masterLogFile, 'runLog');
        runLog(end+1) = newLogEntry; 
    else
        runLog = newLogEntry;
    end

    save(masterLogFile, 'runLog');
    fprintf('Run logged to: %s (Total runs: %d)\n', masterLogFile, length(runLog));
end