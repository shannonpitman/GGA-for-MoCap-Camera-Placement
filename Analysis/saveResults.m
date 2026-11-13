function [matFilename, matFullPath] = saveResults(out, specs, params, elapsedTime, costFunctionType, warmStartUsed, coverageStats)
% Outputs three files: data saved in a .mat, a .txt summary, and NB data
% collated for later comparison in a .mat log file. Also adds the run's best
% solution and the top cfg.ArchiveFraction of its final population to the
% warm-start seed archive (Results/SeedArchive/<preset>/<condition>.mat).
%
% Run files go to Results/<N>Cams/ for the lab preset and to
% Results/<preset>/<N>Cams/ for any other preset (see runResultsDir).

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
    if isfield(out, 'popDiversity')
        % Per-locus genotype diversity (std across population, averaged
        % over loci) — recorded so plotGA_PopulationDiversity can use a
        % real diversity metric instead of the cost-spread proxy.
        saveData.PopulationStdHistory = out.popDiversity;
    end
    saveData.ElapsedTime = elapsedTime;
    saveData.Timestamp = currentDateTime;
    saveData.CoverageStats = coverageStats;

    [preset, hardware] = presetOf(specs);
    saveData.Preset   = preset;
    saveData.Hardware = hardware;

    % What this run started from (warm runs) ...
    saveData.WarmSeeds = specs.warmChromosomes;
    if isfield(specs, 'WarmSeedInfo')
        saveData.WarmSeedInfo = specs.WarmSeedInfo;
    end
    saveData.NumSeedsInjected = size(specs.warmChromosomes, 1) * logical(warmStartUsed);

    % ... and the best part of where it ended (later warm starts draw on it).
    archiveFrac = 0.10;
    if isfield(specs, 'RunConfig') && isfield(specs.RunConfig, 'ArchiveFraction')
        archiveFrac = specs.RunConfig.ArchiveFraction;
    end
    if isfield(out, 'pop') && ~isempty(out.pop)
        nTop = max(1, min(numel(out.pop), round(archiveFrac * numel(out.pop))));
        saveData.FinalTopK.Chromosomes = vertcat(out.pop(1:nTop).Chromosome);
        saveData.FinalTopK.Costs       = [out.pop(1:nTop).Cost].';
    end

    % Per-camera breakdown
    for i = 1:specs.Cams
        chromStart = (i-1)*6 + 1;
        chromEnd = i*6;
        saveData.Cameras(i).Position = out.bestsol.Chromosome(chromStart:chromStart+2);
        saveData.Cameras(i).Orientation = out.bestsol.Chromosome(chromEnd-2:chromEnd);  % rotation vector [rad]
        saveData.Cameras(i).OrientationDegrees = rad2deg(rotm2eul( ...
            genesToRotm(saveData.Cameras(i).Orientation), "XYZ"));                     % XYZ Euler [deg], display
    end

    %% MAT file (written into Results/<N>Cams/, RunFilename stays a basename
    %  so existing log entries that store basenames keep resolving.)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));   % .../<proj>
    runDir      = runResultsDir(preset, specs.Cams);
    if ~isfolder(runDir), mkdir(runDir); end

    matFilename = sprintf('%dCams_Run_%s.mat', specs.Cams, dateTimeStr);
    matFullPath = fullfile(runDir, matFilename);
    save(matFullPath, 'saveData');
    fprintf('Results saved to: %s\n', matFullPath);

    %% Text summary
    txtFilename = sprintf('%dCams_Run_%s.txt', specs.Cams, dateTimeStr);
    txtFullPath = fullfile(runDir, txtFilename);
    fid = fopen(txtFullPath, 'w');

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

    % --- NEW: Log target modality and grid info ---
    targetTypeNames = {'UAV (Full Volume)', 'UGV (Floor Slab)'};
    gridModeNames = {'Uniform Grid', 'Normal Grid'};
    if isfield(specs, 'TargetType') && specs.TargetType >= 1 && specs.TargetType <= 2
        fprintf(fid, 'Target Type: %s\n', targetTypeNames{specs.TargetType});
    end
    if isfield(specs, 'TargetMode') && specs.TargetMode >= 1 && specs.TargetMode <= 2
        fprintf(fid, 'Grid Mode: %s\n', gridModeNames{specs.TargetMode});
    end
    fprintf(fid, 'Number of Target Points: %d\n', specs.NumPoints);

    fprintf(fid, '\nCamera Configurations:\n');
    fprintf(fid, '---------------------\n');
    for i = 1:specs.Cams
        fprintf(fid, '\nCamera %d:\n', i);
        fprintf(fid, '  Position (m):        [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Position);
        fprintf(fid, '  Rotation vector (rad): [%.3f, %.3f, %.3f]\n', saveData.Cameras(i).Orientation);
        fprintf(fid, '  XYZ Euler (deg):       [%.1f, %.1f, %.1f]\n', saveData.Cameras(i).OrientationDegrees);
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
    fprintf('Summary saved to: %s\n', txtFullPath);

    %% Append to master log (Results/Logs/GGA_RunsLog.mat)
    logsDir = fullfile(projectRoot, 'Results', 'Logs');
    if ~isfolder(logsDir), mkdir(logsDir); end
    masterLogFile = fullfile(logsDir, 'GGA_RunsLog.mat');

    newLogEntry.Timestamp = currentDateTime;
    newLogEntry.NumCameras = specs.Cams;
    newLogEntry.BestCost = saveData.BestCost;
    newLogEntry.ElapsedTime = elapsedTime;
    newLogEntry.CostFunctionType = costFunctionType;
    newLogEntry.WarmStart = warmStartUsed;
    newLogEntry.GAParams = params;
    newLogEntry.RunFilename = matFilename;
    
    if isfield(specs, 'TargetType')
        newLogEntry.TargetType = specs.TargetType;
    else
        newLogEntry.TargetType = NaN; % unknown for legacy runs
    end

    if isfield(specs, 'TargetMode')
        newLogEntry.GridMode = specs.TargetMode;
    else
        newLogEntry.GridMode = NaN;
    end

    newLogEntry.Spacing = specs.spacing;
    newLogEntry.NumTargetPoints = specs.NumPoints;
    newLogEntry.Preset   = preset;
    newLogEntry.Hardware = hardware;

    if isfile(masterLogFile)
        load(masterLogFile, 'runLog');
        
        % Backwards compatibility: ensure existing entries have new fields ---
        requiredFields = {'TargetType', 'GridMode', 'Spacing', 'NumTargetPoints'};
        for f = 1:length(requiredFields)
            if ~isfield(runLog, requiredFields{f})
                [runLog.(requiredFields{f})] = deal(NaN);
            end
        end
        % Runs logged before presets existed were all the OptiTrack lab rig
        if ~isfield(runLog, 'Preset'),   [runLog.Preset]   = deal('optitrack_lab'); end
        if ~isfield(runLog, 'Hardware'), [runLog.Hardware] = deal('optitrack');     end
        
        runLog(end+1) = newLogEntry; 
    else
        runLog = newLogEntry;
    end

    save(masterLogFile, 'runLog');
    fprintf('Run logged to: %s (Total runs: %d)\n', masterLogFile, length(runLog));

    %% Warm-start seed archive
    if ~isfield(specs, 'ArchiveRun') || specs.ArchiveRun
        try
            condKey = conditionKey(specs.Cams, costFunctionType, specs.TargetType, ...
                                   specs.TargetMode, specs.spacing);
            if isfield(saveData, 'FinalTopK')
                C = saveData.FinalTopK.Chromosomes;  J = saveData.FinalTopK.Costs;
            else
                C = out.bestsol.Chromosome;  J = out.bestsol.Cost;
            end
            src = repmat({'finalpop'}, size(C, 1), 1);
            src{1} = 'best';
            % make sure row 1 really is the run best
            C(1, :) = out.bestsol.Chromosome;  J(1) = out.bestsol.Cost;
            n = updateSeedArchive(preset, condKey, C, J, src, matFilename, warmStartUsed);
            fprintf('Seed archive %s/%s: +%d solutions\n', preset, condKey, n);
        catch ME
            warning('saveResults:archive', 'Seed archive not updated: %s', ME.message);
        end
    end
end

function [preset, hardware] = presetOf(specs)
    preset = 'optitrack_lab';
    if isfield(specs, 'RunConfig') && isfield(specs.RunConfig, 'Preset')
        preset = specs.RunConfig.Preset;
    end
    hardware = 'optitrack';
    if isfield(specs, 'HardwareProfile')
        hardware = specs.HardwareProfile;
    end
end
