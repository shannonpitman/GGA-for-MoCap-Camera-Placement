function summary = reEvaluateSavedRuns(varargin)
% REEVALUATESAVEDRUNS  Archive existing Results/ and recompute every saved
% chromosome's cost + coverage with the corrected (post-bugfix) code.
%
% PURPOSE
%   The lens-range swap, uv-indexing and isInsidePlanes bugs corrupted the
%   BestCost values stored in every Results/<N>Cams/*.mat file. The
%   chromosomes themselves are still valid camera placements; only the
%   recorded costs/coverage are wrong. This script:
%       1. Copies Results/ to Results_pre-bugfix/ (skip if already present)
%       2. Loads every per-run .mat file
%       3. Re-evaluates each saved BestSolution.Chromosome with the
%          corrected resUncertaintyCost / dynamicOcclusionCost /
%          combinedCostFunction
%       4. Recomputes coverage stats with the corrected findVisibleCameras
%       5. Writes corrected_<...> fields back into the same .mat alongside
%          the originals (so nothing is destroyed)
%       6. Writes a master CSV-style summary table sorted by corrected cost
%
% USAGE
%   reEvaluateSavedRuns()                      % default: archive + re-eval
%   reEvaluateSavedRuns('Archive', false)      % skip the copy step
%   reEvaluateSavedRuns('CamCounts', [7 8])    % only 7C and 8C runs
%   reEvaluateSavedRuns('DryRun', true)        % print what would be done
%
% Returns a struct array `summary` with one entry per re-evaluated run.

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'Archive',   true,           @islogical);
    addParameter(p, 'CamCounts', [6 7 8],        @isnumeric);
    addParameter(p, 'DryRun',    false,          @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    resultsDir  = fullfile(projectRoot, 'Results');
    archiveDir  = fullfile(projectRoot, 'Results_pre-bugfix');

    if ~isfolder(resultsDir)
        error('reEvaluateSavedRuns:NoResults', ...
            'Results/ folder not found at %s', resultsDir);
    end

    %% Step 1 — archive existing Results
    if opts.Archive
        if isfolder(archiveDir)
            fprintf('Archive already exists at %s — skipping copy.\n', archiveDir);
        else
            fprintf('Archiving Results/ -> Results_pre-bugfix/ ...\n');
            if ~opts.DryRun
                copyfile(resultsDir, archiveDir);
            end
            fprintf('Done.\n');
        end
    end

    %% Step 2 — gather per-run .mat files
    runFiles = {};
    for n = opts.CamCounts(:)'
        runDir = fullfile(resultsDir, sprintf('%dCams', n));
        if ~isfolder(runDir), continue; end
        d = dir(fullfile(runDir, sprintf('%dCams_Run_*.mat', n)));
        for k = 1:numel(d)
            runFiles{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
    end
    fprintf('Found %d run files to re-evaluate.\n', numel(runFiles));

    summary = repmat(struct( ...
        'File',            '', ...
        'NumCameras',      NaN, ...
        'OldBestCost',     NaN, ...
        'NewResCost',      NaN, ...
        'NewOccCost',      NaN, ...
        'NewCombCost',     NaN, ...
        'OldZeroPct',      NaN, ...
        'NewZeroPct',      NaN, ...
        'NewMinCov',       NaN, ...
        'NewAvgCov',       NaN, ...
        'CostFuncType',    NaN), 1, numel(runFiles));

    for k = 1:numel(runFiles)
        runFile = runFiles{k};
        fprintf('\n[%d/%d] %s\n', k, numel(runFiles), runFile);

        try
            tmp = load(runFile);
        catch ME
            warning('Could not load %s: %s', runFile, ME.message);
            continue;
        end
        if ~isfield(tmp, 'saveData')
            warning('Skipping (no saveData field): %s', runFile);
            continue;
        end
        sd     = tmp.saveData;
        specs  = sd.Specifications;
        chrom  = sd.BestSolution.Chromosome;

        % Backfill any fields older saves are missing (e.g. FocalWide,
        % PreComputed.maxCameraRangeWide) so the corrected cost functions
        % can evaluate.
        try
            specs = backfillLegacySpecs(specs);
        catch ME
            warning('Backfill failed on %s: %s', runFile, ME.message);
            continue;
        end

        %% Step 3 — recompute costs with the corrected code
        if opts.DryRun
            newRes = NaN; newOcc = NaN; newComb = NaN;
        else
            try
                newRes  = resUncertaintyCost(chrom, specs);
                newOcc  = dynamicOcclusionCost(chrom, specs);
                newComb = combinedCostFunction(chrom, specs);
            catch ME
                warning('Re-evaluation failed on %s: %s', runFile, ME.message);
                continue;
            end
        end

        %% Step 4 — recompute coverage stats with corrected findVisibleCameras
        [cameras, camCenters] = setupCameras(chrom, specs.Cams, ...
            specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint);
        T = specs.Target;
        nP = size(T,1);
        cov = zeros(nP,1);
        for q = 1:nP
            [vis, ~] = findVisibleCameras(T(q,:), cameras, camCenters, ...
                specs.Cams, specs.Resolution, ...
                specs.PreComputed.maxCameraRange, ...
                specs.PreComputed.maxCameraRangeWide, specs.FocalWide);
            cov(q) = numel(vis);
        end
        newCovStats.numPoints           = nP;
        newCovStats.zeroCameras         = sum(cov == 0);
        newCovStats.oneCamera           = sum(cov == 1);
        newCovStats.twoPlusCameras      = sum(cov >= 2);
        newCovStats.zeroCamerasPercent  = 100*newCovStats.zeroCameras/nP;
        newCovStats.oneCameraPercent    = 100*newCovStats.oneCamera/nP;
        newCovStats.twoPlusCamerasPercent = 100*newCovStats.twoPlusCameras/nP;
        newCovStats.avgCoverage         = mean(cov);
        newCovStats.maxCoverage         = max(cov);
        newCovStats.minCoverage         = min(cov);
        newCovStats.medianCoverage      = median(cov);

        %% Step 5 — write corrected fields back (alongside the originals)
        sd.Corrected.ResolutionCost        = newRes;
        sd.Corrected.OcclusionCost         = newOcc;
        sd.Corrected.CombinedCost          = newComb;
        sd.Corrected.CoverageStats         = newCovStats;
        sd.Corrected.RecomputedTimestamp   = datetime('now');
        sd.Corrected.Notes                 = ['Recomputed after fixing ' ...
            'lens-range swap (findVisibleCameras), uv indexing ' ...
            '(computePointUncertainty), and all->any (isInsidePlanes).'];

        if ~opts.DryRun
            saveData = sd; %#ok<NASGU>
            save(runFile, 'saveData');
        end

        %% Build summary row
        oldZeroPct = NaN;
        if isfield(sd, 'CoverageStats') && isfield(sd.CoverageStats, 'zeroCamerasPercent')
            oldZeroPct = sd.CoverageStats.zeroCamerasPercent;
        end
        summary(k).File         = runFile;
        summary(k).NumCameras   = specs.Cams;
        summary(k).OldBestCost  = sd.BestCost;
        summary(k).NewResCost   = newRes;
        summary(k).NewOccCost   = newOcc;
        summary(k).NewCombCost  = newComb;
        summary(k).OldZeroPct   = oldZeroPct;
        summary(k).NewZeroPct   = newCovStats.zeroCamerasPercent;
        summary(k).NewMinCov    = newCovStats.minCoverage;
        summary(k).NewAvgCov    = newCovStats.avgCoverage;
        summary(k).CostFuncType = guessCostType(sd);

        fprintf('  Old BestCost: %.6f | New CF3: %.6f | New min cov: %d | new 0-cov: %.1f%%\n', ...
            sd.BestCost, newComb, newCovStats.minCoverage, newCovStats.zeroCamerasPercent);
    end

    %% Step 6 — write a sortable summary CSV
    summary = summary(arrayfun(@(s) ~isempty(s.File), summary));
    if isempty(summary)
        fprintf('\nNo runs re-evaluated.\n');
        return;
    end

    [~, order] = sort([summary.NewCombCost]);
    summary = summary(order);

    csvPath = fullfile(resultsDir, 'Logs', 'corrected_costs_summary.csv');
    if ~isfolder(fileparts(csvPath)), mkdir(fileparts(csvPath)); end
    if ~opts.DryRun
        fid = fopen(csvPath, 'w');
        fprintf(fid, 'File,NumCameras,CostFuncType,OldBestCost,NewResCost,NewOccCost,NewCombCost,OldZeroPct,NewZeroPct,NewMinCov,NewAvgCov\n');
        for s = summary
            [~,bn,ext] = fileparts(s.File);
            fprintf(fid, '%s%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%d,%.3f\n', ...
                bn, ext, s.NumCameras, s.CostFuncType, s.OldBestCost, ...
                s.NewResCost, s.NewOccCost, s.NewCombCost, ...
                s.OldZeroPct, s.NewZeroPct, s.NewMinCov, s.NewAvgCov);
        end
        fclose(fid);
        fprintf('\nSummary written to: %s\n', csvPath);
    end
end


function ct = guessCostType(sd)
% Recover the cost-function type from saved data, falling back to NaN if
% it isn't recorded (older saves don't store it directly inside saveData).
    ct = NaN;
    if isfield(sd, 'GAParams') && isfield(sd.GAParams, 'CostFunctionType')
        ct = sd.GAParams.CostFunctionType;
    end
end
