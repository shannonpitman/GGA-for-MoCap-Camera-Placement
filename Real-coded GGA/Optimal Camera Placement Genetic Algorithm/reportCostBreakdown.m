function T = reportCostBreakdown(varargin)

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'NumCameras', 7,   @isnumeric);
    addParameter(p, 'Spacing',    1.0, @isnumeric);
    addParameter(p, 'LogFile',    '',  @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    projectRoot = fileparts(mfilename('fullpath'));
    if isempty(opts.LogFile)
        opts.LogFile = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat');
    end

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    rows = cell(0, 11);

    fprintf('\n%-5s %-8s | %-9s %-9s %-9s | %-9s %-9s %-9s | %-8s %-8s %-8s\n', ...
        'TT', 'GM', 'GA Junc', 'GA Jocc', 'GA J', 'AH Junc', 'AH Jocc', 'AH J', ...
        'R_unc', 'R_occ', 'R_total');
    fprintf('%s\n', repmat('-', 1, 108));

    for tt = 1:2
        for gm = 1:2
            %% --- Find the best GA run for this condition -------------
            try
                [log, ~] = loadGARuns('LogFile', opts.LogFile, ...
                    'NumCameras', opts.NumCameras, 'CostFunction', 3, ...
                    'TargetType', tt, 'GridMode', gm, 'Spacing', opts.Spacing);
            catch ME
                fprintf('  [%s / %s] no matching runs (%s)\n', ...
                    ttNames{tt}, gmNames{gm}, ME.message);
                continue;
            end

            [~, bestIdx] = min([log.BestCost]);
            bestEntry = log(bestIdx);

            runFile = resolveRunPath(bestEntry.RunFilename, bestEntry.NumCameras);
            if ~isfile(runFile)
                fprintf('  [%s / %s] could not resolve run file %s\n', ...
                    ttNames{tt}, gmNames{gm}, bestEntry.RunFilename);
                continue;
            end

            tmp = load(runFile, 'saveData');
            sd    = tmp.saveData;
            specs = backfillLegacySpecs(sd.Specifications);
            chrom = sd.BestSolution.Chromosome;

            [cameras, camCenters] = setupCameras(chrom, specs.Cams, ...
                specs.Resolution, specs.Focal, specs.FocalWide, ...
                specs.PrincipalPoint, specs.PixelSize);

            rawUnc = resUncertainty(specs, cameras, camCenters);
            rawOcc = dynamicOcclusion(specs, cameras, camCenters);

            [gaTotal, gaJunc, gaJocc] = cf3Terms(rawUnc, rawOcc, specs);

            adhoc = evaluateOptiTrackCost('TargetType', tt, 'GridMode', gm, ...
                'Spacing', opts.Spacing, ...
                'WeightUnc', specs.WeightUncertainty, ...
                'WeightOcc', specs.WeightOcclusion);

            ratioUnc   = adhoc.WeightedUnc / gaJunc;
            ratioOcc   = adhoc.WeightedOcc / gaJocc;
            ratioTotal = adhoc.CF3 / gaTotal;

            fprintf('%-5s %-8s | %-9.4f %-9.4f %-9.4f | %-9.4f %-9.4f %-9.4f | %6.1fx  %6.1fx  %6.1fx\n', ...
                ttNames{tt}, gmNames{gm}, gaJunc, gaJocc, gaTotal, ...
                adhoc.WeightedUnc, adhoc.WeightedOcc, adhoc.CF3, ...
                ratioUnc, ratioOcc, ratioTotal);

            % Sanity check against the logged / recomputed combined cost
            if abs(gaTotal - bestEntry.BestCost) > 1e-3
                fprintf('    NOTE: recomputed J (%.4f) differs from logged BestCost (%.4f)\n', ...
                    gaTotal, bestEntry.BestCost);
                fprintf('    (expected if this run predates a cost-function bugfix — see reEvaluateSavedRuns.m)\n');
            end

            rows(end+1, :) = {ttNames{tt}, gmNames{gm}, gaJunc, gaJocc, gaTotal, ...
                adhoc.WeightedUnc, adhoc.WeightedOcc, adhoc.CF3, ...
                ratioUnc, ratioOcc, ratioTotal}; %#ok<AGROW>
        end
    end

    fprintf('\n');

    if isempty(rows)
        T = table();
        return;
    end

    T = cell2table(rows, 'VariableNames', ...
        {'TargetType', 'GridMode', 'GA_Junc', 'GA_Jocc', 'GA_Total', ...
         'Adhoc_Junc', 'Adhoc_Jocc', 'Adhoc_Total', ...
         'Ratio_Unc', 'Ratio_Occ', 'Ratio_Total'});
end
