function S = thesisResultsSummary(varargin)
%THESISRESULTSSUMMARY  Every number the write-up quotes, in one pass.
%
%   S = thesisResultsSummary()
%   S = thesisResultsSummary('NumCameras', [6 7], 'Spacing', 1.0)
%
%   Recomputes, from the runs in Results/ and under the CURRENT
%   Results/normTable.mat utopia/nadir normalisation, every quantity the
%   paper and dissertation chapter report:
%
%     0. Normalisation constants actually in force (utopia, nadir, norm)
%     1. GA parameters as they were really run (not as first drafted)
%     2. CF3 best/median/IQR per (cameras x target type x grid mode),
%        split cold vs warm start
%     3. Optimised GA Rig vs Manually Posed Rig: combined J, and the weighted
%        occlusion and resolution terms separately, with ratios
%     4. Coverage: 0 / 1 / 2+ camera percentages and mean cameras per point
%     5. Pairwise baseline angles: fraction inside the [40,140] band
%     6. Convergence: initial / mid / final medians, drop, warm-vs-cold
%     7. Wall-clock time per run
%
%   Everything printed is also written to
%   figures/results_summary.txt, and the machine-readable tables to
%   figures/tables/*.csv, so the chapter can be updated by copy-paste and
%   the CSVs can be re-imported into LaTeX.
%
%   Name-Value parameters:
%     'NumCameras'  camera counts to report.        Default [6 7]
%     'Spacing'     evaluation grid spacing [m].    Default 1.0
%     'CostFunction' cost function for the cost tables. Default 3
%     'OutDir'      where to write.  Default <projectRoot>/figures
%
%   The GA-vs-ad-hoc comparison is only meaningful at 7 cameras (the real
%   rig has seven), so sections 3-5 are restricted to 7 regardless of
%   'NumCameras'.

    projectRoot = addProjectPaths();

    p = inputParser;
    addParameter(p, 'NumCameras',   [6 7], @isnumeric);
    addParameter(p, 'Spacing',      1.0,   @isnumeric);
    addParameter(p, 'CostFunction', 3,     @isnumeric);
    addParameter(p, 'OutDir',       fullfile(projectRoot, 'figures'), @ischar);
    parse(p, varargin{:});
    o = p.Results;

    tableDir = fullfile(o.OutDir, 'tables');
    if ~isfolder(o.OutDir),  mkdir(o.OutDir);  end
    if ~isfolder(tableDir),  mkdir(tableDir);  end

    logFile = fullfile(projectRoot, 'Results', 'Logs', 'GGA_RunsLog.mat');

    txtPath = fullfile(o.OutDir, 'results_summary.txt');
    if isfile(txtPath), delete(txtPath); end
    diary(txtPath); diary on;
    cleanup = onCleanup(@() diary('off'));

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    banner('RESULTS SUMMARY');
    fprintf('Generated:      %s\n', char(datetime('now')));
    fprintf('Project root:   %s\n', projectRoot);
    fprintf('Log file:       %s\n', logFile);
    fprintf('Grid spacing:   %.2f m\n', o.Spacing);
    fprintf('Camera counts:  %s\n', mat2str(o.NumCameras));

    S = struct();

    %% ==================================================================
    %  0. Normalisation constants in force
    %  ==================================================================
    banner('0. NORMALISATION CONSTANTS (Results/normTable.mat)');
    fprintf(['CF3 scales each term as (raw - utopia)/(nadir - utopia), so utopia -> 0\n' ...
             'and nadir -> 1 in each objective BEFORE the 0.5/0.5 weighting. A value\n' ...
             'above 1 in either term therefore means the configuration is worse than\n' ...
             'the nadir reference for that objective.\n\n']);

    ntFile = fullfile(projectRoot, 'Results', 'normTable.mat');
    if isfile(ntFile)
        NT = load(ntFile, 'normTable');
        nt = NT.normTable;
        keep = ismember([nt.NumCameras], o.NumCameras) & ...
               abs([nt.Spacing] - o.Spacing) < 1e-9;
        nt = nt(keep);
        Tnorm = table();
        Tnorm.TargetType = ttNames([nt.TargetType])';
        Tnorm.GridMode   = gmNames([nt.GridMode])';
        Tnorm.NumCameras = [nt.NumCameras]';
        Tnorm.utopiaUnc  = [nt.utopiaUnc]';
        Tnorm.nadirUnc   = [nt.nadirUnc]';
        Tnorm.uncertNorm = [nt.uncertNorm]';
        Tnorm.utopiaOcc  = [nt.utopiaOcc]';
        Tnorm.nadirOcc   = [nt.nadirOcc]';
        Tnorm.occlNorm   = [nt.occlNorm]';
        disp(Tnorm);
        writetable(Tnorm, fullfile(tableDir, 'norm_constants.csv'));
        S.NormConstants = Tnorm;
    else
        fprintf('  normTable.mat not found — CF3 is running on setupCostParams defaults.\n');
    end

    %% ==================================================================
    %  1. GA parameters actually used
    %  ==================================================================
    banner('1. GA PARAMETERS AS RUN');
    L = load(logFile, 'runLog');
    runLog = L.runLog;
    S.NumRunsTotal = numel(runLog);

    gp = [runLog.GAParams];
    pf = fieldnames(gp);
    names  = strings(numel(pf)+1, 1);
    values = strings(numel(pf)+1, 1);
    for i = 1:numel(pf)
        v = [gp.(pf{i})];
        names(i)  = string(pf{i});
        values(i) = strjoin(compose('%g', unique(v(:))'), ' / ');
    end
    names(end)  = "numRunsLogged";
    values(end) = string(numel(runLog));
    Tp = table(names, values, 'VariableNames', {'Parameter', 'Value'});
    disp(Tp);
    fprintf(['Note: population size auto-scales as numCams x 6 x 10, hence the two\n' ...
             'values for nPop (360 at 6 cameras, 420 at 7).\n']);
    writetable(Tp, fullfile(tableDir, 'ga_parameters.csv'));
    S.GAParams = Tp;

    %% ==================================================================
    %  2. CF3 cost distribution per cell
    %  ==================================================================
    banner('2. CF3 COST DISTRIBUTION PER CELL');
    rows = {};
    for nc = o.NumCameras(:)'
        for tt = 1:2
            for gm = 1:2
                for ws = [false true]
                    m = [runLog.NumCameras] == nc & ...
                        [runLog.CostFunctionType] == o.CostFunction & ...
                        [runLog.TargetType] == tt & ...
                        [runLog.GridMode] == gm & ...
                        [runLog.WarmStart] == ws & ...
                        abs([runLog.Spacing] - o.Spacing) < 1e-6;
                    if ~any(m), continue; end
                    c = [runLog(m).BestCost];
                    rows(end+1,:) = {nc, ttNames{tt}, gmNames{gm}, ...
                        ternary(ws,'warm','cold'), numel(c), min(c), ...
                        median(c), quantile(c,0.25), quantile(c,0.75), ...
                        mean(c), std(c)}; %#ok<AGROW>
                end
            end
        end
    end
    Tcost = cell2table(rows, 'VariableNames', {'NumCameras','TargetType', ...
        'GridMode','Start','n','Best','Median','Q25','Q75','Mean','Std'});
    disp(Tcost);
    writetable(Tcost, fullfile(tableDir, 'cf3_cost_distribution.csv'));
    S.CostDistribution = Tcost;

    % Warm vs cold effect, per cell
    banner('2b. WARM-START EFFECT (median final cost, and run-to-run spread)');
    wrows = {};
    for nc = o.NumCameras(:)'
        for tt = 1:2
            for gm = 1:2
                cold = Tcost(Tcost.NumCameras==nc & strcmp(Tcost.TargetType,ttNames{tt}) & ...
                    strcmp(Tcost.GridMode,gmNames{gm}) & strcmp(Tcost.Start,'cold'), :);
                warm = Tcost(Tcost.NumCameras==nc & strcmp(Tcost.TargetType,ttNames{tt}) & ...
                    strcmp(Tcost.GridMode,gmNames{gm}) & strcmp(Tcost.Start,'warm'), :);
                if isempty(cold) || isempty(warm), continue; end
                wrows(end+1,:) = {nc, ttNames{tt}, gmNames{gm}, ...
                    cold.Median, warm.Median, ...
                    100*(warm.Median - cold.Median)/cold.Median, ...
                    cold.Std, warm.Std, ...
                    100*(warm.Std - cold.Std)/cold.Std}; %#ok<AGROW>
            end
        end
    end
    Twarm = cell2table(wrows, 'VariableNames', {'NumCameras','TargetType','GridMode', ...
        'ColdMedian','WarmMedian','MedianChangePct','ColdStd','WarmStd','StdChangePct'});
    disp(Twarm);
    writetable(Twarm, fullfile(tableDir, 'warm_cold_effect.csv'));
    S.WarmCold = Twarm;

    %% ==================================================================
    %  3. Optimised GA Rig vs Manually Posed Rig  (paper Tables II, III, IV)
    %  ==================================================================
    banner('3. GA-BEST vs OPTITRACK AD-HOC  (7 cameras)');
    fprintf(['Junc and Jocc are the weighted terms as they enter J, i.e.\n' ...
             '   Junc = w_res*(rawUnc - utopiaUnc)/uncertNorm,  J = Junc + Jocc.\n\n']);
    Tbreak = reportCostBreakdown('NumCameras', 7, 'Spacing', o.Spacing, ...
                                 'LogFile', logFile);
    writetable(Tbreak, fullfile(tableDir, 'cost_breakdown_7cam.csv'));
    S.CostBreakdown = Tbreak;

    if ~isempty(Tbreak)
        fprintf('\n--- Paper Table II: best returned combined cost ---\n');
        fprintf('%-5s %-8s %10s %10s %8s\n', 'TT','GM','Optimised GA Rig','Ad-hoc','Ratio');
        for i = 1:height(Tbreak)
            fprintf('%-5s %-8s %10.4f %10.4f %7.1fx\n', Tbreak.TargetType{i}, ...
                Tbreak.GridMode{i}, Tbreak.GA_Total(i), Tbreak.Adhoc_Total(i), ...
                Tbreak.Ratio_Total(i));
        end
        fprintf('\n--- Paper Table III: occlusion contribution ---\n');
        fprintf('%-5s %-8s %10s %10s %8s %10s\n', 'TT','GM','GA Jocc','Ad-hoc Jocc','Ratio','GA % of J');
        for i = 1:height(Tbreak)
            fprintf('%-5s %-8s %10.4f %10.4f %7.1fx %9.1f%%\n', Tbreak.TargetType{i}, ...
                Tbreak.GridMode{i}, Tbreak.GA_Jocc(i), Tbreak.Adhoc_Jocc(i), ...
                Tbreak.Ratio_Occ(i), 100*Tbreak.GA_Jocc(i)/Tbreak.GA_Total(i));
        end
        fprintf('\n--- Paper Table IV: resolution-uncertainty contribution ---\n');
        fprintf('%-5s %-8s %10s %10s %8s %10s\n', 'TT','GM','GA Junc','Ad-hoc Junc','Ratio','GA % of J');
        for i = 1:height(Tbreak)
            fprintf('%-5s %-8s %10.4f %10.4f %7.1fx %9.1f%%\n', Tbreak.TargetType{i}, ...
                Tbreak.GridMode{i}, Tbreak.GA_Junc(i), Tbreak.Adhoc_Junc(i), ...
                Tbreak.Ratio_Unc(i), 100*Tbreak.GA_Junc(i)/Tbreak.GA_Total(i));
        end
        fprintf('\nCombined-cost improvement spans %.1fx - %.1fx across the four cells.\n', ...
            min(Tbreak.Ratio_Total), max(Tbreak.Ratio_Total));
        fprintf('Occlusion share of the GA combined cost spans %.1f%% - %.1f%%.\n', ...
            min(100*Tbreak.GA_Jocc./Tbreak.GA_Total), ...
            max(100*Tbreak.GA_Jocc./Tbreak.GA_Total));
    end

    %% ==================================================================
    %  4 & 5. Coverage and baseline angles, Optimised GA Rig vs ad-hoc
    %  ==================================================================
    banner('4/5. COVERAGE AND PAIRWISE BASELINE ANGLES (7 cameras)');
    covRows = {};
    angRows = {};
    optiChrom = buildOptiTrackChromosome();

    for tt = 1:2
        for gm = 1:2
            try
                [chrom, specs, cost] = bestRunFor(runLog, 7, o.CostFunction, tt, gm, o.Spacing);
            catch ME
                fprintf('  [%s/%s] skipped: %s\n', ttNames{tt}, gmNames{gm}, ME.message);
                continue;
            end

            [~, gaCov]   = perTargetCoverage(chrom,     specs);
            [~, optiCov] = perTargetCoverage(optiChrom, specs);

            covRows(end+1,:) = {ttNames{tt}, gmNames{gm}, specs.NumPoints, cost, ...
                gaCov.zeroPct, gaCov.onePct, gaCov.twoPlusPct, gaCov.avg, ...
                optiCov.zeroPct, optiCov.onePct, optiCov.twoPlusPct, optiCov.avg}; %#ok<AGROW>

            angGA   = pairwiseBaselineAngles(chrom,     specs);
            angOpti = pairwiseBaselineAngles(optiChrom, specs);
            lo = specs.PreComputed.minTriangAngle;
            hi = specs.PreComputed.maxTriangAngle;
            inBand = @(a) 100 * sum(a >= lo & a <= hi) / max(numel(a), 1);

            angRows(end+1,:) = {ttNames{tt}, gmNames{gm}, ...
                inBand(angGA), median(angGA), numel(angGA), ...
                inBand(angOpti), median(angOpti), numel(angOpti)}; %#ok<AGROW>

            fprintf('  %s / %s done (N = %d target points)\n', ...
                ttNames{tt}, gmNames{gm}, specs.NumPoints);
        end
    end

    Tcov = cell2table(covRows, 'VariableNames', {'TargetType','GridMode','NumPoints', ...
        'GA_Cost','GA_ZeroPct','GA_OnePct','GA_TwoPlusPct','GA_AvgCams', ...
        'Adhoc_ZeroPct','Adhoc_OnePct','Adhoc_TwoPlusPct','Adhoc_AvgCams'});
    fprintf('\n--- Coverage ---\n'); disp(Tcov);
    writetable(Tcov, fullfile(tableDir, 'coverage_7cam.csv'));
    S.Coverage = Tcov;

    Tang = cell2table(angRows, 'VariableNames', {'TargetType','GridMode', ...
        'GA_InBandPct','GA_MedianDeg','GA_NumPairs', ...
        'Adhoc_InBandPct','Adhoc_MedianDeg','Adhoc_NumPairs'});
    fprintf('\n--- Pairwise baseline angles, in-band = [%d, %d] deg ---\n', 40, 140);
    disp(Tang);
    writetable(Tang, fullfile(tableDir, 'baseline_angles_7cam.csv'));
    S.BaselineAngles = Tang;

    fprintf('Well-conditioned pairs: GA %.1f%% - %.1f%%, ad-hoc %.1f%% - %.1f%%.\n', ...
        min(Tang.GA_InBandPct), max(Tang.GA_InBandPct), ...
        min(Tang.Adhoc_InBandPct), max(Tang.Adhoc_InBandPct));
    fprintf('Single-camera coverage (GA):  %.2f%% - %.2f%% of points seen by >=1 camera.\n', ...
        min(100 - Tcov.GA_ZeroPct), max(100 - Tcov.GA_ZeroPct));
    fprintf('Unobserved points: GA %.2f%% - %.2f%%, ad-hoc %.2f%% - %.2f%%.\n', ...
        min(Tcov.GA_ZeroPct), max(Tcov.GA_ZeroPct), ...
        min(Tcov.Adhoc_ZeroPct), max(Tcov.Adhoc_ZeroPct));
    fprintf('Mean cameras per point: GA %.2f - %.2f, ad-hoc %.2f - %.2f.\n', ...
        min(Tcov.GA_AvgCams), max(Tcov.GA_AvgCams), ...
        min(Tcov.Adhoc_AvgCams), max(Tcov.Adhoc_AvgCams));

    %% ==================================================================
    %  6. Convergence behaviour
    %  ==================================================================
    banner('6. CONVERGENCE BEHAVIOUR (CF3)');
    convRows = {};
    for nc = o.NumCameras(:)'
        for tt = 1:2
            for gm = 1:2
                for ws = [false true]
                    m = find([runLog.NumCameras] == nc & ...
                             [runLog.CostFunctionType] == o.CostFunction & ...
                             [runLog.TargetType] == tt & ...
                             [runLog.GridMode] == gm & ...
                             [runLog.WarmStart] == ws & ...
                             abs([runLog.Spacing] - o.Spacing) < 1e-6);
                    if isempty(m), continue; end

                    H = [];
                    for k = m
                        try
                            f = resolveRunPath(runLog(k).RunFilename, runLog(k).NumCameras);
                            sd = load(f, 'saveData').saveData;
                            h = sd.ConvergenceHistory(:);
                            if isempty(H), H = h; else, H(:, end+1) = h; end %#ok<AGROW>
                        catch
                        end
                    end
                    if isempty(H), continue; end

                    G      = size(H, 1);
                    gMid   = min(50, G);
                    medIni = median(H(1, :));
                    medMid = median(H(gMid, :));
                    medFin = median(H(end, :));

                    convRows(end+1,:) = {nc, ttNames{tt}, gmNames{gm}, ...
                        ternary(ws,'warm','cold'), size(H,2), G, ...
                        medIni, medMid, medFin, ...
                        100*(medMid - medIni)/medIni, ...
                        100*(medFin - medIni)/medIni, ...
                        min(H(end,:)), std(H(end,:))}; %#ok<AGROW>
                end
            end
        end
    end
    Tconv = cell2table(convRows, 'VariableNames', {'NumCameras','TargetType','GridMode', ...
        'Start','n','Generations','MedianGen1','MedianGen50','MedianFinal', ...
        'DropBy50Pct','DropTotalPct','BestFinal','StdFinal'});
    disp(Tconv);
    writetable(Tconv, fullfile(tableDir, 'convergence_stats.csv'));
    S.Convergence = Tconv;

    %% ==================================================================
    %  7. Computation time
    %  ==================================================================
    banner('7. WALL-CLOCK TIME PER RUN');
    trows = {};
    for nc = o.NumCameras(:)'
        for cf = 1:3
            for tt = 1:2
                m = [runLog.NumCameras] == nc & ...
                    [runLog.CostFunctionType] == cf & ...
                    [runLog.TargetType] == tt & ...
                    abs([runLog.Spacing] - o.Spacing) < 1e-6;
                if ~any(m), continue; end
                t  = [runLog(m).ElapsedTime];
                np = [runLog(m).NumTargetPoints];
                trows(end+1,:) = {nc, cf, ttNames{tt}, numel(t), median(np), ...
                    median(t)/60, min(t)/60, max(t)/60}; %#ok<AGROW>
            end
        end
    end
    Ttime = cell2table(trows, 'VariableNames', {'NumCameras','CostFunction', ...
        'TargetType','n','TargetPoints','MedianMin','MinMin','MaxMin'});
    disp(Ttime);
    writetable(Ttime, fullfile(tableDir, 'computation_time.csv'));
    S.ComputationTime = Ttime;

    %% ==================================================================
    banner('WRITTEN');
    fprintf('  %s\n', txtPath);
    d = dir(fullfile(tableDir, '*.csv'));
    for i = 1:numel(d)
        fprintf('  %s\n', fullfile(tableDir, d(i).name));
    end

    diary off;
end


%% ---- Local helpers ------------------------------------------------------

function banner(txt)
    fprintf('\n\n%s\n', repmat('=', 1, 78));
    fprintf('  %s\n', txt);
    fprintf('%s\n\n', repmat('=', 1, 78));
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function [chrom, specs, cost] = bestRunFor(runLog, nc, cf, tt, gm, sp)
    m = [runLog.NumCameras] == nc & [runLog.CostFunctionType] == cf & ...
        [runLog.TargetType] == tt & [runLog.GridMode] == gm & ...
        abs([runLog.Spacing] - sp) < 1e-6;
    cand = runLog(m);
    if isempty(cand)
        error('no runs for %dC CF%d TT%d GM%d sp%.2f', nc, cf, tt, gm, sp);
    end
    [cost, i] = min([cand.BestCost]);
    sd = load(resolveRunPath(cand(i).RunFilename, nc), 'saveData').saveData;
    chrom = sd.BestSolution.Chromosome;
    specs = backfillLegacySpecs(sd.Specifications);
end
