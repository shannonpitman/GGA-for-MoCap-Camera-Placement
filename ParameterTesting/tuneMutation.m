function result = tuneMutation(varargin)
%TUNEMUTATION  Tune the GA mutation rate and step sizes.
%
%   result = tuneMutation()            % full study (simulation machine)
%   result = tuneMutation('Quick', true)   % 2-minute smoke test
%
%   Three parameters are tuned together (default ranges narrowed from the
%   2026-09-30 Mac pilot, Results/Tuning/MacPilot):
%     MutationRate      mu       per-gene mutation probability   [0.02, 0.3]
%     MutationSigmaPos  s_pos    position step s.d. [m]          [0.1, 0.5]
%     MutationSigmaRot  s_rot    rotation-vector step s.d.       [1, 10] deg
%   The 1/L rule (mu = 1/numGenes = 0.024 at 7 cameras) lies inside the
%   mu range.
%
%   STAGE 1 - search. Bayesian optimisation (bayesopt, expected improvement)
%   over the three parameters. Each evaluation runs the GA on one scenario
%   at a reduced budget (Generations) for NumSeeds seeds and scores the
%   MEDIAN best combined cost J. Seeds are shared by every candidate
%   (common random numbers), so differences come from the parameters, not
%   the draw.
%
%   STAGE 2 - validation. The best ValidateTop candidates from stage 1, the
%   ExtraCandidates (default: the pilot winner) and the current runConfig
%   setting are re-run at the full budget
%   (FullGenerations) for ValidateSeeds seeds each. The recommendation is
%   the candidate with the lowest median final J; the IQR is reported so a
%   lucky median is visible.
%
%   The combined cost uses fixed normalisation (UseNormTable = false,
%   setupCostParams defaults) so tuning does not depend on a utopia/nadir
%   table that is itself rebuilt from GA runs.
%
%   Results (every evaluation, validation runs, recommendation) are saved to
%   Results/Tuning/tuneMutation_<timestamp>.mat and a .txt summary. Apply the
%   recommendation by setting MutationRate / MutationSigmaPos /
%   MutationSigmaRot in Setup/runConfig.m.
%
%   Name-value parameters
%     'Preset'          runConfig preset             'optitrack_lab'
%     'NumCameras'      7
%     'TargetType'      1 (UAV)       'GridMode' 1   'Spacing' 1.0
%     'Generations'     stage-1 budget                40
%     'NumSeeds'        stage-1 seeds per candidate    3
%     'MaxEvaluations'  bayesopt evaluations          30
%     'ValidateTop'     candidates re-run in stage 2   3
%     'ValidateSeeds'   seeds per stage-2 candidate    5
%     'FullGenerations' stage-2 budget               100
%     'RateRange'       mu search range                [0.02 0.3]
%     'SigmaPosRange'   s_pos search range [m]         [0.1 0.5]
%     'SigmaRotDegRange' s_rot search range [deg]      [1 10]
%     'ExtraCandidates' struct array (Label, MutationRate, MutationSigmaPos,
%                       SigmaRotDeg) validated in stage 2; default = pilot winner
%     'Quick'           tiny smoke-test settings     false

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'Preset',          'optitrack_lab', @ischar);
    addParameter(p, 'NumCameras',      7,    @isnumeric);
    addParameter(p, 'TargetType',      1,    @isnumeric);
    addParameter(p, 'GridMode',        1,    @isnumeric);
    addParameter(p, 'Spacing',         1.0,  @isnumeric);
    addParameter(p, 'Generations',     40,   @isnumeric);
    addParameter(p, 'NumSeeds',        3,    @isnumeric);
    addParameter(p, 'MaxEvaluations',  30,   @isnumeric);
    addParameter(p, 'ValidateTop',     3,    @isnumeric);
    addParameter(p, 'ValidateSeeds',   5,    @isnumeric);
    addParameter(p, 'FullGenerations', 100,  @isnumeric);
    addParameter(p, 'PopulationSize',  [],   @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'RateRange',        [0.02 0.3], @isnumeric);
    addParameter(p, 'SigmaPosRange',    [0.1 0.5],  @isnumeric);
    addParameter(p, 'SigmaRotDegRange', [1 10],     @isnumeric);
    addParameter(p, 'ExtraCandidates', struct('Label', 'pilot', 'MutationRate', 0.114, ...
        'MutationSigmaPos', 0.406, 'SigmaRotDeg', 4.2), @isstruct);
    addParameter(p, 'Quick',           false, @islogical);
    addParameter(p, 'OutputDir',       '',   @ischar);
    parse(p, varargin{:});
    o = p.Results;

    if o.Quick
        o.Generations = 2;  o.NumSeeds = 1;  o.MaxEvaluations = 4;
        o.ValidateTop = 1;  o.ValidateSeeds = 1;  o.FullGenerations = 2;
        o.PopulationSize = 12;  o.NumCameras = 4;  o.TargetType = 2;  o.Spacing = 2.0;
    end
    if isempty(o.OutputDir)
        o.OutputDir = fullfile(addProjectPaths(), 'Results', 'Tuning');
    end
    if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end

    base = runConfig(o.Preset, 'PopulationSize', o.PopulationSize);
    stamp = string(datetime('now'), 'yyyyMMdd_HHmmss');

    fprintf('\n  tuneMutation — %s, %d cams, TT%d GM%d sp=%.2f\n', o.Preset, ...
        o.NumCameras, o.TargetType, o.GridMode, o.Spacing);
    fprintf('  Stage 1: %d bayesopt evaluations x %d seeds x %d generations\n', ...
        o.MaxEvaluations, o.NumSeeds, o.Generations);
    fprintf('  Search ranges: mu [%g %g], s_pos [%g %g] m, s_rot [%g %g] deg\n', ...
        o.RateRange, o.SigmaPosRange, o.SigmaRotDegRange);
    fprintf('  Stage 2: top %d + %d extra + current setting x %d seeds x %d generations\n\n', ...
        o.ValidateTop, numel(o.ExtraCandidates), o.ValidateSeeds, o.FullGenerations);

    %% Stage 1: Bayesian optimisation
    vars = [ ...
        optimizableVariable('MutationRate',     o.RateRange,        'Transform', 'log'), ...
        optimizableVariable('MutationSigmaPos', o.SigmaPosRange,    'Transform', 'log'), ...
        optimizableVariable('SigmaRotDeg',      o.SigmaRotDegRange, 'Transform', 'log')];

    evalLog = struct('MutationRate', {}, 'MutationSigmaPos', {}, 'SigmaRotDeg', {}, ...
                     'MedianJ', {}, 'Costs', {}, 'FinalDiversity', {}, 'Seconds', {});
    stage1Seeds = 1:o.NumSeeds;

    function objective = scoreCandidate(x)
        [costs, div, secs] = runSeeds(base, o, x.MutationRate, x.MutationSigmaPos, ...
            deg2rad(x.SigmaRotDeg), o.Generations, stage1Seeds);
        objective = median(costs);
        evalLog(end+1) = struct('MutationRate', x.MutationRate, ...
            'MutationSigmaPos', x.MutationSigmaPos, 'SigmaRotDeg', x.SigmaRotDeg, ...
            'MedianJ', objective, 'Costs', costs, 'FinalDiversity', div, 'Seconds', secs);
        fprintf('  eval %2d: mu=%.3f  s_pos=%.3f m  s_rot=%.1f deg  -> median J %.5f  (%.1f min)\n', ...
            numel(evalLog), x.MutationRate, x.MutationSigmaPos, x.SigmaRotDeg, ...
            objective, secs/60);
        save(fullfile(o.OutputDir, sprintf('tuneMutation_%s.mat', stamp)), 'evalLog', 'o');
    end

    bo = bayesopt(@scoreCandidate, vars, ...
        'MaxObjectiveEvaluations', o.MaxEvaluations, ...
        'NumSeedPoints', min(4, o.MaxEvaluations), ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'IsObjectiveDeterministic', false, ...
        'PlotFcn', [], 'Verbose', 0);

    %% Stage 2: validation at full budget
    [~, order] = sort([evalLog.MedianJ]);
    top = order(1:min(o.ValidateTop, numel(order)));
    cands = struct('Label', {}, 'MutationRate', {}, 'MutationSigmaPos', {}, 'SigmaRotDeg', {});
    for k = 1:numel(top)
        e = evalLog(top(k));
        cands(end+1) = struct('Label', sprintf('tuned #%d', k), 'MutationRate', e.MutationRate, ...
            'MutationSigmaPos', e.MutationSigmaPos, 'SigmaRotDeg', e.SigmaRotDeg); %#ok<AGROW>
    end
    for k = 1:numel(o.ExtraCandidates)
        x = o.ExtraCandidates(k);
        cands(end+1) = struct('Label', x.Label, 'MutationRate', x.MutationRate, ...
            'MutationSigmaPos', x.MutationSigmaPos, 'SigmaRotDeg', x.SigmaRotDeg); %#ok<AGROW>
    end
    cands(end+1) = struct('Label', 'current', 'MutationRate', base.MutationRate, ...
        'MutationSigmaPos', base.MutationSigmaPos, 'SigmaRotDeg', rad2deg(base.MutationSigmaRot));

    valSeeds = 1000 + (1:o.ValidateSeeds);      % unseen by stage 1
    fprintf('\n  Stage 2 validation (%d generations):\n', o.FullGenerations);
    for k = 1:numel(cands)
        c = cands(k);
        [costs, div] = runSeeds(base, o, c.MutationRate, c.MutationSigmaPos, ...
            deg2rad(c.SigmaRotDeg), o.FullGenerations, valSeeds);
        cands(k).Costs = costs;
        cands(k).MedianJ = median(costs);
        cands(k).IQR = iqr(costs);
        cands(k).FinalDiversity = div;
        fprintf('  %-9s mu=%.3f s_pos=%.3f m s_rot=%5.1f deg | median J %.5f  IQR %.5f\n', ...
            c.Label, c.MutationRate, c.MutationSigmaPos, c.SigmaRotDeg, ...
            cands(k).MedianJ, cands(k).IQR);
    end

    [~, best] = min([cands.MedianJ]);
    rec = cands(best);
    fprintf('\n  RECOMMENDED (%s): MutationRate = %.3f, MutationSigmaPos = %.3f, MutationSigmaRot = %.4f (%.1f deg)\n\n', ...
        rec.Label, rec.MutationRate, rec.MutationSigmaPos, deg2rad(rec.SigmaRotDeg), rec.SigmaRotDeg);

    result.Options = o;
    result.Stage1 = evalLog;
    result.BayesOpt = bo;
    result.Validation = cands;
    result.Recommended = rec;

    outFile = fullfile(o.OutputDir, sprintf('tuneMutation_%s.mat', stamp));
    save(outFile, 'result', 'evalLog', 'o');
    writeSummary(strrep(outFile, '.mat', '.txt'), result);
    fprintf('  Saved %s\n', outFile);
end

%% ------------------------------------------------------------------------
function [costs, finalDiv, secs] = runSeeds(base, o, mu, sPos, sRot, gens, seeds)
    cfg = base;
    cfg.MutationRate = mu;
    cfg.MutationSigmaPos = sPos;
    cfg.MutationSigmaRot = sRot;
    cfg.MaxGenerations = gens;
    [specs, problem, params] = buildRunSpecs(cfg, o.NumCameras, 3, ...
        o.TargetType, o.GridMode, o.Spacing, 'UseNormTable', false);

    costs = zeros(1, numel(seeds));
    finalDiv = zeros(1, numel(seeds));
    t0 = tic;
    for k = 1:numel(seeds)
        rng(seeds(k), 'twister');                % common random numbers
        [~, out] = evalc('RunGA(problem, params, specs)');
        costs(k) = out.bestsol.Cost;
        finalDiv(k) = out.popDiversity(end);
    end
    secs = toc(t0);
end

function writeSummary(file, r)
    fid = fopen(file, 'w');
    c = onCleanup(@() fclose(fid));
    o = r.Options;
    fprintf(fid, 'tuneMutation summary\n====================\n');
    fprintf(fid, 'Preset %s, %d cameras, TT%d GM%d spacing %.2f m\n', o.Preset, ...
        o.NumCameras, o.TargetType, o.GridMode, o.Spacing);
    fprintf(fid, 'Stage 1: %d evaluations x %d seeds x %d generations\n', ...
        numel(r.Stage1), o.NumSeeds, o.Generations);
    fprintf(fid, 'Stage 2: %d seeds x %d generations\n\n', o.ValidateSeeds, o.FullGenerations);
    fprintf(fid, 'Stage 1 (sorted by median J)\n');
    [~, ord] = sort([r.Stage1.MedianJ]);
    for k = ord
        e = r.Stage1(k);
        fprintf(fid, '  mu=%.3f  s_pos=%.3f m  s_rot=%5.1f deg  median J %.5f\n', ...
            e.MutationRate, e.MutationSigmaPos, e.SigmaRotDeg, e.MedianJ);
    end
    fprintf(fid, '\nStage 2 validation\n');
    for k = 1:numel(r.Validation)
        v = r.Validation(k);
        fprintf(fid, '  %-9s mu=%.3f s_pos=%.3f m s_rot=%5.1f deg | median J %.5f  IQR %.5f\n', ...
            v.Label, v.MutationRate, v.MutationSigmaPos, v.SigmaRotDeg, v.MedianJ, v.IQR);
    end
    rec = r.Recommended;
    fprintf(fid, '\nRECOMMENDED (%s): MutationRate = %.3f, MutationSigmaPos = %.3f, MutationSigmaRot = %.4f rad (%.1f deg)\n', ...
        rec.Label, rec.MutationRate, rec.MutationSigmaPos, deg2rad(rec.SigmaRotDeg), rec.SigmaRotDeg);
end
