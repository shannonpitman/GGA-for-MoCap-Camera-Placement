function AB = abParameterisation(varargin)
%ABPARAMETERISATION  A/B test: rotation-vector vs XYZ-Euler orientation genes.
%
%   AB = abParameterisation()               % 5 paired repeats, reduced budget
%   AB = abParameterisation('Quick', true)  % smoke test
%
%   Stand-alone study; the production GA only uses rotation vectors. This
%   script runs the same GA loop as GA_Core/RunGA (tournament selection,
%   block crossover, Gaussian mutation, repair, elitism) twice per seed:
%
%     rotvec  genes 4-6 = rotation vector, wrapped to |r| <= pi
%     euler   genes 4-6 = XYZ Euler [alpha beta gamma], alpha/gamma wrapped
%             to (-pi, pi], beta clamped to [-pi/2, pi/2] (the pre-2026-09
%             bounds). Costs are evaluated after converting to rotation
%             vectors, so both arms use identical cost functions.
%
%   Both arms of a repeat start from the same physical initial population
%   (generated once, then encoded per arm) and use the same random seed for
%   the GA loop, so the comparison is paired. Mutation uses the same
%   numbers (runConfig MutationRate / MutationSigmaPos / MutationSigmaRot)
%   in each arm's own gene units, which is the point of the test.
%
%   Reported: final best J per repeat and arm, median / IQR, wins per
%   repeat, Wilcoxon signed-rank p-value, and median convergence curves.
%   Saved to Results/Tuning/ABParameterisation/.

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'Preset',         'optitrack_lab', @ischar);
    addParameter(p, 'NumCameras',     6,    @isnumeric);
    addParameter(p, 'TargetType',     1,    @isnumeric);
    addParameter(p, 'GridMode',       1,    @isnumeric);
    addParameter(p, 'Spacing',        1.0,  @isnumeric);
    addParameter(p, 'PopulationSize', 120,  @isnumeric);
    addParameter(p, 'Generations',    30,   @isnumeric);
    addParameter(p, 'Repeats',        5,    @isnumeric);
    addParameter(p, 'Quick',          false, @islogical);
    addParameter(p, 'OutputDir',      '',   @ischar);
    parse(p, varargin{:});
    o = p.Results;
    if o.Quick
        o.NumCameras = 4; o.Spacing = 2.0; o.PopulationSize = 8; o.Generations = 2; o.Repeats = 1;
    end
    if isempty(o.OutputDir)
        o.OutputDir = fullfile(addProjectPaths(), 'Results', 'Tuning', 'ABParameterisation');
    end
    if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end

    cfg = runConfig(o.Preset, 'PopulationSize', o.PopulationSize, 'MaxGenerations', o.Generations);
    [specs, problem, params] = buildRunSpecs(cfg, o.NumCameras, 3, o.TargetType, ...
        o.GridMode, o.Spacing, 'UseNormTable', false);

    arms = {'rotvec', 'euler'};
    finalJ = zeros(o.Repeats, 2);
    curves = nan(o.Generations, o.Repeats, 2);
    fprintf('\n  A/B parameterisation: %d cams, pop %d, %d generations, %d repeats\n', ...
        o.NumCameras, o.PopulationSize, o.Generations, o.Repeats);

    for s = 1:o.Repeats
        rng(s, 'twister');
        init = zeros(params.nPop, problem.nVar);
        for i = 1:params.nPop
            init(i, :) = initialPopulation(problem.VarMin, problem.VarMax, ...
                specs.SectionCentres, o.NumCameras, specs.MountRegions);
        end
        for a = 1:2
            t0 = tic;
            [best, curve] = runArm(arms{a}, init, problem, params, specs, 1000 + s);
            finalJ(s, a) = best;
            curves(:, s, a) = curve;
            fprintf('  repeat %d  %-6s  best J = %.5f  (%.1f min)\n', s, arms{a}, best, toc(t0)/60);
        end
    end

    AB.Options = o;
    AB.Arms = arms;
    AB.FinalJ = finalJ;
    AB.Curves = curves;
    AB.Median = median(finalJ, 1);
    AB.IQR = iqr(finalJ, 1);
    AB.RotvecWins = sum(finalJ(:,1) < finalJ(:,2));
    if o.Repeats > 1
        AB.SignRankP = signrank(finalJ(:,1), finalJ(:,2));
    else
        AB.SignRankP = NaN;
    end

    fprintf('\n  rotvec median J %.5f (IQR %.5f) | euler median J %.5f (IQR %.5f)\n', ...
        AB.Median(1), AB.IQR(1), AB.Median(2), AB.IQR(2));
    fprintf('  rotvec better in %d of %d repeats; Wilcoxon signed-rank p = %.3g\n', ...
        AB.RotvecWins, o.Repeats, AB.SignRankP);

    stamp = string(datetime('now'), 'yyyyMMdd_HHmmss');
    outFile = fullfile(o.OutputDir, sprintf('abParameterisation_%s.mat', stamp));
    save(outFile, 'AB');
    fid = fopen(strrep(outFile, '.mat', '.txt'), 'w');
    fprintf(fid, 'A/B parameterisation: %d cams, pop %d, %d generations, %d repeats, TT%d GM%d sp %.2f\n', ...
        o.NumCameras, o.PopulationSize, o.Generations, o.Repeats, o.TargetType, o.GridMode, o.Spacing);
    for s = 1:o.Repeats
        fprintf(fid, '  repeat %d: rotvec %.5f | euler %.5f\n', s, finalJ(s,1), finalJ(s,2));
    end
    fprintf(fid, 'rotvec median %.5f (IQR %.5f) | euler median %.5f (IQR %.5f)\n', ...
        AB.Median(1), AB.IQR(1), AB.Median(2), AB.IQR(2));
    fprintf(fid, 'rotvec better in %d of %d; Wilcoxon signed-rank p = %.3g\n', ...
        AB.RotvecWins, o.Repeats, AB.SignRankP);
    fclose(fid);
    fprintf('  Saved %s\n', outFile);
end

%% ------------------------------------------------------------------------
function [bestCost, curve] = runArm(arm, initRotvec, problem, params, specs, seed)
% Same loop as RunGA, with encode/decode around the orientation genes.
    numCams = specs.Cams;
    nPop = params.nPop;
    nC = round(params.pC*nPop/2)*2;
    regions = specs.MountRegions;

    encode = @(c) c;  decode = @(c) c;
    if strcmp(arm, 'euler')
        encode = @toEuler;  decode = @eulerChromToRotvec;
    end

    pop = zeros(nPop, problem.nVar);
    for i = 1:nPop
        pop(i, :) = encode(initRotvec(i, :));
    end
    costs = evaluate(pop, decode, problem, specs);

    rng(seed, 'twister');
    curve = nan(params.MaxIt, 1);
    popStruct = struct('Chromosome', num2cell(pop, 2), 'Cost', num2cell(costs));
    for it = 1:params.MaxIt
        parents = Tournament(popStruct, problem.nVar, nC, params.Tournamentsize);
        kids = zeros(nC, problem.nVar);
        for k = 1:nC/2
            [kids(2*k-1, :), kids(2*k, :)] = DoublePointCrossover(parents(2*k-1, :), parents(2*k, :), numCams);
        end
        for l = 1:nC
            y = Mutate(kids(l, :), params.mu, params.sigma);
            if mod(it, 10) == 0 || rand < 0.2
                y = encode(fixPoorCameras(decode(y), specs, 0.05));
            end
            y = normalise(y, arm, problem);
            y = enforcePositions(y, regions);
            kids(l, :) = y;
        end
        kidCosts = evaluate(kids, decode, problem, specs);

        all = [pop; kids];  allC = [costs; kidCosts];
        [allC, ord] = sort(allC);
        pop = all(ord(1:nPop), :);  costs = allC(1:nPop);
        popStruct = struct('Chromosome', num2cell(pop, 2), 'Cost', num2cell(costs));
        curve(it) = costs(1);
    end
    bestCost = costs(1);
end

function c = evaluate(pop, decode, problem, specs)
    n = size(pop, 1);
    c = zeros(n, 1);
    f = problem.CostFunction;
    parfor i = 1:n
        c(i) = f(decode(pop(i, :)), specs);
    end
end

function y = normalise(y, arm, problem)
    if strcmp(arm, 'rotvec')
        y = wrapChromRotvecs(y);
        y = min(max(y, problem.VarMin), problem.VarMax);
    else
        for s = 0:6:numel(y) - 6
            y(s+4) = wrapAnglePi(y(s+4));
            y(s+5) = min(max(y(s+5), -pi/2), pi/2);   % old beta bounds
            y(s+6) = wrapAnglePi(y(s+6));
        end
        posMin = problem.VarMin;  posMax = problem.VarMax;
        for s = 0:6:numel(y) - 6
            y(s+1:s+3) = min(max(y(s+1:s+3), posMin(s+1:s+3)), posMax(s+1:s+3));
        end
    end
end

function y = enforcePositions(y, regions)
    for s = 0:6:numel(y) - 6
        y(s+1:s+3) = projectToMount(y(s+1:s+3), regions);
    end
end

function c = toEuler(c)
    for s = 0:6:numel(c) - 6
        c(s+4:s+6) = rotm2eul(genesToRotm(c(s+4:s+6)), "XYZ");
    end
end
