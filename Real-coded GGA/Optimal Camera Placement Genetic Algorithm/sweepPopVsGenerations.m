%% sweepPopVsGenerations.m
% Population-size vs generation-count study for the camera-placement GGA.

clc; clear; close all;
addProjectPaths();

%% USER CONFIG  %%
% Configurations to compare: {label, nPop, MaxIt}
configs = {
    '420x50', 420, 50    
    '420x100',  420,  100
};

nSeeds = 3;  % independent repeats per config (>=3 recommended)
seeds  = 1:nSeeds; % rng seeds

% Elitism handling (see header note)
ELITISM_MODE = 'fixed';    % 'scaled' -> elitismDelay = round(frac*MaxIt)
                            % 'fixed'  -> keep setupGAparams value (0)
elitismFrac  = 1/3;

% Problem definition (mirrors runCameraOptimiser.m -- keep in sync)
numCams          = 7;
volume           = [-4 4; -4 4; 0 4];
cameraLowerBounds = [-5 -4.5 0  -pi -pi -pi];
cameraUpperBounds = [ 5  4.5 4.8 pi  pi  pi];
targetType       = 2;       % 1=UAV, 2=UGV
targetMode       = 1;       % 1=uniform grid, 2=normalised grid
spacing          = 1;
UGV_maxHeight    = 0.5;
UGV_zSpacing     = 0.25;
costFunctionType = 3;       % 3 = combined
weightUncertainty = 0.5;
weightOcclusion   = 0.5;
%% --------------------------------------------------------------------- %%

% ---- Build target spacing (mirrors runCameraOptimiser) ----
if targetType == 2
    volume(3,:)   = [0, UGV_maxHeight];
    targetSpacing = [spacing, spacing, min(UGV_zSpacing, UGV_maxHeight)];
else
    targetSpacing = spacing;
end

% ---- Problem ----
problem = setupProblem(numCams, costFunctionType, cameraUpperBounds, cameraLowerBounds);

% Specs
specs = setupHardwareSpecs(numCams);
specs.warmStart = false;
specs.warmChromosomes = [];
specs.WeightUncertainty = weightUncertainty;
specs.WeightOcclusion = weightOcclusion;
specs.TargetType = targetType;
specs.TargetMode = targetMode;
specs.Target = generateTargetSpace(volume, targetMode, targetSpacing);
specs.NumPoint = size(specs.Target,1);
specs.spacing = spacing;
if targetType == 2
    specs.spacingZ = UGV_zSpacing;
end
specs.SectionCentres = generateSectionCentres(numCams, volume);
specs = setupCostParams(specs);

% ---- Output folder ----
stamp   = datestr(now,'yyyymmdd_HHMMSS');
outDir  = fullfile('Sweep_PopVsGen', ['sweep_' stamp]);
if ~exist(outDir,'dir'); mkdir(outDir); end

nConfig = size(configs,1);
R = struct();   % results container, one element per config

fprintf('=== Population vs Generations sweep ===\n');
fprintf('numCams=%d, costFn=%d, nSeeds=%d, elitism=%s\n\n', ...
    numCams, costFunctionType, nSeeds, ELITISM_MODE);

tSweep = tic;
for c = 1:nConfig
    label = configs{c,1};
    nPop  = configs{c,2};
    MaxIt = configs{c,3};

    params = setupGAparams(MaxIt, nPop);
    switch ELITISM_MODE
        case 'scaled'
            params.elitismDelay = round(elitismFrac*MaxIt);
        case 'fixed'
            % keep setupGAparams default
        otherwise
            error('Unknown ELITISM_MODE');
    end

    budget = nPop*(1+MaxIt);   % initial pop + nPop offspring per gen (pC=1)
    fprintf('[%d/%d] %s  (nPop=%d, MaxIt=%d, elitismDelay=%d, budget~%d evals)\n', ...
        c, nConfig, label, nPop, MaxIt, params.elitismDelay, budget);

    bestHist = nan(MaxIt, nSeeds);   % best cost per generation, per seed
    divHist  = nan(MaxIt, nSeeds);   % diversity per generation, per seed
    finalBest = nan(1, nSeeds);      % final best cost per seed
    runTime   = nan(1, nSeeds);

    for s = 1:nSeeds
        rng(seeds(s));               % controls initial pop + operators
        tRun = tic;
        out  = RunGA(problem, params, specs);
        runTime(s) = toc(tRun);

        bestHist(:,s) = out.bestcost(:);
        divHist(:,s)  = out.popDiversity(:);
        finalBest(s)  = out.bestsol.Cost;
        fprintf('    seed %d/%d: final best = %.6f (%.1f s)\n', ...
            s, nSeeds, finalBest(s), runTime(s));
    end

    R(c).label     = label;
    R(c).nPop      = nPop;
    R(c).MaxIt     = MaxIt;
    R(c).budget    = budget;
    R(c).evalAxis  = (nPop + nPop*(1:MaxIt))';  % initial pop + cumulative offspring evals
    R(c).bestHist  = bestHist;
    R(c).divHist   = divHist;
    R(c).finalBest = finalBest;
    R(c).runTime   = runTime;
    R(c).elitismDelay = params.elitismDelay;
end
fprintf('\nSweep complete in %.1f min.\n\n', toc(tSweep)/60);

%% ---------------------------- SUMMARY --------------------------------- %%
fprintf('%-10s %6s %6s %10s %10s %10s %10s\n', ...
    'config','nPop','MaxIt','budget','meanBest','stdBest','minBest');
for c = 1:nConfig
    fb = R(c).finalBest;
    fprintf('%-10s %6d %6d %10d %10.5f %10.5f %10.5f\n', ...
        R(c).label, R(c).nPop, R(c).MaxIt, R(c).budget, ...
        mean(fb), std(fb), min(fb));
end
fprintf('\n(lower cost = better)\n');

%%  PLOTS
cmap = lines(nConfig);

% (1) Convergence vs generation
fig1 = figure('Name','Convergence vs generation','Color','w','Position',[80 80 720 480]);
hold on;
h = gobjects(nConfig,1);
for c = 1:nConfig
    m = mean(R(c).bestHist,2);
    sd = std(R(c).bestHist,0,2);
    g = (1:R(c).MaxIt)';
    fill([g;flipud(g)], [m-sd;flipud(m+sd)], cmap(c,:), ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    h(c) = plot(g, m, '-', 'Color', cmap(c,:), 'LineWidth', 1.8);
end
grid on; box on;
xlabel('Generation'); ylabel('Best cost (mean \pm std)');
title('Convergence vs generation');
legend(h, {R.label}, 'Location','northeast');
saveas(fig1, fullfile(outDir,'convergence_vs_generation.png'));

% (2) Convergence vs evaluation budget (fair-budget comparison)
fig2 = figure('Name','Convergence vs evaluations','Color','w','Position',[120 120 720 480]);
hold on;
for c = 1:nConfig
    m = mean(R(c).bestHist,2);
    plot(R(c).evalAxis, m, '-', 'Color', cmap(c,:), 'LineWidth', 1.8);
end
grid on; box on;
xlabel('Cumulative offspring evaluations'); ylabel('Best cost (mean)');
title('Convergence vs evaluation budget');
legend({R.label}, 'Location','northeast');
saveas(fig2, fullfile(outDir,'convergence_vs_evaluations.png'));

% (3) Diversity vs generation
fig3 = figure('Name','Diversity vs generation','Color','w','Position',[160 160 720 480]);
hold on;
for c = 1:nConfig
    m = mean(R(c).divHist,2);
    plot((1:R(c).MaxIt)', m, '-', 'Color', cmap(c,:), 'LineWidth', 1.8);
end
grid on; box on;
xlabel('Generation'); ylabel('Population diversity (mean per-locus std)');
title('Population diversity vs generation');
legend({R.label}, 'Location','northeast');
saveas(fig3, fullfile(outDir,'diversity_vs_generation.png'));

% (4) Boxplot of final best cost per config
fig4 = figure('Name','Final best cost','Color','w','Position',[200 200 720 480]);
allFinal = []; grp = {};
for c = 1:nConfig
    allFinal = [allFinal, R(c).finalBest];
    grp = [grp, repmat({R(c).label}, 1, numel(R(c).finalBest))];
end
boxplot(allFinal, grp);
grid on; box on;
ylabel('Final best cost'); xlabel('Configuration');
title(sprintf('Final best cost across %d seeds (lower = better)', nSeeds));
saveas(fig4, fullfile(outDir,'final_cost_boxplot.png'));

%% SAVE
meta = struct('numCams',numCams,'costFunctionType',costFunctionType, ...
    'nSeeds',nSeeds,'seeds',seeds,'ELITISM_MODE',ELITISM_MODE, ...
    'elitismFrac',elitismFrac,'timestamp',stamp);
save(fullfile(outDir,'sweep_results.mat'), 'R', 'meta', 'configs');
fprintf('\nSaved results + figures to: %s\n', outDir);