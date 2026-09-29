function tests = testRunConfig
% runConfig + buildRunSpecs must reproduce the specs the batch scripts
% built by hand before the refactor, and reject bad configuration.
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    addProjectPaths();
end

function testLabPresetMatchesLegacySpecs(tc)
    cfg = runConfig('optitrack_lab');
    for tt = 1:2
        for gm = 1:2
            [specs, problem, params] = buildRunSpecs(cfg, 7, 3, tt, gm, 1.0, 'UseNormTable', false);
            legacy = legacySpecs(7, tt, gm, 1.0);
            for f = {'Target', 'SectionCentres', 'PreComputed', 'Resolution', 'PixelSize', ...
                     'Focal', 'FocalWide', 'Range', 'RangeWide', 'PrincipalPoint', ...
                     'WeightUncertainty', 'WeightOcclusion', 'NumPoints'}
                verifyEqual(tc, specs.(f{1}), legacy.(f{1}), sprintf('TT%d GM%d %s', tt, gm, f{1}));
            end
            verifyEqual(tc, problem.VarMin, repmat([-5 -4.5 0 -pi -pi/2 -pi], 1, 7));
            verifyEqual(tc, problem.VarMax, repmat([ 5  4.5 4.8 pi pi/2 pi], 1, 7));
            verifyEqual(tc, params.nPop, 420);
            verifyEqual(tc, [params.mu, params.sigma, params.Tournamentsize, params.pC], [0.5 0.1 3 1]);
        end
    end
end

function testOverridesApply(tc)
    cfg = runConfig('optitrack_lab', 'Volume', [-3 3; -3 3; 0 3], 'MutationRate', 0.3);
    [specs, ~, params] = buildRunSpecs(cfg, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyLessThanOrEqual(tc, max(abs(specs.Target(:, 1:2)), [], 'all'), 3);
    verifyEqual(tc, params.mu, 0.3);
    verifyEqual(tc, params.nPop, 6*6*10);
end

function testUnknownFieldErrors(tc)
    verifyError(tc, @() runConfig('optitrack_lab', 'Volumee', 1), 'runConfig:UnknownField');
    verifyError(tc, @() runConfig('nope'), 'runConfig:UnknownPreset');
end

function testLowcostNeedsHardwareValues(tc)
    cfg = runConfig('lowcost_tripod');
    verifyError(tc, @() buildRunSpecs(cfg, 7, 3, 1, 1, 1.0), 'setupHardwareSpecs:Incomplete');
end

function testNormTablePerPreset(tc)
    verifyTrue(tc, endsWith(normTableFile('optitrack_lab'), fullfile('Results', 'normTable.mat')));
    verifyTrue(tc, endsWith(normTableFile('lowcost_tripod'), 'normTable_lowcost_tripod.mat'));
end

function testBatchRunGADryRunTakesPresetAndOverrides(tc)
    args = {'DryRun', true, 'CameraRange', 7, 'CostFunctions', 3, 'TargetTypes', 1, ...
            'GridModes', 1, 'NumRepeats', 1};
    verifyWarningFree(tc, @() quiet(@() batchRunGA(args{:}, 'MutationRate', 0.3)));
    verifyError(tc, @() quiet(@() batchRunGA(args{:}, 'Bogus', 1)), 'runConfig:UnknownField');
end

function testSmokeRunGA(tc)
    cfg = runConfig('optitrack_lab', 'MaxGenerations', 2, 'PopulationSize', 12);
    [specs, problem, params] = buildRunSpecs(cfg, 4, 3, 2, 1, 2.0, 'UseNormTable', false);
    out = evalcOut(@() RunGA(problem, params, specs));
    verifyTrue(tc, isfinite(out.bestsol.Cost));
    verifyEqual(tc, numel(out.bestsol.Chromosome), problem.nVar);
end

%% Helpers
function quiet(fn) %#ok<INUSD>
    evalc('fn()');
end

function out = evalcOut(fn)
    [~, out] = evalc('fn()');
end

function specs = legacySpecs(numCams, tt, gm, sp)
% Spec construction exactly as batchRunGA did it before runConfig.
    volume = [-4 4; -4 4; 0 4];
    if tt == 2
        volume(3,:) = [0, 0.5];
        zSp = min(0.25, 0.5);
        targetSpacing = [sp, sp, zSp];
    else
        targetSpacing = sp;
    end
    specs = setupHardwareSpecs(numCams);
    specs.WeightUncertainty = 0.5;
    specs.WeightOcclusion   = 0.5;
    specs.TargetType = tt;
    specs.TargetMode = gm;
    specs.Target = generateTargetSpace(volume, gm, targetSpacing);
    specs.NumPoints = size(specs.Target, 1);
    specs.spacing = sp;
    specs.SectionCentres = generateSectionCentres(numCams, volume);
    specs.UseNormTable = false;
    specs = setupCostParams(specs);
end
