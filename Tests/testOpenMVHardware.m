function tests = testOpenMVHardware
% The OpenMV RT1062 profile must give the camera the GA is meant to model:
% full-sensor pinhole with the measured focal length, one lens, 8 m range.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    addProjectPaths();
    tc.TestData.specs = setupHardwareSpecs(7, 'openmv_rt1062');
end

function testFieldOfView(tc)
    s = tc.TestData.specs;
    fpx = s.Focal / s.PixelSize;
    hfov = 2*atand(s.Resolution(1)/2 / fpx);
    vfov = 2*atand(s.Resolution(2)/2 / fpx);
    verifyEqual(tc, hfov, 58.2, 'AbsTol', 0.1);
    verifyEqual(tc, vfov, 45.3, 'AbsTol', 0.1);
    verifyEqual(tc, s.PrincipalPoint, [1296 972]);
end

function testSingleLens(tc)
    s = tc.TestData.specs;
    verifyTrue(tc, s.SingleLens);
    verifyEqual(tc, s.FocalWide, s.Focal);
    verifyEqual(tc, s.RangeWide, s.Range);
    cams = makeCams(s, zeros(1, 6*7));
    verifyEqual(tc, cellfun(@(c) c.f, cams), repmat(s.Focal, 7, 1));
end

function testProjection(tc)
    s = tc.TestData.specs;
    cam = makeCams(s, zeros(1, 6));            % at origin, looking along +z
    cam = cam{1};
    uv = cam.project([0 0 5]);
    verifyEqual(tc, uv(:).', [1296 972], 'AbsTol', 1e-6);
    fpx = s.Focal / s.PixelSize;
    x = 5 * (1296 - 1) / fpx;                  % 1 px inside the right edge
    uv = cam.project([x 0 5]);
    verifyEqual(tc, uv(1), 2591, 'AbsTol', 1e-3);
end

function testRangeCutoff(tc)
    s = tc.TestData.specs;
    [cams, cc] = setupCameras(zeros(1, 6), 1, s.Resolution, s.Focal, s.FocalWide, ...
        s.PrincipalPoint, s.PixelSize);
    P = [0 0 7.9; 0 0 8.1; 0.5 0.5 3];
    vis = projectVisibilityOcclusion(cams, P, cc, s.Resolution, s.Range, s.RangeWide, s.FocalWide);
    verifyEqual(tc, vis(:).', [true false true]);
end

function testFocalMatchesCalibration(tc)
    d = fullfile(addProjectPaths(), '..', '..', 'Motion-Capture', 'MATLAB', 'calibration', 'intrinsics');
    f = dir(fullfile(d, 'intrinsics_cam*.mat'));
    assumeNotEmpty(tc, f, 'Motion-Capture intrinsics not on this machine');
    fx = zeros(numel(f), 1);
    for k = 1:numel(f)
        S = load(fullfile(f(k).folder, f(k).name), 'intr');
        fx(k) = S.intr.focalLength(1);
        verifyEqual(tc, S.intr.imageSize, [1944 2592], f(k).name);
    end
    s = tc.TestData.specs;
    verifyEqual(tc, s.Focal / s.PixelSize, mean(fx), 'RelTol', 0.01);
end

function testPresetBuildsAndRuns(tc)
    cfg = runConfig('openmv_lab', 'MaxGenerations', 2, 'PopulationSize', 12);
    verifyEqual(tc, cfg.CameraRange, 6:10);
    verifyEqual(tc, [cfg.TargetTypes cfg.GridModes], [1 1]);
    [specs, problem, params] = buildRunSpecs(cfg, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyEqual(tc, specs.HardwareProfile, 'openmv_rt1062');
    [~, out] = evalc('RunGA(problem, params, specs)');
    verifyTrue(tc, isfinite(out.bestsol.Cost));
end

function testTripodMounting(tc)
    cfg = runConfig('openmv_lab', 'MaxGenerations', 3, 'PopulationSize', 12);
    verifyEqual(tc, cfg.Mount.Model, 'tripod');
    verifyEqual(tc, cfg.Mount.TripodHeight, [0.75 1.75]);
    [specs, problem, params] = buildRunSpecs(cfg, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    [~, out] = evalc('RunGA(problem, params, specs)');
    P = reshape(out.bestsol.Chromosome, 6, []).';
    verifyGreaterThanOrEqual(tc, P(:, 3), 0.75 - 1e-9);
    verifyLessThanOrEqual(tc, P(:, 3), 1.75 + 1e-9);
    outside = abs(P(:, 1)) >= 4 - 1e-9 | abs(P(:, 2)) >= 4 - 1e-9;   % not in the 8x8 footprint
    verifyTrue(tc, all(outside));
end

function testOptiTrackStillAlternates(tc)
    s = setupHardwareSpecs(4, 'optitrack');
    verifyFalse(tc, s.SingleLens);
    cams = makeCams(s, zeros(1, 24));
    verifyEqual(tc, cellfun(@(c) c.f, cams), [s.Focal; s.FocalWide; s.Focal; s.FocalWide]);
end

function cams = makeCams(s, chrom)
    cams = setupCameras(chrom, numel(chrom)/6, s.Resolution, s.Focal, s.FocalWide, ...
        s.PrincipalPoint, s.PixelSize);
end
