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
    cfg = runConfig('openmv_2m', 'MaxGenerations', 2, 'PopulationSize', 12);
    verifyEqual(tc, cfg.CameraRange, 6:12);
    verifyEqual(tc, cfg.Volume, [-2.5 2.5; -3 3; 0.2 2.0]);
    verifyEqual(tc, [cfg.TargetTypes cfg.GridModes], [1 1]);
    [specs, problem, params] = buildRunSpecs(cfg, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyEqual(tc, specs.HardwareProfile, 'openmv_rt1062');
    [~, out] = evalc('RunGA(problem, params, specs)');
    verifyTrue(tc, isfinite(out.bestsol.Cost));
end

function testTripodMounting(tc)
    cfg = runConfig('openmv_2m', 'MaxGenerations', 3, 'PopulationSize', 12);
    verifyEqual(tc, cfg.Mount.Model, 'tripod');
    verifyEqual(tc, cfg.Mount.TripodHeight, [0.75 1.75]);
    [specs, problem, params] = buildRunSpecs(cfg, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    [~, out] = evalc('RunGA(problem, params, specs)');
    P = reshape(out.bestsol.Chromosome, 6, []).';
    verifyGreaterThanOrEqual(tc, P(:, 3), 0.75 - 1e-9);
    verifyLessThanOrEqual(tc, P(:, 3), 1.75 + 1e-9);
    outside = abs(P(:, 1)) >= 2.5 - 1e-9 | abs(P(:, 2)) >= 3 - 1e-9;  % not in the 5 x 6 m flight space
    verifyTrue(tc, all(outside));
    verifyTrue(tc, all(abs(P(:, 1)) <= 4 + 1e-9 & abs(P(:, 2)) <= 4.5 + 1e-9));  % within the 1.5 m band
end

function testFourMetrePreset(tc)
    % first sweep: full 4 m volume, 6-10 cameras; old name maps to it
    for name = {'openmv_4m', 'openmv_lab'}
        cfg = runConfig(name{1});
        verifyEqual(tc, cfg.Preset, 'openmv_4m');
        verifyEqual(tc, cfg.Volume(3, :), [0 4]);
        verifyEqual(tc, cfg.CameraRange, 6:10);
    end
    verifyTrue(tc, endsWith(runResultsDir('openmv_4m', 7), fullfile('Results', 'openmv_4m', '7Cams')));
    verifyTrue(tc, endsWith(normTableFile('openmv_2m'), 'normTable_openmv_2m.mat'));
end

function testFlightSpaceGrid(tc)
    % 5 x 6 m flight space, 0.2-2 m: 6 x 7 x 3 points, same for both hardware
    for name = {'openmv_2m', 'optitrack_tripod'}
        specs = buildRunSpecs(runConfig(name{1}), 7, 3, 1, 1, 1.0, 'UseNormTable', false);
        verifyEqual(tc, specs.NumPoints, 126, name{1});
        verifyEqual(tc, unique(specs.Target(:, 3)).', [0.2 1.1 2.0], 'AbsTol', 1e-9);
        verifyEqual(tc, [min(specs.Target(:, 1:2)); max(specs.Target(:, 1:2))], [-2.5 -3; 2.5 3], 'AbsTol', 1e-9);
    end
    o = runConfig('optitrack_tripod');
    verifyEqual(tc, o.Hardware, 'optitrack');
    verifyEqual(tc, o.Mount.Model, 'tripod');
    verifyEqual(tc, o.CameraRange, [6 7]);
    verifyTrue(tc, endsWith(runResultsDir('optitrack_tripod', 7), fullfile('Results', 'optitrack_tripod', '7Cams')));
    % lab preset untouched
    lab = buildRunSpecs(runConfig('optitrack_lab'), 7, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyEqual(tc, lab.NumPoints, 405);
end

function testCabinetPresets(tc)
    % cabinet + floor + table-end mounting, 2 m and 4 m volumes; floor-only control
    c2 = runConfig('openmv_cab_2m');  c4 = runConfig('openmv_cab_4m');  f2 = runConfig('openmv_floor_2m');
    verifyEqual(tc, [c2.CameraRange; c4.CameraRange; f2.CameraRange], repmat(6:12, 3, 1));
    s2 = buildRunSpecs(c2, 7, 3, 1, 1, 1.0, 'UseNormTable', false);
    s4 = buildRunSpecs(c4, 7, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyEqual(tc, s2.NumPoints, 126);
    verifyEqual(tc, s4.NumPoints, 210);
    verifyEqual(tc, unique(s4.Target(:, 3)).', [0.2 1.15 2.1 3.05 4.0], 'AbsTol', 1e-9);
    reg = mountRegions(c2);
    verifyEqual(tc, numel(reg), 6);
    cab = reg(startsWith({reg.Name}, 'Cabinet'));
    verifyEqual(tc, [cab(2).Lo; cab(2).Hi], [4 -2.75 1.4; 4.75 2.75 2.9], 'AbsTol', 1e-9);
    ends = reg(startsWith({reg.Name}, 'End'));
    verifyEqual(tc, [ends(1).Lo(3) ends(1).Hi(3); ends(2).Lo(3) ends(2).Hi(3)], [0.5 2.75; 0.5 2.75], 'AbsTol', 1e-9);
    fr = mountRegions(f2);
    verifyEqual(tc, {fr.Type}, {'ring'});
    verifyEqual(tc, [fr.Lo(3) fr.Hi(3)], [0.5 2.0]);
    % a short GA keeps every camera on a mount and out of the flight space
    cfg = runConfig('openmv_cab_4m', 'MaxGenerations', 3, 'PopulationSize', 12);
    [specs, problem, params] = buildRunSpecs(cfg, 8, 3, 1, 1, 1.0, 'UseNormTable', false);
    [~, out] = evalc('RunGA(problem, params, specs)');
    P = reshape(out.bestsol.Chromosome, 6, []).';
    for k = 1:size(P, 1)
        inR = arrayfun(@(r) all(P(k, 1:3) >= r.Lo - 1e-9) && all(P(k, 1:3) <= r.Hi + 1e-9), specs.MountRegions);
        verifyTrue(tc, any(inR), sprintf('camera %d not on a mount', k));
    end
    verifyTrue(tc, all(abs(P(:, 1)) >= 2.5 - 1e-9 | abs(P(:, 2)) >= 3 - 1e-9));
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
