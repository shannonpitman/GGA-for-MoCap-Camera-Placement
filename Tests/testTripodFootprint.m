function tests = testTripodFootprint
% Tripod footprint (*_fp presets): tripod legs stay clear of walls and
% cabinet edges, no two tripods stand closer than MinTripodSpacing, and no
% camera looks below MinElevation (OpenMV ball head: level or up).
    tests = functiontests(localfunctions);
end

function setupOnce(tc) %#ok<INUSD>
    addProjectPaths();
end

function testDefaultsUnchanged(tc)
    % every existing preset keeps point tripods and no spacing rule
    for p = {'optitrack_lab', 'openmv_2m', 'optitrack_tripod', 'openmv_cab_2m', 'openmv_cab_4m', 'openmv_floor_2m'}
        c = runConfig(p{1});
        verifyEqual(tc, [c.Mount.TripodFootRadius c.Mount.MinTripodSpacing c.Mount.MinElevation], [0 0 -90], p{1});
    end
    reg = mountRegions(runConfig('openmv_cab_2m'));
    cab = reg(strcmp({reg.Name}, 'Cabinet +x'));
    verifyEqual(tc, [cab.Lo; cab.Hi], [4 -2.75 1.4; 4.75 2.75 2.9], 'AbsTol', 1e-9);
end

function testFootprintPresets(tc)
    c = runConfig('openmv_cab_2m_fp');
    verifyEqual(tc, c.Preset, 'openmv_cab_2m_fp');               % own results folder
    verifyEqual(tc, [c.Mount.TripodFootRadius c.Mount.MinTripodSpacing c.Mount.MinElevation], [0.3 0.6 0]);
    verifyEqual(tc, buildRunSpecs(c, 6, 3, 1, 1, 1.0, 'UseNormTable', false).MinElevation, 0);
    verifyEqual(tc, c.Volume, runConfig('openmv_cab_2m').Volume);
    verifyEqual(tc, runConfig('openmv_cab_4m_fp').Volume(3, :), [0.2 4.0]);
    verifyEqual(tc, runConfig('openmv_floor_2m_fp').Mount.Model, 'tripod');
    reg = mountRegions(c);
    get = @(n) reg(strcmp({reg.Name}, n));
    cab = get('Cabinet +x');
    verifyEqual(tc, [cab.Lo(1:2); cab.Hi(1:2)], [4.3 -2.45; 4.45 2.45], 'AbsTol', 1e-9);
    band = get('Floor band +x');
    verifyEqual(tc, [band.Lo(1:2); band.Hi(1:2)], [2.5 -4.2; 3.7 4.2], 'AbsTol', 1e-9);
    e = get('End +y (floor/table)');
    verifyEqual(tc, [e.Lo(2) e.Hi(2)], [3 4.2], 'AbsTol', 1e-9);  % flight-space edge not inset
    ring = mountRegions(runConfig('openmv_floor_2m_fp'));
    verifyEqual(tc, [ring.Lo(1:2); ring.Hi(1:2)], [-3.7 -4.2; 3.7 4.2], 'AbsTol', 1e-9);
end

function testTooWideForCabinet(tc)
    c = runConfig('openmv_cab_2m', 'Mount', struct('Model', 'tripod_cabinets', 'TripodHeight', [0.5 2], ...
        'Regions', [], 'TripodFootRadius', 0.4, 'MinTripodSpacing', 0, 'EndHeight', [0.5 2.75], ...
        'Cabinet', struct('XDepth', [4 4.75], 'YSpan', [-2.75 2.75], 'Height', 0.9)));
    verifyError(tc, @() mountRegions(c), 'mountRegions:TripodTooWide');
end

function testSeparateCoincident(tc)
    % two cameras on the same spot on the cabinet end up 0.6 m apart, still on it
    c = runConfig('openmv_cab_2m_fp');  reg = mountRegions(c);
    chrom = [4.4 0 2.9 0 0 0, 4.4 0 2.9 0 0 0];
    out = enforceMount(chrom, reg, 0.6);
    P = reshape(out, 6, []).';
    verifyGreaterThanOrEqual(tc, norm(P(1, 1:2) - P(2, 1:2)), 0.6 - 1e-6);
    verifyTrue(tc, all(arrayfun(@(k) onMount(P(k, 1:3), reg), 1:2)));
    verifyEqual(tc, out([4:6 10:12]), chrom([4:6 10:12]));    % orientation untouched
end

function testSeparateRandomRigs(tc)
    % 12 cameras sampled anywhere on the mounts are always repaired
    rng(3);
    for p = {'openmv_cab_2m_fp', 'openmv_floor_2m_fp'}
        c = runConfig(p{1});  reg = mountRegions(c);
        for t = 1:30
            chrom = zeros(1, 72);
            for k = 1:12, chrom(6*k-5:6*k-3) = sampleMount(reg); end
            [out, nViol] = separateTripods(enforceMount(chrom, reg), reg, 0.6);
            P = reshape(out, 6, []).';
            verifyEqual(tc, nViol, 0, sprintf('%s trial %d', p{1}, t));
            verifyGreaterThanOrEqual(tc, min(pairDistances(P(:, 1:2))), 0.6 - 1e-6);
            verifyTrue(tc, all(arrayfun(@(k) onMount(P(k, 1:3), reg), 1:12)));
        end
    end
end

function testGAKeepsFootprint(tc)
    cfg = runConfig('openmv_cab_2m_fp', 'MaxGenerations', 3, 'PopulationSize', 12);
    [specs, problem, params] = buildRunSpecs(cfg, 10, 3, 1, 1, 1.0, 'UseNormTable', false);
    verifyEqual(tc, specs.MinTripodSpacing, 0.6);
    [~, out] = evalc('RunGA(problem, params, specs)');
    P = reshape(out.bestsol.Chromosome, 6, []).';
    verifyGreaterThanOrEqual(tc, min(pairDistances(P(:, 1:2))), 0.6 - 1e-6);
    verifyTrue(tc, all(arrayfun(@(k) onMount(P(k, 1:3), specs.MountRegions), 1:size(P, 1))));
    verifyGreaterThanOrEqual(tc, min(elevations(out.bestsol.Chromosome)), -1e-6);
    for k = 1:numel(out.pop)                    % whole final population, not only the best
        verifyGreaterThanOrEqual(tc, min(elevations(out.pop(k).Chromosome)), -1e-6);
    end
end

function testElevationLimit(tc)
    % random orientations: low cameras end exactly level with the same
    % azimuth and horizon roll; legal cameras are untouched
    rng(5);
    n = 200;
    chrom = zeros(1, 6*n);
    for k = 1:n
        v = randn(1, 3);  chrom(6*k-2:6*k) = v / norm(v) * pi * rand;
    end
    chrom(4:6) = [pi 0 0];                      % straight down
    [out, nCh] = limitElevation(chrom, 0);
    e0 = elevations(chrom);  e1 = elevations(out);
    verifyEqual(tc, nCh, nnz(e0 < 0));
    verifyGreaterThanOrEqual(tc, min(e1), -1e-9);
    low = e0 < 0;
    verifyEqual(tc, e1(low), zeros(nnz(low), 1), 'AbsTol', 1e-9);
    verifyEqual(tc, out(repelem(~low, 6)'), chrom(repelem(~low, 6)'));
    for k = find(low(:)).'
        if k == 1, continue; end                % straight down has no azimuth
        R0 = genesToRotm(chrom(6*k-2:6*k));  R1 = genesToRotm(out(6*k-2:6*k));
        a0 = atan2(R0(2,3), R0(1,3));  a1 = atan2(R1(2,3), R1(1,3));
        verifyEqual(tc, wrapToPi(a1 - a0), 0, 'AbsTol', 1e-9);
        % pure pitch: camera axes keep their components along the
        % horizontal pitch axis, so the roll relative to the horizon is kept
        a = [-sin(a0); cos(a0); 0];
        verifyEqual(tc, a.' * R1, a.' * R0, 'AbsTol', 1e-9);
    end
    verifyEqual(tc, limitElevation(chrom, -90), chrom);
    verifyEqual(tc, out(1:3), chrom(1:3));      % positions untouched
end

function testLegaliseRig(tc)
    c = runConfig('openmv_cab_2m');
    specs = buildRunSpecs(c, 6, 3, 1, 1, 1.0, 'UseNormTable', false);
    % two tripods 0.2 m apart on the cabinet's back edge, the rest at the ends
    chrom = [4.75 0 2.9 0 0 0, 4.75 0.2 2.9 0 0 0, -4.75 0 2.9 0 0 0, ...
             0 4.5 1 0 0 0, 0 -4.5 1 0 0 0, 2 4.5 1 0 0 0];
    [out, info] = legaliseRig(chrom, specs, 'openmv_cab_2m_fp');
    P = reshape(out, 6, []).';
    verifyEqual(tc, info.Violations, 0);
    verifyLessThanOrEqual(tc, max(abs(P(1:3, 1))), 4.45 + 1e-9);   % off the back edge
    verifyLessThanOrEqual(tc, max(abs(P(4:6, 2))), 4.2 + 1e-9);    % off the end walls
    verifyEqual(tc, size(info.RawOcc_deg), [1 2]);
end

function e = elevations(chrom)
    P = reshape(chrom, 6, []).';
    e = zeros(size(P, 1), 1);
    for k = 1:size(P, 1)
        R = genesToRotm(P(k, 4:6));  e(k) = asind(max(-1, min(1, R(3, 3))));
    end
end

function tf = onMount(p, reg)
    tf = any(arrayfun(@(r) all(p >= r.Lo - 1e-9) && all(p <= r.Hi + 1e-9), reg));
end
