function tests = testMount
% Cameras may only sit in mountable regions (walls, ceiling, tripod ring),
% never inside the capture footprint.
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    addProjectPaths();
end

function testFootprintPositionIsProjectedOut(tc)
    regions = mountRegions(runConfig('optitrack_lab'));
    p = projectToMount([0.5 -1 1.2], regions);
    verifyTrue(tc, inAnyRegion(p, regions));
    verifyFalse(tc, inFootprintBelowCeiling(p));
end

function testProjectionIsIdempotent(tc)
    regions = mountRegions(runConfig('optitrack_lab'));
    rng(3);
    for k = 1:200
        p = projectToMount([-6 -5 -1] + rand(1,3).*[12 10 7], regions);
        verifyEqual(tc, projectToMount(p, regions), p, 'AbsTol', 1e-12);
    end
end

function testSamplesLieInRegions(tc)
    for model = {'walls_ceiling_tripod', 'walls_ceiling', 'tripod'}
        cfg = labWithMount(model{1});
        regions = mountRegions(cfg);
        rng(4);
        for k = 1:500
            verifyTrue(tc, inAnyRegion(sampleMount(regions), regions), model{1});
        end
    end
end

function testTripodOnlyNeverOnWallOrCeiling(tc)
    cfg = labWithMount('tripod');
    regions = mountRegions(cfg);
    rng(5);
    for k = 1:300
        p = projectToMount([-6 -5 -1] + rand(1,3).*[12 10 7], regions);
        verifyGreaterThanOrEqual(tc, p(3), cfg.Mount.TripodHeight(1));
        verifyLessThanOrEqual(tc, p(3), cfg.Mount.TripodHeight(2));
        verifyFalse(tc, inFootprintBelowCeiling(p));
    end
end

function testBoxModelIsUnconstrained(tc)
    regions = mountRegions(labWithMount('box'));
    verifyEmpty(tc, regions);
    verifyEqual(tc, projectToMount([0 0 1], regions), [0 0 1]);
end

function testGAKeepsCamerasMounted(tc)
    for model = {'walls_ceiling_tripod', 'tripod'}
        cfg = labWithMount(model{1});
        cfg.MaxGenerations = 3;
        cfg.PopulationSize = 16;
        [specs, problem, params] = buildRunSpecs(cfg, 4, 3, 2, 1, 2.0, 'UseNormTable', false);
        [~, out] = evalc('RunGA(problem, params, specs)');
        for i = 1:numel(out.pop)
            ch = out.pop(i).Chromosome;
            for s = 0:6:numel(ch) - 6
                verifyTrue(tc, inAnyRegion(ch(s+1:s+3), specs.MountRegions), ...
                    sprintf('%s: individual %d camera %d', model{1}, i, s/6 + 1));
            end
        end
    end
end

%% Helpers
function cfg = labWithMount(model)
    cfg = runConfig('optitrack_lab');
    cfg.Mount.Model = model;
end

function tf = inAnyRegion(p, regions)
    tol = 1e-9;
    tf = false;
    for k = 1:numel(regions)
        r = regions(k);
        inBox = all(p >= r.Lo - tol) && all(p <= r.Hi + tol);
        if inBox && strcmp(r.Type, 'ring')
            inBox = ~(p(1) > r.InnerLo(1) + tol && p(1) < r.InnerHi(1) - tol && ...
                      p(2) > r.InnerLo(2) + tol && p(2) < r.InnerHi(2) - tol);
        end
        tf = tf || inBox;
    end
end

function tf = inFootprintBelowCeiling(p)
    tf = abs(p(1)) < 4 && abs(p(2)) < 4 && p(3) < 4.8 - 1e-9;
end
