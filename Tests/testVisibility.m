function tests = testVisibility
% Both cost terms, the coverage statistics and the repair operator must use
% the same visibility test (FOV + in front + effective range).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    addProjectPaths();
    specs = setupHardwareSpecs(7);
    specs.WeightUncertainty = 0.5;
    specs.WeightOcclusion   = 0.5;
    specs.TargetType = 1;
    specs.TargetMode = 1;
    specs.Target     = generateTargetSpace([-4 4; -4 4; 0 4], 1, 1.0);
    specs.NumPoints  = size(specs.Target, 1);
    specs.spacing    = 1.0;
    specs.UseNormTable = false;
    specs = setupCostParams(specs);
    tc.TestData.specs = specs;
    tc.TestData.chrom = buildOptiTrackChromosome();
end

function testCostPenaltiesMatchCoverage(tc)
    specs = tc.TestData.specs;
    [cams, cc] = setupCameras(tc.TestData.chrom, specs.Cams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    [~, unc] = resUncertainty(specs, cams, cc);
    [~, occ] = dynamicOcclusion(specs, cams, cc);
    cov = perTargetCoverage(tc.TestData.chrom, specs);
    pen = specs.PreComputed.penaltyUncertainty;

    verifyTrue(tc, any(cov == 0) && any(cov == 1), 'fixture should include 0- and 1-camera points');
    verifyEqual(tc, unc(cov == 0), repmat(pen, nnz(cov == 0), 1));
    verifyEqual(tc, unc(cov == 1), repmat(0.5*pen, nnz(cov == 1), 1));
    verifyEqual(tc, occ(cov == 0), repmat(720, nnz(cov == 0), 1));
    verifyEqual(tc, occ(cov == 1), repmat(360, nnz(cov == 1), 1));
end

function testRepairCountsMatchCoverage(tc)
    specs = tc.TestData.specs;
    [cams, cc] = setupCameras(tc.TestData.chrom, specs.Cams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    visMask = projectVisibilityOcclusion(cams, specs.Target, cc, specs.Resolution, ...
        specs.PreComputed.maxCameraRange, specs.PreComputed.maxCameraRangeWide, specs.FocalWide);
    verifyEqual(tc, cameraCoverageCounts(tc.TestData.chrom, specs), sum(visMask, 1).');
end

function testRangeAffectsResolutionCost(tc)
    specs = tc.TestData.specs;
    [cams, cc] = setupCameras(tc.TestData.chrom, specs.Cams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    full = resUncertainty(specs, cams, cc);
    specs.PreComputed.maxCameraRange     = 3;
    specs.PreComputed.maxCameraRangeWide = 3;
    short = resUncertainty(specs, cams, cc);
    verifyGreaterThan(tc, short, full);
end
