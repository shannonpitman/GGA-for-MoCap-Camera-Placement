function tests = testRotvec
% Exponential-map (rotation-vector) orientation genes.
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    addProjectPaths();
end

function testRoundTrip(tc)
    rng(11);
    for k = 1:500
        R = randomRotm();
        verifyEqual(tc, genesToRotm(rotmToGenes(R)), R, 'AbsTol', 1e-9);
    end
end

function testNearPiAndIdentity(tc)
    for theta = [0, 1e-10, pi - 1e-8, pi]
        for axis = {[1 0 0], [0 1 0], [0 0 1], [1 1 1]/sqrt(3)}
            R = axang2rotm([axis{1}, theta]);
            verifyEqual(tc, genesToRotm(rotmToGenes(R)), R, 'AbsTol', 1e-7);
            verifyLessThanOrEqual(tc, norm(rotmToGenes(R)), pi + 1e-12);
        end
    end
end

function testWrapKeepsRotation(tc)
    rng(12);
    for k = 1:300
        r = (rand(1,3) - 0.5) * 4*pi;          % |r| up to ~2*sqrt(3)*pi
        w = rotvecWrap(r);
        verifyLessThanOrEqual(tc, norm(w), pi + 1e-12);
        verifyEqual(tc, genesToRotm(w), genesToRotm(r), 'AbsTol', 1e-9);
    end
end

function testMatchesToolboxAxang(tc)
    rng(13);
    for k = 1:100
        r = rotvecWrap((rand(1,3) - 0.5) * 2*pi);
        th = norm(r);
        verifyEqual(tc, genesToRotm(r), axang2rotm([r/th, th]), 'AbsTol', 1e-12);
    end
end

function testAimPointsOpticalAxis(tc)
    rng(14);
    dirs = [randn(50, 3); 0 0 1; 0 0 -1];
    for k = 1:size(dirs, 1)
        d = dirs(k, :) / norm(dirs(k, :));
        R = genesToRotm(aimRotvec(d));
        verifyEqual(tc, R(:, 3).', d, 'AbsTol', 1e-12);
    end
end

function testLegacyEulerRigScoresIdentically(tc)
% The same physical rig must cost the same whether stored as legacy XYZ
% Euler genes (converted) or built directly from the rotation matrices.
    specs = setupHardwareSpecs(7);
    specs.WeightUncertainty = 0.5; specs.WeightOcclusion = 0.5;
    specs.TargetType = 1; specs.TargetMode = 1;
    specs.Target = generateTargetSpace([-4 4; -4 4; 0 4], 1, 2.0);
    specs.NumPoints = size(specs.Target, 1); specs.spacing = 2.0;
    specs.UseNormTable = false;
    specs = setupCostParams(specs);

    rv = buildOptiTrackChromosome();           % rotation-vector genes
    eulChrom = rv;
    for s = 0:6:numel(rv) - 6
        eulChrom(s+4:s+6) = rotm2eul(genesToRotm(rv(s+4:s+6)), "XYZ");
    end
    converted = eulerChromToRotvec(eulChrom);

    [c1, cc1] = setupCameras(rv, 7, specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    [c2, cc2] = setupCameras(converted, 7, specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    verifyEqual(tc, resUncertainty(specs, c2, cc2), resUncertainty(specs, c1, cc1), 'AbsTol', 1e-9);
    verifyEqual(tc, dynamicOcclusion(specs, c2, cc2), dynamicOcclusion(specs, c1, cc1), 'AbsTol', 1e-9);
end

function testLegacyRunLoadsAsRotvec(tc)
    rv = buildOptiTrackChromosome();
    eulChrom = rv;
    for s = 0:6:numel(rv) - 6
        eulChrom(s+4:s+6) = rotm2eul(genesToRotm(rv(s+4:s+6)), "XYZ");
    end
    sd.BestSolution.Chromosome = eulChrom;
    sd.Specifications = struct('Cams', 7);          % legacy: no Parameterisation
    got = runChromosome(sd);
    for s = 0:6:numel(rv) - 6
        verifyEqual(tc, genesToRotm(got(s+4:s+6)), genesToRotm(rv(s+4:s+6)), 'AbsTol', 1e-9);
    end
    sd.BestSolution.Chromosome = rv;
    sd.Specifications.Parameterisation = 'rotvec';
    verifyEqual(tc, runChromosome(sd), rv);
end

function testUprightRollKeepsOpticalAxis(tc)
    rv = buildOptiTrackChromosome();
    [up, rep] = uprightCameras(rv);
    info = cameraOrientationInfo(up, 7);
    verifyFalse(tc, any(info.Inverted & ~info.Degenerate));
    for s = 0:6:numel(rv) - 6
        R0 = genesToRotm(rv(s+4:s+6));  R1 = genesToRotm(up(s+4:s+6));
        verifyEqual(tc, R1(:, 3), R0(:, 3), 'AbsTol', 1e-9);   % same optical axis
    end
    verifyGreaterThanOrEqual(tc, rep.NumChanged, 0);
end

function testSnapProducesGridAngles(tc)
    rv = buildOptiTrackChromosome();
    snapped = snapChromosome(rv, 'OrientationStepDeg', 10, 'KeepUpright', false);
    eul = chromToEulerDeg(snapped);
    verifyEqual(tc, mod(eul + 1e-6, 10), zeros(size(eul)), 'AbsTol', 1e-4);
end

%% Helpers
function R = randomRotm()
    q = randn(1, 4); q = q / norm(q);
    R = quat2rotm(q);
end
