%% testProjectAllPointsEquivalence.m
% Verifies the batched projectAllPoints path produces IDENTICAL per-point
% uncertainties to the original per-point cameras{i}.project path, then
% times the projection speedup.
%
% Equivalence is checked by running computePointUncertainty two ways on the
% same points:
%   NEW: precomputed U/V/visMask from projectAllPoints (trailing args)
%   OLD: no trailing args -> function projects each point itself
% If these agree to ~1e-9 across many random camera layouts, the batched
% projection is a drop-in replacement.

clc; clear; close all;
addProjectPaths();

%% ---- Problem setup (mirror runCameraOptimiser: 7-cam UGV, CF3) ----
numCams          = 7;
volume           = [-4 4; -4 4; 0 4];
cameraLowerBounds = [-5 -4.5 0  -pi -pi -pi];
cameraUpperBounds = [ 5  4.5 4.8 pi  pi  pi];
targetMode       = 1;
spacing          = 1;
UGV_maxHeight    = 0.5;
UGV_zSpacing     = 0.25;

volume(3,:)   = [0, UGV_maxHeight];
targetSpacing = [spacing, spacing, min(UGV_zSpacing, UGV_maxHeight)];

specs = setupHardwareSpecs(numCams);
specs.TargetType = 2;
specs.TargetMode = targetMode;
specs.Target     = generateTargetSpace(volume, targetMode, targetSpacing);
specs.NumPoints  = size(specs.Target,1);
specs.spacing    = spacing;
specs.spacingZ   = UGV_zSpacing;
specs.SectionCentres = generateSectionCentres(numCams, volume);
specs = setupCostParams(specs);

Target     = specs.Target;
numPoints  = specs.NumPoints;
resolution = specs.Resolution;
adj = specs.PreComputed.adjacentSurfaces;
du  = specs.PreComputed.du;
dv  = specs.PreComputed.dv;
pen = specs.PreComputed.penaltyUncertainty;
w2  = specs.PreComputed.w2;

VarMin = repmat(cameraLowerBounds,1,numCams);
VarMax = repmat(cameraUpperBounds,1,numCams);

fprintf('Target points: %d,  cameras: %d\n\n', numPoints, numCams);

%% ---- Equivalence check across random camera layouts ----
nTest  = 25;
maxErr = 0;
nMism  = 0;
worst  = struct('err',0,'visNew',NaN,'visOld',NaN,'t',NaN,'p',NaN);
rng(1);
for t = 1:nTest
    chrom = VarMin + rand(1,6*numCams).*(VarMax - VarMin);
    [cameras, CamCenters] = setupCameras(chrom, numCams, resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    [U, V, visMask] = projectAllPoints(cameras, Target, CamCenters, resolution);

    % Ground-truth visibility from per-point project() (matches original)
    visOldAll = false(numPoints, numCams);
    for p = 1:numPoints
        for i = 1:numCams
            uv = cameras{i}.project(Target(p,:));
            visOldAll(p,i) = uv(1) >= 1 && uv(1) <= resolution(1) && ...
                             uv(2) >= 1 && uv(2) <= resolution(2);
        end
    end

    for p = 1:numPoints
        aNew = computePointUncertainty(Target(p,:), cameras, CamCenters, numCams, ...
            adj, du, dv, pen, w2, resolution, U(p,:), V(p,:), visMask(p,:));
        aOld = computePointUncertainty(Target(p,:), cameras, CamCenters, numCams, ...
            adj, du, dv, pen, w2, resolution);   % original per-point path
        e = abs(aNew - aOld);
        maxErr = max(maxErr, e);
        if e > 1e-9
            nMism = nMism + 1;
            if e > worst.err
                worst = struct('err',e,'visNew',sum(visMask(p,:)), ...
                    'visOld',sum(visOldAll(p,:)),'t',t,'p',p);
            end
        end
    end
end
fprintf('Equivalence: max |new - old| = %.3e over %d layouts x %d points\n', ...
    maxErr, nTest, numPoints);
fprintf('Mismatching points (>1e-9): %d of %d\n', nMism, nTest*numPoints);
if nMism > 0
    fprintf('  worst: err=%.3e  visible cams new=%d old=%d  (layout %d, point %d)\n', ...
        worst.err, worst.visNew, worst.visOld, worst.t, worst.p);
end
if maxErr < 1e-9
    fprintf('  PASS -- batched projection matches per-point projection.\n\n');
else
    fprintf('  ** MISMATCH -- investigate before using. **\n\n');
end

%% ---- Timing: per-point project vs batched projectAllPoints ----
chrom = VarMin + rand(1,6*numCams).*(VarMax - VarMin);
[cameras, CamCenters] = setupCameras(chrom, numCams, resolution, ...
    specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

% warm-up (JIT)
projectAllPoints(cameras, Target, CamCenters, resolution);
cameras{1}.project(Target(1,:));

nRep = 50;
tOld = 0; tNew = 0;
for r = 1:nRep
    tic;
    for p = 1:numPoints
        for i = 1:numCams
            cameras{i}.project(Target(p,:));
        end
    end
    tOld = tOld + toc;

    tic;
    [~,~,~] = projectAllPoints(cameras, Target, CamCenters, resolution);
    tNew = tNew + toc;
end
fprintf('Projection time over %d reps:\n', nRep);
fprintf('  per-point  : %.4f s\n', tOld);
fprintf('  batched    : %.4f s\n', tNew);
fprintf('  speedup    : %.1fx (projection stage only)\n\n', tOld / tNew);

%% ---- Sanity: full resUncertainty still runs ----
u = resUncertainty(specs, cameras, CamCenters);
fprintf('resUncertainty() returns %.6f (finite: %d)\n', u, isfinite(u));
