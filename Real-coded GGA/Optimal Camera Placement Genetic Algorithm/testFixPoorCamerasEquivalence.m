%% testFixPoorCamerasEquivalence.m
% Verifies the vectorised camera-coverage counts (cameraCoverageCounts, used
% by the new fixPoorCameras) match the original object-based, per-point
% method (setupCameras + findVisibleCameras) EXACTLY, and times the speedup.
%
% Coverage counts are integers, so a correct vectorisation must match to the
% exact count -- any intrinsics/rotation mismatch shows up immediately.

clc; clear; close all;
addProjectPaths();

%% ---- Problem setup (7-cam UGV, CF3) ----
numCams          = 7;
volume           = [-4 4; -4 4; 0 4];
cameraLowerBounds = [-5 -4.5 0  -pi -pi -pi];
cameraUpperBounds = [ 5  4.5 4.8 pi  pi  pi];
spacing          = 1;
UGV_maxHeight    = 0.5;
UGV_zSpacing     = 0.25;

volume(3,:)   = [0, UGV_maxHeight];
targetSpacing = [spacing, spacing, min(UGV_zSpacing, UGV_maxHeight)];

specs = setupHardwareSpecs(numCams);
specs.TargetType = 2;
specs.TargetMode = 1;
specs.Target     = generateTargetSpace(volume, 1, targetSpacing);
specs.NumPoints  = size(specs.Target,1);
specs.spacing    = spacing;
specs.spacingZ   = UGV_zSpacing;
specs.SectionCentres = generateSectionCentres(numCams, volume);
specs = setupCostParams(specs);

Target     = specs.Target;
numPoints  = specs.NumPoints;
resolution = specs.Resolution;
maxR   = specs.PreComputed.maxCameraRange;
maxRW  = specs.PreComputed.maxCameraRangeWide;
fW     = specs.FocalWide;

VarMin = repmat(cameraLowerBounds,1,numCams);
VarMax = repmat(cameraUpperBounds,1,numCams);

fprintf('Target points: %d, cameras: %d\n\n', numPoints, numCams);

%% ---- Equivalence: coverage counts (new vs original) ----
nTest = 50;
maxErr = 0;
nLayoutMism = 0;
rng(1);
for t = 1:nTest
    chrom = VarMin + rand(1,6*numCams).*(VarMax - VarMin);

    % NEW: vectorised, object-free
    covNew = cameraCoverageCounts(chrom, specs);

    % OLD: setupCameras + per-point findVisibleCameras
    [cameras, camCenters] = setupCameras(chrom, numCams, resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    covOld = zeros(numCams,1);
    for p = 1:numPoints
        [vis,~] = findVisibleCameras(Target(p,:), cameras, camCenters, numCams, ...
            resolution, maxR, maxRW, fW);
        for k = 1:numel(vis)
            covOld(vis(k)) = covOld(vis(k)) + 1;
        end
    end

    e = max(abs(covNew - covOld));
    maxErr = max(maxErr, e);
    if e > 0, nLayoutMism = nLayoutMism + 1; end
end
fprintf('Coverage equivalence: max |new - old| = %d count(s), layouts differing = %d of %d\n', ...
    maxErr, nLayoutMism, nTest);
if maxErr == 0
    fprintf('  PASS -- vectorised coverage matches object-based coverage exactly.\n\n');
else
    fprintf('  ** MISMATCH -- investigate before using. **\n\n');
end

%% ---- Equivalence: full fixPoorCameras output ----
% (reorientation logic is unchanged, so identical coverage => identical y)
maxYErr = 0;
rng(7);
for t = 1:nTest
    chrom = VarMin + rand(1,6*numCams).*(VarMax - VarMin);
    yNew = fixPoorCameras(chrom, specs, 0.05);
    yRef = fixPoorCamerasOld(chrom, specs, 0.05);   % local copy below
    maxYErr = max(maxYErr, max(abs(yNew - yRef)));
end
fprintf('fixPoorCameras output: max |new - old| = %.3e\n', maxYErr);
if maxYErr < 1e-9
    fprintf('  PASS -- identical corrected chromosomes.\n\n');
else
    fprintf('  ** MISMATCH -- investigate. **\n\n');
end

%% ---- Timing ----
nRep = 200;
chrom = VarMin + rand(1,6*numCams).*(VarMax - VarMin);
fixPoorCameras(chrom, specs, 0.05);          % warm-up
fixPoorCamerasOld(chrom, specs, 0.05);

tic; for r = 1:nRep, fixPoorCameras(chrom, specs, 0.05);    end; tNew = toc;
tic; for r = 1:nRep, fixPoorCamerasOld(chrom, specs, 0.05); end; tOld = toc;
fprintf('fixPoorCameras timing over %d calls:\n', nRep);
fprintf('  original (objects + per-point): %.4f s\n', tOld);
fprintf('  vectorised                    : %.4f s\n', tNew);
fprintf('  speedup                       : %.1fx\n', tOld / tNew);

%% ---- Local copy of the ORIGINAL fixPoorCameras (reference) ----
function y = fixPoorCamerasOld(x, specs, coverageThreshold)
    numCams = specs.Cams;
    TargetSpace = specs.Target;
    numPoints = size(TargetSpace, 1);
    minPointsRequired = ceil(coverageThreshold * numPoints);
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;
    maxRange = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    [cameras, camCenters] = setupCameras(x, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint, specs.PixelSize);
    cameraCoverage = zeros(numCams, 1);
    for p = 1:numPoints
        point = TargetSpace(p, :);
        [visibleCams, ~] = findVisibleCameras(point, cameras, camCenters, numCams, resolution, maxRange, maxRangeWide, focalLengthWide);
        for k = 1:length(visibleCams)
            cameraCoverage(visibleCams(k)) = cameraCoverage(visibleCams(k)) + 1;
        end
    end
    y = x;
    for c = 1:numCams
        if cameraCoverage(c) < minPointsRequired
            chromStart = (c-1)*6 + 1;
            chromEnd = c*6;
            camPos = x(chromStart:chromStart+2);
            distances = vecnorm(specs.SectionCentres - camPos, 2, 2);
            [~, closestIdx] = min(distances);
            targetPoint = specs.SectionCentres(closestIdx, :);
            directionVec = targetPoint - camPos;
            directionUnit = directionVec / norm(directionVec);
            Z_axis = [0, 0, 1];
            rotAng = acos(dot(Z_axis, directionUnit));
            rotAxis = cross(Z_axis, directionUnit);
            if norm(rotAxis) > 0
                rotAxis = rotAxis / norm(rotAxis);
                q = quaternion([cos(rotAng/2), rotAxis * sin(rotAng/2)]);
                q = normalize(q);
                q2Eul = euler(q, 'XYZ', 'point');
                y(chromEnd-2:chromEnd) = q2Eul;
            end
        end
    end
end
