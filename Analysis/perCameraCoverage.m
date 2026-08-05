function [camCov, stats] = perCameraCoverage(chromosome, specs)
% PERCAMERACOVERAGE  Per-camera grid-coverage for a configuration.
%
%   [camCov, stats] = perCameraCoverage(chromosome, specs)
%
%   Counterpart to perTargetCoverage: instead of asking "how many cameras
%   see each point?", it asks "how many points does each camera see?".
%   Uses the IDENTICAL visibility test (FOV + range + in-front) so the
%   numbers are consistent with perTargetCoverage and with the coverage
%   check inside fixPoorCameras. A camera sees a point iff ALL THREE hold:
%       (a) FOV check      - perspective projection lands inside [1,W]x[1,H]
%       (b) Range check    - distance <= effective range (narrow/wide lens)
%       (c) In-front check - depth in camera frame is positive (z_cam > 0)
%
%   OUTPUTS
%     camCov - numCams x 1 vector, number of grid points visible per camera
%     stats  - struct with fields:
%               .pct       numCams x 1, %% of grid points each camera sees
%               .avgPct    mean per-camera coverage (%%)
%               .minPct    lowest per-camera coverage (%%)
%               .maxPct    highest per-camera coverage (%%)
%               .numCams   camera count
%               .nPts      number of target points
%
%   Use this to justify the fixPoorCameras 5%% repair threshold: compare a
%   working camera's coverage (e.g. avgPct / maxPct) against the floor.
%
%   See also: perTargetCoverage, fixPoorCameras, findVisibleCameras, surveyGACoverage

    numCams         = specs.Cams;
    resolution      = specs.Resolution;
    focalLength     = specs.Focal;
    focalLengthWide = specs.FocalWide;
    principalPoint  = specs.PrincipalPoint;
    T               = specs.Target;
    nPts            = size(T, 1);

    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    [cameras, camCenters] = setupCameras(chromosome, numCams, resolution, ...
        focalLength, focalLengthWide, principalPoint, specs.PixelSize);

    % World-frame optical axis (+Z column of each camera's rotation matrix).
    opticAxes = zeros(3, numCams);
    for c = 1:numCams
        Rcw = cameras{c}.T.rotm;
        opticAxes(:, c) = Rcw(:, 3);
    end

    camCov = zeros(numCams, 1);
    for c = 1:numCams
        for pt = 1:nPts
            point = T(pt, :);

            % (c) In front of camera
            viewVec = point(:) - camCenters(:, c);
            depth   = dot(viewVec, opticAxes(:, c));
            if depth <= 0
                continue;
            end

            % (a) Inside image plane
            uv = cameras{c}.project(point);
            if ~(uv(1) >= 1 && uv(1) <= resolution(1) && ...
                 uv(2) >= 1 && uv(2) <= resolution(2))
                continue;
            end

            % (b) Within effective range
            distance = norm(viewVec);
            if cameras{c}.f == focalWide
                effRange = maxRangeWide;
            else
                effRange = maxRange;
            end
            if distance > effRange || distance <= 0
                continue;
            end

            camCov(c) = camCov(c) + 1;
        end
    end

    if nargout > 1
        stats.pct     = 100 * camCov / nPts;
        stats.avgPct  = mean(stats.pct);
        stats.minPct  = min(stats.pct);
        stats.maxPct  = max(stats.pct);
        stats.numCams = numCams;
        stats.nPts    = nPts;
    end
end
