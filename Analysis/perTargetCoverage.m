function [cov, stats] = perTargetCoverage(chromosome, specs)
% PERTARGETCOVERAGE  Per-point camera-visibility count for a configuration.
%
%   [cov, stats] = perTargetCoverage(chromosome, specs)
%
%   Standalone, callable version of the local perPointVisibility helper
%   inside plotHeatmap_GAvsOptiTrack.m. Counts, for every target point in
%   specs.Target, how many cameras can actually SEE it. A camera sees a
%   point iff ALL THREE hold:
%       (a) FOV check     - perspective projection lands inside [1,W]x[1,H]
%       (b) Range check   - distance <= effective range (narrow/wide lens)
%       (c) In-front check - depth in camera frame is positive (z_cam > 0)
%
%   OUTPUTS
%     cov   - nPts x 1 vector, cameras visible per target point
%     stats - struct with fields:
%               .avg        mean cameras per point
%               .zeroPct    %% of points seen by 0 cameras
%               .onePct     %% of points seen by exactly 1 camera
%               .twoPlusPct %% of points seen by >= 2 cameras (triangulable)
%               .numCams    camera count
%               .nPts       number of target points
%
%   See also: plotHeatmap_GAvsOptiTrack, surveyGACoverage, findVisibleCameras

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

    % Shared FOV + in-front + range test (same as both cost functions)
    visMask = projectVisibilityOcclusion(cameras, T, camCenters, resolution, ...
        maxRange, maxRangeWide, focalWide);
    cov = sum(visMask, 2);

    if nargout > 1
        stats.avg        = mean(cov);
        stats.zeroPct    = 100 * sum(cov == 0) / nPts;
        stats.onePct     = 100 * sum(cov == 1) / nPts;
        stats.twoPlusPct = 100 * sum(cov >= 2) / nPts;
        stats.numCams    = numCams;
        stats.nPts       = nPts;
    end
end
