function angles = pairwiseBaselineAngles(chrom, specs)
%PAIRWISEBASELINEANGLES  Convergence angles of every visible camera pair.
%
%   angles = pairwiseBaselineAngles(chrom, specs)
%
%   For each target point in specs.Target, finds every camera that can see
%   it and computes the convergence (baseline) angle at the point for every
%   pair of those cameras. Returns one concatenated column of angles in
%   degrees over all (point, pair) combinations.
%
%   The triangulability band of Section III-B accepts a pair only when this
%   angle lies in [specs.PreComputed.minTriangAngle,
%   specs.PreComputed.maxTriangAngle] — so the fraction of this vector
%   inside that band is the "well-conditioned pairs" statistic.
%
%   findVisibleCameras returns unit vectors pointing from each camera
%   towards the point. The angle between two such vectors equals the angle
%   subtended at the point between the two camera lines of sight, so
%   acosd(dot(.,.)) is the triangulation angle directly.
%
%   Lifted out of plotBaselineAngles_GAvsOptiTrack.m so the figure and the
%   reported numbers cannot drift apart.
%
%   See also: plotBaselineAngles_GAvsOptiTrack, findVisibleCameras,
%             perTargetCoverage.

    numCams = specs.Cams;
    [cameras, camCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, ...
        specs.PrincipalPoint, specs.PixelSize);

    resolution         = specs.Resolution;
    TargetSpace        = specs.Target;
    maxCameraRange     = specs.PreComputed.maxCameraRange;
    maxCameraRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide          = specs.FocalWide;

    nPts = size(TargetSpace, 1);

    % Worst-case upper bound on pairs so the pool can be preallocated.
    maxPairsPerPoint = numCams * (numCams - 1) / 2;
    pool = nan(nPts * maxPairsPerPoint, 1);
    head = 0;

    for pt = 1:nPts
        point = TargetSpace(pt, :);
        [visCams, viewVecs] = findVisibleCameras(point, cameras, camCenters, ...
            numCams, resolution, maxCameraRange, maxCameraRangeWide, focalWide);
        nv = numel(visCams);
        if nv < 2
            continue;
        end
        for i = 1:(nv-1)
            for j = (i+1):nv
                c = dot(viewVecs(:,i), viewVecs(:,j));
                c = max(min(c, 1), -1);
                head = head + 1;
                pool(head) = acosd(c);
            end
        end
    end

    angles = pool(1:head);
end
