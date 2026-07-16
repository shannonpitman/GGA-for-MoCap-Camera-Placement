function [visMask, viewUnit] = projectVisibilityOcclusion(cameras, TargetSpace, camCentres, resolution, maxRange, maxRangeWide, FocalWide)
%PROJECTVISIBILITYOCCLUSION  Batched visibility + view vectors for the
%   dynamic-occlusion cost. Vectorised replacement for the per-point
%   projection/range work inside findVisibleCameras.
%
%   [visMask, viewUnit] = projectVisibilityOcclusion(cameras, TargetSpace, ...
%                             camCentres, resolution, maxRange, maxRangeWide, FocalWide)
%
%   A (point, camera) pair is visible when the point projects inside the
%   image AND is in front of the camera (via projectAllPoints) AND lies
%   within the camera's effective range (wide-angle lenses use the shorter
%   maxRangeWide, narrow-angle use maxRange).
%
%   Outputs:
%     visMask  - N x numCams logical visibility
%     viewUnit - N x numCams x 3 unit vectors from camera centre to point
%                (zero where dist == 0; only used where visMask is true)

    % Image-plane + in-front visibility (shared with the uncertainty path).
    [~, ~, inImage] = projectAllPoints(cameras, TargetSpace, camCentres, resolution);

    numCams = numel(cameras);
    N       = size(TargetSpace, 1);
    visMask  = false(N, numCams);
    viewUnit = zeros(N, numCams, 3);

    for i = 1:numCams
        c    = camCentres(:, i).';          % 1 x 3 camera centre
        dvec = TargetSpace - c;             % N x 3  point - camera centre
        dist = vecnorm(dvec, 2, 2);         % N x 1

        if cameras{i}.f == FocalWide
            effRange = maxRangeWide;
        else
            effRange = maxRange;
        end
        inRange = dist <= effRange & dist > 0;

        visMask(:, i) = inImage(:, i) & inRange;

        unit = dvec ./ dist;               % N x 3 (guard dist==0 below)
        unit(~isfinite(unit)) = 0;
        viewUnit(:, i, 1) = unit(:, 1);
        viewUnit(:, i, 2) = unit(:, 2);
        viewUnit(:, i, 3) = unit(:, 3);
    end
end
