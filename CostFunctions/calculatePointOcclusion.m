function Q = calculatePointOcclusion(visibleCams, camViewVectors, minAngle, maxAngle)
% Occlusion error for a single target point (Rahimian & Kearney 2017, Eq. 1)
%   E = 360 + u  if no camera sees the point (u = 360)
%   E = 360      if exactly one camera sees it
%   E = Q        otherwise, where Q is the sum of occluder orientations for
%                which the point is not visible in two triangulable views
%
% visibleCams    - indices of cameras with the point in FOV and range
% camViewVectors - 3 x k unit view vectors, column i belongs to visibleCams(i)

    numVisible = length(visibleCams);

    if numVisible < 1
        Q = 720;   % 360 + u, u = 360 (paper, Sec. 3.4)
        return;
    elseif numVisible < 2
        Q = 360;
        return;
    end

    % Horizontal direction of each view vector [deg]. Angles stay in the
    % same column order as camViewVectors so indices always match.
    horiz = camViewVectors(1:2, :);
    viewAngles = mod(atan2d(horiz(2, :), horiz(1, :)), 360);

    % A camera straight above/below the point lies on the occluder's axis
    % and is never occluded (paper, Sec. 3.4).
    onAxis = vecnorm(horiz) < 1e-9;

    Q = calculateOccludedSections(viewAngles, onAxis, camViewVectors, minAngle, maxAngle);
end
