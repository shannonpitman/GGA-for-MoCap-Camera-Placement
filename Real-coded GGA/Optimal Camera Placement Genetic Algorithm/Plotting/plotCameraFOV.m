function plotCameraFOV(camera, camCentre, effectiveRange, varargin)
% PLOTCAMERAFOV  Draw a camera's true field-of-view frustum.
%
% plotCameraFOV(camera, camCentre, effectiveRange) plots the pyramid that
% extends from the camera's perspective centre out to its effective range,
% with the apex at camCentre and the base sized by the camera's intrinsic
% matrix and image resolution. The drawn frustum therefore matches what
% findVisibleCameras / computePointUncertainty actually use as the camera's
% visibility region.
%
% Optional name/value pairs:
%   'Color'      - RGB triplet (default [0 0.4 0.7])
%   'EdgeAlpha'  - 0..1 (default 0.7)
%   'FaceAlpha'  - 0..1 (default 0.08)
%   'LineWidth'  - default 0.8
%   'Label'      - text label drawn at the apex (default '')
%
% camera must be a Peter-Corke CentralCamera built with the focal length
% and resolution actually used by setupCameras (so K is meaningful).

    p = inputParser;
    addParameter(p, 'Color',     [0 0.4 0.7]);
    addParameter(p, 'EdgeAlpha', 0.7);
    addParameter(p, 'FaceAlpha', 0.08);
    addParameter(p, 'LineWidth', 0.8);
    addParameter(p, 'Label',     '');
    parse(p, varargin{:});
    opts = p.Results;

    camCentre = camCentre(:);                         % 3x1

    % Image-plane corners in pixels (1..W, 1..H).
    res = camera.npix;                                 % [W H] from CentralCamera
    if isempty(res)
        error('plotCameraFOV:NoResolution', ...
            'Camera object must be constructed with a resolution.');
    end
    W = res(1); H = res(2);
    cornersPix = [1  W  W  1;
                  1  1  H  H;
                  1  1  1  1];                         % 3x4 homogeneous pixels

    % Pixel -> ray in camera frame, then -> world frame.
    K = camera.K;                                     % 3x3 intrinsic matrix
    raysCam   = K \ cornersPix;                       % 3x4 directions, camera frame
    R         = camera.T.rotm;                        % camera-to-world rotation
    raysWorld = R * raysCam;                          % 3x4 world directions

    % Normalise and scale to the effective range.
    rayLens   = vecnorm(raysWorld, 2, 1);
    rayUnit   = raysWorld ./ rayLens;
    farCorners = camCentre + rayUnit * effectiveRange; % 3x4 far-plane corners

    % --- Draw apex point + label ---
    hold('on');
    plot3(camCentre(1), camCentre(2), camCentre(3), '.', ...
        'MarkerSize', 14, 'Color', opts.Color);
    if ~isempty(opts.Label)
        text(camCentre(1), camCentre(2), camCentre(3), ...
            ['  ' opts.Label], 'Color', opts.Color, 'FontSize', 8);
    end

    % --- Draw four edges from apex to far corners ---
    for c = 1:4
        plot3([camCentre(1) farCorners(1,c)], ...
              [camCentre(2) farCorners(2,c)], ...
              [camCentre(3) farCorners(3,c)], ...
              '-', 'Color', [opts.Color opts.EdgeAlpha], ...
              'LineWidth', opts.LineWidth);
    end

    % --- Draw far-plane quad as a translucent patch ---
    farX = [farCorners(1,:) farCorners(1,1)];
    farY = [farCorners(2,:) farCorners(2,1)];
    farZ = [farCorners(3,:) farCorners(3,1)];
    patch(farX, farY, farZ, opts.Color, ...
        'FaceAlpha', opts.FaceAlpha, ...
        'EdgeColor', opts.Color, ...
        'EdgeAlpha', opts.EdgeAlpha, ...
        'LineWidth', opts.LineWidth);
end
