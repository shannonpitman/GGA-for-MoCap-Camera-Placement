function plotCameraFOV(camera, camCentre, effectiveRange, varargin)
% PLOTCAMERAFOV  Draw a camera's true field-of-view frustum as a solid 3D pyramid.
%
% plotCameraFOV(camera, camCentre, effectiveRange) plots the pyramid that
% extends from the camera's perspective centre out to its effective range,
% with the apex at camCentre and the base sized by the camera's intrinsic
% matrix and image resolution. The drawn frustum therefore matches what
% findVisibleCameras / computePointUncertainty actually use as the camera's
% visibility region.
%
% The pyramid is built from FIVE translucent patches (4 triangular side
% faces + 1 quadrilateral far face) plus the 4 apex-to-corner edge lines
% drawn on top, so the volume reads as unambiguously 3D from any view
% angle. A small wedge representing the camera body is drawn at the apex
% (use Corke's RVC3 plot_camera if available, falling back to a hand-
% drawn marker).
%
% Optional name/value pairs:
%   'Color'         - RGB triplet (default [0 0.4 0.7])
%   'EdgeAlpha'     - 0..1 (default 0.85)
%   'FaceAlpha'     - 0..1 (default 0.22) — alpha of side + far faces
%   'LineWidth'     - default 1.0
%   'Label'         - text label drawn at the apex (default '')
%   'ShowBody'      - true to draw a small camera-body wedge at the apex.
%                     Default true.
%   'BodyScale'     - scale factor for the body wedge (m). Default 0.25.
%   'ShowOpticalAxis' - true to draw a dashed line along the optical
%                       axis from apex to far-plane centre. Default true.
%   'VisualRange'   - if non-empty, cap the drawn pyramid length at
%                     this value (m). Useful when effectiveRange is much
%                     larger than the scene (e.g. 16 m range in an 8 m
%                     room). The actual effectiveRange is unchanged for
%                     cost calculations — this only shortens what is
%                     DRAWN. Default [] (no cap).
%
% camera must be a Peter-Corke CentralCamera built with the focal length
% and resolution actually used by setupCameras (so K is meaningful).

    p = inputParser;
    addParameter(p, 'Color',           [0 0.4 0.7]);
    addParameter(p, 'EdgeAlpha',       0.85);
    addParameter(p, 'FaceAlpha',       0.22);
    addParameter(p, 'LineWidth',       1.0);
    addParameter(p, 'Label',           '');
    addParameter(p, 'ShowBody',        true,    @islogical);
    addParameter(p, 'BodyScale',       0.25,    @isnumeric);
    addParameter(p, 'ShowOpticalAxis', true,    @islogical);
    addParameter(p, 'VisualRange',     [],      @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});
    opts = p.Results;

    % Cap the drawn pyramid length without affecting the caller's notion
    % of effectiveRange (which is used for visibility cost calculations).
    drawRange = effectiveRange;
    if ~isempty(opts.VisualRange) && opts.VisualRange > 0
        drawRange = min(drawRange, opts.VisualRange);
    end

    camCentre = camCentre(:);                         % 3x1

    % Image-plane corners in pixels (1..W, 1..H), ordered CCW so adjacent
    % columns of cornersPix are adjacent corners of the image plane.
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

    % Normalise and scale to the (possibly capped) drawing range.
    rayLens    = vecnorm(raysWorld, 2, 1);
    rayUnit    = raysWorld ./ rayLens;
    farCorners = camCentre + rayUnit * drawRange;     % 3x4 far-plane corners

    % Optical axis (camera +Z in world frame) — used for body + axis line.
    opticAxisWorld = R(:, 3);                         % 3x1 unit
    farCentre      = camCentre + opticAxisWorld * drawRange;

    hold('on');

    % --- (1) Side faces as triangular patches (apex + two adj. corners) ---
    %     This is what makes the frustum read as a solid pyramid.
    for c = 1:4
        c2 = mod(c, 4) + 1;     % next corner index (wrap around)
        triX = [camCentre(1), farCorners(1,c), farCorners(1,c2)];
        triY = [camCentre(2), farCorners(2,c), farCorners(2,c2)];
        triZ = [camCentre(3), farCorners(3,c), farCorners(3,c2)];
        patch(triX, triY, triZ, opts.Color, ...
            'FaceAlpha', opts.FaceAlpha, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    % --- (2) Far-plane quad as patch (slightly stronger alpha) ---
    farX = farCorners(1, :);
    farY = farCorners(2, :);
    farZ = farCorners(3, :);
    patch(farX, farY, farZ, opts.Color, ...
        'FaceAlpha', min(opts.FaceAlpha * 1.5, 1), ...
        'EdgeColor', opts.Color, ...
        'EdgeAlpha', opts.EdgeAlpha, ...
        'LineWidth', opts.LineWidth, ...
        'HandleVisibility', 'off');

    % --- (3) Apex-to-corner edges (drawn on top of patches for clarity) ---
    for c = 1:4
        plot3([camCentre(1) farCorners(1,c)], ...
              [camCentre(2) farCorners(2,c)], ...
              [camCentre(3) farCorners(3,c)], ...
              '-', 'Color', [opts.Color opts.EdgeAlpha], ...
              'LineWidth', opts.LineWidth, ...
              'HandleVisibility', 'off');
    end

    % --- (4) Optical axis (dashed line from apex to far-plane centre) ---
    if opts.ShowOpticalAxis
        plot3([camCentre(1) farCentre(1)], ...
              [camCentre(2) farCentre(2)], ...
              [camCentre(3) farCentre(3)], ...
              ':', 'Color', [opts.Color 0.55], ...
              'LineWidth', max(opts.LineWidth - 0.2, 0.6), ...
              'HandleVisibility', 'off');
    end

    % --- (5) Camera body at apex ---
    %     A small hand-drawn wedge oriented along the optical axis. We
    %     deliberately do NOT call camera.plot_camera() from RVC3 here:
    %     that method has side effects on axis limits, view angle, and
    %     axis-equal mode, which would clobber the parent subplot.
    if opts.ShowBody
        % Big filled marker at the apex (the lens centre).
        plot3(camCentre(1), camCentre(2), camCentre(3), 'o', ...
            'MarkerSize', 8, 'MarkerFaceColor', opts.Color, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.7, ...
            'HandleVisibility', 'off');
        % Short stub along the optical axis — visually "the lens".
        stub = camCentre + opticAxisWorld * (opts.BodyScale * 0.7);
        plot3([camCentre(1) stub(1)], ...
              [camCentre(2) stub(2)], ...
              [camCentre(3) stub(3)], ...
              '-', 'Color', opts.Color, ...
              'LineWidth', opts.LineWidth + 0.6, ...
              'HandleVisibility', 'off');
    end

    % --- (6) Label drawn last so it's not occluded by patches ---
    if ~isempty(opts.Label)
        text(camCentre(1), camCentre(2), camCentre(3), ...
            ['  ' opts.Label], 'Color', opts.Color, ...
            'FontSize', 9, 'FontWeight', 'bold');
    end
end
