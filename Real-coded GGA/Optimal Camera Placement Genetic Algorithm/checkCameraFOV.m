function checkCameraFOV()
%CHECKCAMERAFOV  Diagnostic — does the simulated camera FOV match the
% physical OptiTrack lens specs?
%
% Prints the implied horizontal and vertical FOV for both lens types
% under three regimes:
%   1. Old code path  : pixelSize ignored (CentralCamera default 10 um)
%   2. Wrong PixelSize: 1.4e-6 m as in the original setupHardwareSpecs
%   3. Corrected      : 4.57e-6 m derived from datasheet (current setting)
%
% Datasheet target  : 56° x 46° (narrow), 82° x 60° (wide).
%
% Use after the pixel-size fix to confirm cameras{1}.fov() matches
% physical specs before re-running the GA batch.

    addProjectPaths();

    res = [1280 1024];
    fNarrow = 0.0055;
    fWide   = 0.0035;
    pp      = [res(1)/2, res(2)/2];

    % Identity-pose chromosome for two cams: narrow at index 1, wide at 2
    chrom = [0 0 0  0 0 0,    0 0 0  0 0 0];

    rhoOld     = 1.0e-5;     % Corke default when rho not supplied
    rhoOldSet  = 1.4e-6;     % the pre-fix setupHardwareSpecs value
    rhoFixed   = 4.57e-6;    % current (correct) value

    fprintf('\n==========================================================\n');
    fprintf(' Camera FOV diagnostic — datasheet target: 56°x46° / 82°x60°\n');
    fprintf('==========================================================\n');

    runCase('Old (rho omitted; Corke default 10 µm)', chrom, res, ...
            fNarrow, fWide, pp, rhoOld);
    runCase('Wrong PixelSize (1.4 µm)',             chrom, res, ...
            fNarrow, fWide, pp, rhoOldSet);
    runCase('Corrected (4.57 µm)',                  chrom, res, ...
            fNarrow, fWide, pp, rhoFixed);

    fprintf('\nReminder: existing GA results in Results/ were computed\n');
    fprintf('under the OLD path (~99°/122° FOV). Re-run batchRunGA to\n');
    fprintf('optimise against the corrected camera model.\n\n');
end


function runCase(label, chrom, res, fN, fW, pp, rho)
    fprintf('\n--- %s ---\n', label);

    if isempty(rho) || isnan(rho)
        cams = setupCameras(chrom, 2, res, fN, fW, pp);
    else
        cams = setupCameras(chrom, 2, res, fN, fW, pp, rho);
    end

    [hNarrow, vNarrow] = computeFOV(cams{1}, res);
    [hWide,   vWide]   = computeFOV(cams{2}, res);

    fprintf('  Narrow lens (f=%.2f mm): HFOV %5.1f°  VFOV %5.1f°\n', ...
        cams{1}.f * 1e3, hNarrow, vNarrow);
    fprintf('  Wide   lens (f=%.2f mm): HFOV %5.1f°  VFOV %5.1f°\n', ...
        cams{2}.f * 1e3, hWide, vWide);
end


function [hfov, vfov] = computeFOV(cam, res)
% Compute HFOV and VFOV in degrees from camera intrinsics K.
    K = cam.K;
    fx = K(1,1);
    fy = K(2,2);
    hfov = 2 * atand(res(1) / (2 * fx));
    vfov = 2 * atand(res(2) / (2 * fy));
end
