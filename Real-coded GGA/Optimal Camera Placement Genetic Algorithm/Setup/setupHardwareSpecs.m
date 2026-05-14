function specs = setupHardwareSpecs(numCams)
% Taken from the OptiTrack datasheet.
%
% PixelSize note (corrected 2026-05-14):
%   The OptiTrack datasheet does not report sensor pixel pitch directly;
%   only focal length and FOV are listed. Pixel size therefore has to be
%   back-solved from f, FOV, and image width:
%
%     HFOV = 2 * atan( (W/2) * rho / f )
%     => rho = (f / (W/2)) * tan(HFOV/2)
%
%   With f = 5.5 mm, W = 1280, HFOV = 56°:
%     rho = (5.5e-3 / 640) * tan(28°) = 4.57e-6 m
%
%   The narrow-lens VFOV (46°, H=1024) gives 4.56e-6 m — internally
%   consistent. The wide-lens FOV gives 4.75e-6 m on H, 3.95e-6 m on V
%   (these are slightly off the rectilinear projection, expected for a
%   wide-angle lens with mild distortion). We use 4.57e-6 m as the
%   nominal value because the narrow lens is closer to a textbook
%   rectilinear projection.
%
%   The previous PixelSize = 1.4e-6 m, when ignored as it was here, left
%   CentralCamera using its default rho = 10 um, which gave a SIMULATED
%   FOV of ~99°/122° instead of the physical 56°/82°. Cost-function
%   visibility was therefore over-stated and the GA was optimising
%   against an unphysical camera. setupCameras.m now passes PixelSize
%   through, so this value matters.
    specs.Cams = numCams;
    specs.Resolution = [1280 1024];
    specs.PixelSize = 4.57e-6;  % Square pixel size [m]; back-solved
                                % from datasheet (f=5.5mm, HFOV=56°,
                                % W=1280). See header note above.
    specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
    specs.Focal = 0.0055;       % Narrow-angle focal length [m]
    specs.FocalWide = 0.0035;   % Wide-angle focal length [m]
    specs.Range = 16;           % Narrow-angle effective range [m] for
                                % passive markers (800 exposure, gain 6,
                                % lowest f-stop)
    specs.RangeWide = 9;        % Wide-angle effective range [m]
end