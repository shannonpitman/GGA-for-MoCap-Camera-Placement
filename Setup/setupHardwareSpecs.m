function specs = setupHardwareSpecs(numCams, profile)
%SETUPHARDWARESPECS  Camera intrinsics and ranges for a hardware profile.
%
%   specs = setupHardwareSpecs(numCams)               % 'optitrack'
%   specs = setupHardwareSpecs(numCams, 'openmv_rt1062')
%
%   Odd camera slots use the narrow lens (Focal, Range), even slots the wide
%   lens (FocalWide, RangeWide); see setupCameras. Single-lens hardware sets
%   SingleLens = true with FocalWide = Focal and RangeWide = Range, so every
%   slot gets the same camera.

    if nargin < 2 || isempty(profile)
        profile = 'optitrack';
    end

    specs.Cams = numCams;
    specs.HardwareProfile = profile;

    switch lower(profile)
        case 'optitrack'
            % Taken from the OptiTrack datasheet.
            specs.Resolution = [1280 1024];
            specs.PixelSize  = 4.57e-6; % Square pixel size [m]; back-solved
                                        % from datasheet (f=5.5mm, HFOV=56°,
                                        % W=1280).
            specs.Focal      = 0.0055;  % Narrow-angle focal length [m]
            specs.FocalWide  = 0.0035;  % Wide-angle focal length [m]
            specs.Range      = 16;      % Narrow-angle effective range [m] for
                                        % passive markers (800 exposure, gain 6,
                                        % lowest f-stop)
            specs.RangeWide  = 9;       % Wide-angle effective range [m]

        case 'openmv_rt1062'
            % OpenMV Cam RT1062, OV5640 sensor, stock M12 lens + 850 nm filter.
            % Modelled at the full 2592 x 1944 active array: the tracker
            % reports centroids in sensor pixels whatever framesize/window it
            % reads out (Motion-Capture/MATLAB/core/+rig/fromIntrinsicsFile).
            % f = 2327 px is the mean fx of the cam 1-11 intrinsics
            % (Nov 2026 calibrations, range 2317-2339 px), giving an
            % FOV of 58.2 x 45.3 deg (69.6 deg diagonal). Lens distortion
            % (k1 ~ -0.47, barrel) is not modelled, as for OptiTrack.
            specs.Resolution = [2592 1944];
            specs.PixelSize  = 1.4e-6;          % OV5640 pixel pitch [m]
            specs.Focal      = 2327 * 1.4e-6;   % 3.258 mm
            specs.FocalWide  = specs.Focal;     % single lens
            specs.Range      = 8;               % measured IR marker range [m]
            specs.RangeWide  = specs.Range;
            specs.SingleLens = true;

        case 'lowcost'
            % TODO(Shannon): fill in from the low-cost camera datasheet /
            % measurements. NaN values stop a run before it starts.
            specs.Resolution = [NaN NaN];   % [W H] pixels
            specs.PixelSize  = NaN;         % [m]
            specs.Focal      = NaN;         % narrow-slot focal length [m]
            specs.FocalWide  = NaN;         % wide-slot focal length [m] (= Focal if one lens type)
            specs.Range      = NaN;         % narrow-slot effective range [m]
            specs.RangeWide  = NaN;         % wide-slot effective range [m]

        otherwise
            error('setupHardwareSpecs:UnknownProfile', ...
                'Unknown hardware profile "%s".', profile);
    end

    vals = [specs.Resolution, specs.PixelSize, specs.Focal, specs.FocalWide, ...
            specs.Range, specs.RangeWide];
    if any(isnan(vals))
        error('setupHardwareSpecs:Incomplete', ...
            'Hardware profile "%s" has unfilled values. Edit Setup/setupHardwareSpecs.m.', profile);
    end
    specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
    if ~isfield(specs, 'SingleLens')
        specs.SingleLens = false;
    end
end
