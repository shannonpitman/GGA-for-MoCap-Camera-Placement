function specs = setupHardwareSpecs(numCams, profile)
%SETUPHARDWARESPECS  Camera intrinsics and ranges for a hardware profile.
%
%   specs = setupHardwareSpecs(numCams)               % 'optitrack'
%   specs = setupHardwareSpecs(numCams, 'lowcost')
%
%   Odd camera slots use the narrow lens (Focal, Range), even slots the wide
%   lens (FocalWide, RangeWide); see setupCameras.

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
end
