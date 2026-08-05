function specs = setupHardwareSpecs(numCams)
% Taken from the OptiTrack datasheet.
    specs.Cams = numCams;
    specs.Resolution = [1280 1024];
    specs.PixelSize = 4.57e-6;  % Square pixel size [m]; back-solved
                                % from datasheet (f=5.5mm, HFOV=56°,
                                % W=1280).
    specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
    specs.Focal = 0.0055;       % Narrow-angle focal length [m]
    specs.FocalWide = 0.0035;   % Wide-angle focal length [m]
    specs.Range = 16;           % Narrow-angle effective range [m] for
                                % passive markers (800 exposure, gain 6,
                                % lowest f-stop)
    specs.RangeWide = 9;        % Wide-angle effective range [m]
end