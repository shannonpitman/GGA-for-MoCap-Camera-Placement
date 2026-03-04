function specs = setupHardwareSpecs(numCams)
% Taken from the OptiTrack datasheet
    specs.Cams = numCams;
    specs.Resolution = [1280 1024];
    specs.PixelSize = 1.4e-6; % Square pixel size [m]
    specs.PrincipalPoint = [specs.Resolution(1)/2, specs.Resolution(2)/2];
    specs.Focal = 0.0055;  % Narrow-angle focal length [m], OpenMV = 0.0028
    specs.FocalWide = 0.0035;  % Wide-angle focal length [m]
    specs.Range = 16;  % Narrow-angle effective range [m] for passive markers for 800 exposure, gain of 6 and lowest f-stop
    specs.RangeWide = 9;  % Wide-angle effective range [m]
end