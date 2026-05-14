function cost = resUncertaintyCost(cameraChromosome, specs)
%For cost function type = 1
    numCams = specs.Cams;
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;

    [cameras, CamCenters] = setupCameras(cameraChromosome, numCams, ...
        resolution, focalLength, focalLengthWide, PrincipalPoint, specs.PixelSize);

    % Call the original resolution uncertainty function
    cost = resUncertainty(specs, cameras, CamCenters);
end