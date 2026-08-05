function cost = dynamicOcclusionCost(cameraChromosome, specs)
% For cost function type 3

    numCams = specs.Cams;
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;

    [cameras, CamCenters] = setupCameras(cameraChromosome, numCams, ...
        resolution, focalLength, focalLengthWide, PrincipalPoint, specs.PixelSize);

    cost = dynamicOcclusion(specs, cameras, CamCenters);
end