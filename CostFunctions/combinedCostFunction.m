function totalCost = combinedCostFunction(cameraChromosome, specs)
%Combined cost function using weighted approach between resolution
%uncertainty and dynamic occlusion.
%
% Component scaling and weighting live in cf3Terms.m (single source of
% truth): each objective is min-max scaled as (raw - utopia)/norm so that
% utopia -> 0 and nadir -> 1, then weighted and summed.

    % Initialize cameras
    numCams = specs.Cams;
    resolution = specs.Resolution;
    focalLength = specs.Focal;
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint;

    [cameras, CamCenters] = setupCameras(cameraChromosome, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint, specs.PixelSize);

    % Raw cost components
    costUnc = resUncertainty(specs, cameras, CamCenters);
    costOcc = dynamicOcclusion(specs, cameras, CamCenters);

    % Weighted, utopia-shifted combined cost
    totalCost = cf3Terms(costUnc, costOcc, specs);
end
