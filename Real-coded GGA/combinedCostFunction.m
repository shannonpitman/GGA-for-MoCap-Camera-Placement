function totalCost = combinedCostFunction(cameraChromosome, specs)
%Combined cost function using weighted approach between resolution 
%uncertainty and dynamic occlusion
    
    %Adjustable Weight parameters
    w_uncertainty = specs.WeightUncertainty; 
    w_occlusion = specs.WeightOcclusion; 

    % Initialize cameras
    numCams = specs.Cams;
    resolution = specs.Resolution;
    pixelSize = specs.PixelSize; 
    focalLength = specs.Focal; 
    PrincipalPoint = specs.PrincipalPoint; 
    Npix = specs.npix;
    [cameras, CamCenters] = setupCameras(cameraChromosome, numCams, resolution, pixelSize, focalLength, PrincipalPoint,Npix);
    % Calculate individual cost components
    uncertaintyCost = resUncertainty(specs, cameras, CamCenters);
    occlusionCost = dynamicOcclusion(specs, cameras, CamCenters);
    
    % Combined weighted cost
    totalCost = w_uncertainty * uncertaintyCost + w_occlusion * occlusionCost;
end