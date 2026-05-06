function totalCost = combinedCostFunction(cameraChromosome, specs)
%Combined cost function using weighted approach between resolution 
%uncertainty and dynamic occlusion
    
    %Adjustable Weight parameters
    w_uncertainty = specs.WeightUncertainty; 
    w_occlusion = specs.WeightOcclusion; 

    % Initialize cameras
    numCams = specs.Cams;
    resolution = specs.Resolution;
    % pixelSize = specs.PixelSize; 
    focalLength = specs.Focal; 
    focalLengthWide = specs.FocalWide;
    PrincipalPoint = specs.PrincipalPoint; 
    % Npix = specs.npix;
    uncertNorm = specs.PreComputed.uncertNorm;
    occlNorm = specs.PreComputed.occlNorm;

    [cameras, CamCenters] = setupCameras(cameraChromosome, numCams, resolution, focalLength, focalLengthWide, PrincipalPoint);
    
    % Calculate individual cost components
    uncertaintyCost = resUncertainty(specs, cameras, CamCenters)/ uncertNorm;
    occlusionCost = dynamicOcclusion(specs, cameras, CamCenters)/ occlNorm;
    
    % Combined weighted cost
    totalCost = w_uncertainty*uncertaintyCost + w_occlusion*occlusionCost;
end