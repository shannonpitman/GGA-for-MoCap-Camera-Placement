function totalCost = combinedCostFunction(cameraChromosome, specs)
%Combined cost function using weighted approach between resolution 
%uncertainty and dynamic occlusion
    
    %Adjustable Weight parameters
    w_uncertainty = specs.WeightUncertainty; 
    w_occlusion = specs.WeightOcclusion; 
    
    % Calculate individual cost components
    uncertaintyCost = resUncertainty(cameraChromosome, specs);
    occlusionCost = dynamicOcclusion(cameraChromosome, specs);
    
    % Normalize costs to similar scales before combining
    % (both are already in similar scale - uncertainty in volume, occlusion in degrees)
    % 
    
    % Combined weighted cost
    totalCost = w_uncertainty * uncertaintyCost + w_occlusion * occlusionCost;
end