function specs = setupCostParams(specs)

% Pre-compute data for resolution uncertainty
specs.PreComputed.adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
specs.PreComputed.du = 0.5; %pixels
specs.PreComputed.dv = 0.5; %pixels
specs.PreComputed.penaltyUncertainty = 100; % penalty for unseen points
specs.PreComputed.w2 = 0.2; %weight for ellipsoid optimization

% Pre-compute data for dynamic occlusion
specs.PreComputed.minTriangAngle = 40; % min convergence angle [degrees]
specs.PreComputed.maxTriangAngle = 140; % max convergence angle [degrees]
specs.PreComputed.maxCameraRange = specs.Range; %m narrowFOV effective range
specs.PreComputed.maxCameraRangeWide = specs.RangeWide; %m wideFOV effective range

% Pre-compute target space in homogeneous coordinates for batch processing
specs.TargetHomogeneous = [specs.Target'; ones(1, specs.NumPoints)];

%Normalisation constants for combined cost function
specs.PreComputed.uncertNorm = 100;  
specs.PreComputed.occlNorm = 440;  