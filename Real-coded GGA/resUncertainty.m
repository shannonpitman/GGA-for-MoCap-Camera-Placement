function errorVolume = resUncertainty(cameraChromosome, specs)
%Computes total uncertainty due to image quantisation over a 3D target
%Space for a given chromosome (camera arrangement) 

numCams = specs.Cams;
w2 = 0.2; %weight for second energy function
%Camera Parameters
resolution = specs.Resolution;
pixelSize = specs.PixelSize;
focalLength = specs.Focal;
PrincipalPoint = specs.PrincipalPoint;

TargetSpace = specs.Target;

cameras = cell(numCams,1);
uv = cell(numCams, 2);
cameraCentre = cell(numCams, 3);
worldPoints = cell(3*numCams,4);
planes = cell(numCams,1);

adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];
du = 0.5; %pixels
dv = 0.5; %pixels
penaltyUncertainty = 100; %high uncertainty for points camera cannot see 

%Compute Camera Transforms
for i = 1:numCams
    chromStartIdx = (i-1)*6+1;
    chromEndIdx = i*6;
    camPositions = cameraChromosome(chromStartIdx: chromStartIdx+2);
    camOrientations = cameraChromosome(chromEndIdx-2: chromEndIdx);

    T = se3(eul2rotm(camOrientations, "XYZ"), camPositions); %camera to world cTw
    cameras{i} = CentralCamera(name="cam"+i,resolution= resolution, pixel= pixelSize, focal= focalLength, pose=T, center = PrincipalPoint); 
end    
    

numPoints = size(TargetSpace,1);
uncertainties = zeros(numPoints,1);

for  p =1:numPoints
    point = TargetSpace(p,:);
    isVisible = false(1,numCams); %logical vector to check wehther camera can see point
    
    for i = 1:numCams
        uv{i} = cameras{i}.project(point); %uv projected coordinates 
        u = uv{i}(1);
        v = uv{i}(2);
        if (u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
            cameraCentre{i} = cameras{i}.center().'; %world location of Camera center
            worldPoints{i} = quantToWorld(cameras{i}, uv{i}, du, dv, cameraCentre{i});
            planes{i} = buildPyramidSurf(cameraCentre{i}, worldPoints{i}, adjacentSurfaces);
            isVisible(i) = true;
        end
    end

    visibleIdx = find(isVisible);
    numVisible = numel(visibleIdx);


    if numVisible ==0 %if there are less than 2 cameras then triangulation is not possible 
        uncertainties(p) = penaltyUncertainty; %point is not visible in any cameras FOV
    elseif numVisible==1 %if at least one camera sees the point
        uncertainties(p) = 0.5*penaltyUncertainty;
    else
        vertices = calcVertices(numVisible, adjacentSurfaces, planes);
        unique_vertices = unique(round(vertices,6),'rows', 'stable');
        % size(unique_vertices);
        V_ = unique_vertices-point; %centre the points at the mean 
        C0v_vert = cov(V_); %covariance matrix
        [A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
        transformed_v = A*V_.';%rotated points 

        %Initial guess
        axes0 = max(abs(transformed_v.'), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 

        fun = @(axes) optimiseEllipsoid(axes, transformed_v, axes0, w2);
        x0 = axes0;
        options = optimset('Display', 'notify', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
        [axes_opt, ~] = fminsearch(fun, x0, options);

        %optimised axes [mm]
        a = axes_opt(1)*1000;
        b = axes_opt(2)*1000;
        c = axes_opt(3)*1000;

        uncertainties(p) = 4/3*pi*a*b*c;%mm^3
    end
end
errorVolume = mean(uncertainties);