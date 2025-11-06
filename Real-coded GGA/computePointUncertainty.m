function [uncertainty, colourPoint] = computePointUncertainty(point, cameras, numCams, resolution, adjacentSurfaces, du,dv, penaltyUncertainty, w2)
    uv = zeros(numCams,2);
    planes = cell(1, numCams);
    
    for i = 1:numCams
        uv(i,:) = cameras{i}.project(point); %uv projected coordinates 
    end

    isVisible = (uv(:,1) >= 1) & (uv(:,1) <= resolution(1)) & (uv(:,2) >= 1) & (uv(:,2) <= resolution(2));
    visibleIdx = find(isVisible);
    numVisible = length(visibleIdx);
    
    %early exit if trinagulation is not possible
    if numVisible == 0 
        uncertainty = penaltyUncertainty; %point is not visible in any cameras FOV
        colourPoint = [1,0,0]; %red 
        return;
    elseif numVisible== 1 %if at least one camera sees the point
        uncertainty = 0.5*penaltyUncertainty;
        colourPoint = [0,1,1]; %cyan
        return;
    end

   for i = visibleIdx'
       u = uv(i,1);
       v = uv(i,2);
       cameraCentres = cameras{i}.center().'; %world location of Camera center
       worldPoints = quantToWorld(cameras{i}, u,v, du, dv, cameraCentres);
       planes{i} = buildPyramidSurf(cameraCentres, worldPoints, adjacentSurfaces);
   end

    vertices = calcVertices(numVisible, visibleIdx, adjacentSurfaces, planes);
    unique_vertices = unique(round(vertices,6),'rows', 'stable');

    V_ = unique_vertices-point; %centre the points at the mean 
    C0v_vert = cov(V_); %covariance matrix
    [A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
    
    
    transformed_v = (A*V_.').';%rotated points 

    %Initial guess
    axes0 = max(abs(transformed_v), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 

    fun = @(axes) optimiseEllipsoid(axes, transformed_v.', axes0, w2);
    options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
    [axes_opt, ~] = fminsearch(fun, axes0, options);

    %optimised axes [m]
    a = abs(axes_opt(1));
    b = abs(axes_opt(2));
    c = abs(axes_opt(3));

    uncertainty = 4/3*pi*a*b*c;%m^3
    colourPoint = [0,1,0]; %green 
end

