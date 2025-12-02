function uncertainty = computePointUncertainty(point, cameras, cameraCentres, numCams, adjacentSurfaces, du,dv, penaltyUncertainty, w2, Resolution)
    uv = zeros(numCams,2);
    planes = cell(1, numCams);
    visIdx = false(numCams,1); %Logical operator to see if camera lies within bounds of image plane and is in front of the camera
    
    % Vectorized visibility check for all cameras
    for i = 1:numCams
        uv(i,:) = cameras{i}.project(point); %uv projected coordinates
        u = uv(1);
        v = uv(2);
        if (u >= 1 && u <= Resolution(1) && v >= 1 && v <= Resolution(2))
            worldPoints = quantToWorld(cameras{i}, u, v, du, dv, cameraCentres(:,i));
            planes{i} = buildPyramidSurf(cameraCentres(:,i), worldPoints, adjacentSurfaces);
            visIdx(i) = true;
        end
    end
    
    visibleIdx = find(visIdx);
    numVisible = length(visibleIdx);
    
    %early exit if triangulation is not possible
    if numVisible == 0 
        uncertainty = penaltyUncertainty; %point is not visible in any cameras FOV
        return;
    elseif numVisible== 1 %if at least one camera sees the point
        uncertainty = 0.5*penaltyUncertainty;
        return;
    end

    vertices = calcVertices(numVisible, visibleIdx, adjacentSurfaces, planes);
    unique_vertices = unique(round(vertices,6),'rows', 'stable');

    V_centered = unique_vertices-point; %centre the data at the point (zero mean) -> standardize the data
    C0v_vert = cov(V_centered); %covariance matrix
    
    Eigs = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
    uncertainty = sum(sqrt(abs(Eigs))); % uncertainty based on the trace of the cov matrix 

    % Old code that found error volume 
    % [A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix] -> eigenvectors: new axes which our data lies on, eigenvalues = variances along these new axes
    % 
    % 
    % transformed_v = (A*V_centered.').';%rotated points 
    % 
    % 
    % 
    % %Initial guess
    % axes0 = max(abs(transformed_v), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 
    % 
    % fun = @(axes) optimiseEllipsoid(axes, transformed_v.', axes0, w2);
    % options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
    % [axes_opt, ~] = fminsearch(fun, axes0, options);
    % 
    % %optimised axes [m]
    % a = abs(axes_opt(1));
    % b = abs(axes_opt(2));
    % c = abs(axes_opt(3));
    % 
    % uncertainty = 4/3*pi*a*b*c;%m^3
end

