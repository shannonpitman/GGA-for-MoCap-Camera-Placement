function vertices = calcVertices(numVisible, visibleIdx, adj, planes)
% Summary of this function goes here
%   Detailed explanation goes here
    maxVertices = 16*nchoosek(numVisible,2); %16 = 4 pairs of adj surface of one cam intersected with 4 surfaces on other pyrsamid 
    vertices = zeros(maxVertices,3);
    idx = 1;
    tol = -1e-9;

    
    %Plane params for cameras that can see the point 
    visiblePlanes = cell(numVisible,1);
    for i = 1:numVisible
        visiblePlanes{i} = planes{visibleIdx(i)};
    end

    camPairs = nchoosek(1:numVisible,2);
    numPairs = size(camPairs,1);

    % Camera Pairs
    for pairIdx = 1:numPairs
        i = camPairs(pairIdx, 1);
        k = camPairs(pairIdx, 2);

        %index plane array at visible camera 
        cam_surfacesI = planes{visibleIdx(i)}; 
        cam_surfacesK = planes{visibleIdx(k)}; 
        
        for j =1:4
            n1 = cam_surfacesI(1:3,adj(j,1));
            d1 = cam_surfacesI(4, adj(j,1));
            n2 = cam_surfacesI(1:3,adj(j,2));
            d2 = cam_surfacesI(4,adj(j,2));
            %Process at all other pyramids 
            n3All = cam_surfacesK(1:3,:);
            d3All = cam_surfacesK(4,:);
            
            for l = 1:size(n3All,2)
                n3 = n3All(:,l);
                d3 = d3All(:,l);
                
                A = [n1.'; n2.'; n3.'];
                   
                if abs(det(A)) < -tol %check that planes aren't parraell and thus illconditioned -> no intesection 
                     continue;
                end
                %least squares solve -> finds error region around point its trying to locate
                b = [d1; d2; d3];
                x = A\b; 
                    
                 if isInsidePlanes(x, visiblePlanes, tol, numVisible)
                    vertices(idx, :) = x.'; %row x,y,z
                    idx = idx +1;
                 end
             end
        end
    end
    vertices = vertices(1:idx-1,:);
end