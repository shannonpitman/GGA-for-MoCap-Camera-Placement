function vertices = calcVertices(numVisible, visibleIdx, adj, planes)
% Summary of this function goes here
%   Detailed explanation goes here
    vertices = zeros(16*numVisible, 3); %16 = 4 pairs of adj surface of one cam intersected with 4 surfaces on other pyrsamid 
    idx = 1;
    tol = -1e-9;
    for i = 1:numVisible
        cam_surfaces = planes{visibleIdx(i)}; %index plane array at visible camera 
        
        for j =1:4
            n1 = cam_surfaces(1:3,adj(j,1));
            d1 = cam_surfaces(4, adj(j,1));
            n2 = cam_surfaces(1:3,adj(j,2));
            d2 = cam_surfaces(4,adj(j,2));

            for k = 1:numVisible
                if k == i, continue; 
                end %look at all other pyramids 
                other_surfaces = planes{visibleIdx(k)};
                
                %Loop over each surface from other pyramid
                for l = 1:length(other_surfaces)
                    n3 = other_surfaces(1:3,l);
                    d3 = other_surfaces(4,l);
                    
                    A = [n1.'; n2.'; n3.'];
                    if abs(det(A)) < 1e-9 %check that planes aren't parraell and thus illconditioned -> no intesection 
                         continue;
                    end
    
                    b = [d1; d2; d3];
                    x = A\b; %least squares solve -> finds error region around point its trying to locate
                    
                     if isInsidePlanes(x, cam_surfaces, tol)
                        vertices(idx, :) = x.'; %row x,y,z
                        idx = idx +1;
                     end
                end
            end
        end
    end
    vertices = vertices(1:idx-1,:);
end