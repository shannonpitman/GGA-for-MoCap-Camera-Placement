function vertices = calcVertices(numVisible, visibleIdx, adj, planes)
% Summary of this function goes here
%   Detailed explanation goes here
vertices = zeros(16*numVisible, 3); %16 = 4 pairs of adj surface of one cam intersected with 4 surfaces on other pyrsamid 
idx = 1;
tol = -1e-9;
for i = 1:numVisible
    planeIndex_i = visibleIdx(i);
    cam_surfaces = planes{planeIndex_i};
    for j = 1:4
        n1 = cam_surfaces(adj(j,1)).n;
        d1 = cam_surfaces(adj(j,1)).d;
        n2 = cam_surfaces(adj(j,2)).n;
        d2 = cam_surfaces(adj(j,2)).d;
        for k = 1:numVisible
            if k == i, continue; 
            end %look at all other pyramids 
            planeIndex_k = planes{k};
            other_surfaces = planes{planeIndex_k};
            
            %Loop over each surface from other pyramid
            for l = 1:length(other_surfaces)
                n3 = other_surfaces(l).n;
                d3 = other_surfaces(l).d;
                
                A = [n1.'; n2.'; n3.'];
                b = [d1; d2; d3];
                if abs(det(A)) < 1e-9 %check that planes aren't parraell and thus illconditioned -> no intesection 
                     continue;
                end
                x = A\b; %least squares solve -> finds error region around point its trying to locate
                
                 if isInsidePlanes(x, planes, tol)
                    vertices(idx, :) = x.'; %row x,y,z
                    idx = idx +1;
                 end
            end
        end
    end
end
vertices = vertices(1:idx-1,:);
end