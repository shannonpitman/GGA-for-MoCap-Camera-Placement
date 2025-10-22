function planes = buildPyramidSurf(cameraCentre, worldPoint, adj)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
corner1 = worldPoint(:,adj(:,1));
corner2 = worldPoint(:,adj(:,2));

vector1 = corner1 - cameraCentre; 
vector2 = corner2 - cameraCentre;
planes.n = cross(vector2, vector1, 1);  %normal to the plane of surface 
planes.n = surfNormal./vecnorm(surfNormal,2,1); %normalise to only extract direction
planes.d = sum(planes.n .* cameraCentre,1); %constant of plane calculated by dot product of vector2 (corner 2) with a point on the surface (perspective centre)



%planes = struct('n', [], 'd', []);
% for a = 1:4
%     c1 = worldPoint(:,adj(a,1));
%     c2 = worldPoint(:,adj(a,2));
%     [n, d] = calcPyramidSurf(cameraCentre, c1, c2);
%     planes(a).n = n; %normal
%     planes(a).d = d; %constant 
% end
% end