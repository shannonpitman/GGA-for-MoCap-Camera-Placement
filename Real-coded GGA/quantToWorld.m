function worldPoints = quantToWorld(camera,uv, du,dv, cameraCentre)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
uCorners = [uv(1)-du; uv(1)-du; uv(1)+du; uv(1)+du];
vCorners = [uv(2)-dv; uv(2)+dv; uv(2)+dv; uv(2)-dv];
worldPoints = zeros(3,4);
for j = 1:4
    homogPix = [uCorners(j); vCorners(j); 1];
    camRay = camera.K\homogPix;
    worldPoints(:,j) = camera.T.rotm*camRay + cameraCentre;
end
end 