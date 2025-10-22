function worldPoints = quantToWorld(camera,uv, du,dv, cameraCentre)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    uCorners = [uv(1)-du; uv(1)-du; uv(1)+du; uv(1)+du];
    vCorners = [uv(2)-dv; uv(2)+dv; uv(2)+dv; uv(2)-dv];
    
    homogPix = [uCorners; vCorners; ones(1,4)];
    camRays = camera.K \ homogPix;
    worldPoints = camera.T.rotm*camRays + cameraCentre;

end 