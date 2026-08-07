function worldPoints = quantToWorld(camera,u,v, du,dv, cameraCentre)
%Finds the world points for all corners of the uncertain quantisation region 
    uCorners = [u-du, u-du, u+du, u+du];
    vCorners = [v-dv, v+dv, v+dv, v-dv];
    
    homogPix = [uCorners; vCorners; ones(1,4)];
    worldPoints = camera.T.rotm*(camera.K \ homogPix) + cameraCentre;

end 