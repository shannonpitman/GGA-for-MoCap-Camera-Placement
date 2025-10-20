function [surfNormal,surfConst] = calcPyramidSurf(cameraCentre,corner1, corner2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    vector1 = corner1-cameraCentre; 
    vector2 = corner2-cameraCentre;
    surfNormal= cross((vector2), (vector1));  %normal to the plane of surface 
    surfNormal = surfNormal/norm(surfNormal); %normalise to only extract direction
    surfConst= dot(surfNormal, cameraCentre); %constant of plane calculated by dot product of vector2 (corner 2) with a point on the surface (perspective centre)
end