function surfStruct = buildPyramidSurf(cameraCentre,worldPoint, adj)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
surfStruct = struct('n', {}, 'd', {});

for a = 1:4
    c1 = worldPoint(:,adj(a,1));
    c2 = worldPoint(:,adj(a,2));
    [n, d] = calcPyramidSurf(cameraCentre, c1, c2);
    surfStruct(a).n = n; %normal
    surfStruct(a).d = d; %constant 
end
end