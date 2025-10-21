function errorVolume = resUncertainty(cameraChromosome, point)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

resolution = [2624 1964]; %720p 60fps
pixelSize = [1.4e-6 1.4e-6]; %Square
pp = [res(1)/2 res(2)/2];
focalLength = 0.0028; %focal length [m]

numCams = cameraChromosome/6;
cameras = cell(numCams,1);
uv = cell(numCams, 2);
cameraCentre = cell(numCams, 3);
worldPoints = cell(3*numCams,4);
planes = cell(numCams,1);
adjacentSurfaces = [1 2; 2 3; 3 4; 4 1];

du = 0.5; %pixels
dv = 0.5; %pixels

for i = 1:numCams
    chromStartIdx = (i-1)*6+1;
    chromEndIdx = i*6;
    camPositions = cameraChromosome(chromStartIdx: ChromStartIdx+2);
    camOrientations = cameraChromosome(chromEndIdx-2: chromEndIdx);

    T = se3(eul2rotm(camOrientations, "XYZ"), camPositions); %camera to world cTw
    cameras{i} = CentralCamera(name="cam"+i,resolution= resolution, pixel= pixelSize, focal= focalLength, pose=T, center = pp);
    uv{i} = cameras{i}.project(point); %uv projected coordinates 
    u = uv{i}(1);
    v = uv{i}(2);
    if ~(u >= 1 && u <= resolution(1) && v >= 1 && v <= resolution(2))
        errorVolume = NaN; %point is not visible in camera FOV
        return; %exit res calc function 
    end
    cameraCentre{i} = cameras{i}.center().'; %world location of Camera center
    worldPoints{i} = quantToWorld(cameras{i}, uv{i}, du, dv, cameraCentre{i});
    planes{i} = buildPyramidSurf(cameraCentre{i}, worldPoints{i}, adjacentSurfaces);
end

vertices = calcVertices(numCams, adjacentSurfaces, planes);
unique_vertices = unique(round(vertices,6),'rows', 'stable');
% size(unique_vertices);
V_ = unique_vertices-point; %centre the points at the mean 
C0v_vert = cov(V_); %covariance matrix
[A,~] = eig(C0v_vert); %[eigenvectors, eigenvalues in a diagonal matrix]
transformed_v = A*V_.';%rotated points 

%Initial guess
axes0 = max(abs(transformed_v.'), [], 1).' ; %1x3 [a0,b0,c0], max magnitude of the x,y,z components of v 

w2 = 0.2; %weight for second energy function
fun = @(axes) optimiseEllipsoid(axes, transformed_v, axes0, w2);
x0 = axes0;
options = optimset('Display', 'notify', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000, 'MaxFunEvals', 5000);
[axes_opt, ~] = fminsearch(fun, x0, options);

%optimised axes [mm]
a = axes_opt(1)*1000;
b = axes_opt(2)*1000;
c = axes_opt(3)*1000;

errorVolume= 4/3*pi*a*b*c;%mm^3


end