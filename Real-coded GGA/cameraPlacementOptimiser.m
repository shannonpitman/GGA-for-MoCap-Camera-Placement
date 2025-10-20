onedrive = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Cape Town');
addpath(fullfile(onedrive, 'MSc', 'MATLAB', 'RVC3-MATLAB'));
addpath(fullfile(onedrive, 'MSc', 'MATLAB', 'Simulation', 'Optimising Camera Placement'));
%Camera Intrinsics
res = [2624 1964]; %720p 60fps
pixsize = [1.4e-6 1.4e-6]; %Square
pp = [res(1)/2 res(2)/2];
f = 0.0028; %focal length [m]
% Error Due to Camera Geometry 
%% 
% # Error Due to Change in Baseline in a parallel camera system
% # Error  Due to change in Tilt of cameras in a converging system 

P= [0 0 1]; %world point
num_cams = 2;

% System 1: Baseline Change, parallel cameras
baseline =linspace(0.1,2);
baseline_cm = baseline*100;
orientations = [[0, 0, 0]; [0, 0, 0]];
errorVolumeBase = zeros(size(baseline));

for i = 1:length(baseline)
    xL = -baseline(i)/2;
    xR = baseline(i)/2;
    positions = [[xL, 0, 0]; [xR, 0, 0]];
    errorVolumeBase(i) = resUncertainty(positions, orientations, num_cams, P, res, pixsize, f, pp);
end 
errorVolumeBase;
figure;
plot(baseline_cm, errorVolumeBase)
xlabel('Base Length (cm)')
ylabel('Error Volume (mm^3)')
grid on
title('Error Volume versus baseline length (parallel cameras)')
% System 2: Convergence Angle Change, constant baseline
positions = [[-1, 0, 0]; [1, 0, 0]];
tiltAngle = linspace(0, 90);
errorVolumeTilt = zeros(size(tiltAngle));
for j = 1:length(tiltAngle)
    orientationL = deg2rad(tiltAngle(j));
    orientationR = -deg2rad(tiltAngle(j));
    orientations = [[0, orientationL, 0];[0, orientationR, 0]];
    errorVolumeTilt(j) = resUncertainty(positions, orientations, num_cams, P, res, pixsize, f, pp);
end
errorVolumeTilt;
figure;
plot(tiltAngle, errorVolumeTilt)
xlabel('Tilt Angle (deg)')
ylabel('Error Volume (mm^3)')
grid on
title('Error Volume vs Tilt Angle (Converging Cameras)')
% Error Due to Point Position
%% 
% # Point moving away from a parallel camera system 
% # A matrix of points lying in a plane directly in front of the cameras in 
% (2.1) a parallel configuration and (2.2) a converging system. 

num_cams = 2;
positions = [[-1, 0, 0]; [1, 0, 0]];
%System 1: Point moving further away from parallel camera system
orientations = [[0, 0, 0]; [0, 0, 0]];
Points = linspace(0.5, 5.5);
P_cm = Points*100;
errorVolumePointDistance= zeros(size(Points));
for i = 1:numel(Points)
    Pz = Points(i);
    P = [0, 0, Pz];
    errorVolumePointDistance(i) = resUncertainty(positions, orientations, num_cams, P, res, pixsize, f,pp);
end 
errorVolumePointDistance;
figure;
plot(P_cm, errorVolumePointDistance)
xlabel('Point Distance in z-direction (cm)')
ylabel('Error Volume (mm^3)')
grid on
title('Error Volume vs 3D point moving away (Parallel Cameras)')
%%
num_cams = 2;
positions = [[-0.25, 0, 0]; [0.25, 0, 0]];

%{
 System 2: Matrix of Points lying in a plane directly in front of the cameras
%}
xPoints = -1:0.05:1;
yPoints = -1:0.05:1;
xPoints_cm= xPoints*100;
yPoints_cm = yPoints*100;
[X,Y] = meshgrid(xPoints, yPoints);

% 2.1 Parallel Configuration 
orientations = [[0, 0, 0]; [0, 0, 0]];

errorVolumePlanePar = zeros(length(xPoints), length(yPoints));
for i= 1:length(xPoints)
    for j = 1:length(yPoints)
        pX = xPoints(i);
        pY = yPoints(j);
        P = [pX, pY, 1];
        errorVolumePlanePar(j,i)= resUncertainty(positions, orientations, num_cams, P, res, pixsize, f,pp); %rows =y, cols =x
    end
end
errorVolumePlanePar;
cameras = cell(2,1);

for i = 1: 2
    T = se3(eul2rotm([0,0,0], "XYZ"), positions(i,:)); %camera to world cTw
    cameras{i} = CentralCamera(name="Camera "+i, pose=T);
end
cam1 = cameras{1};
cam2 = cameras{2};
% Boiler plate code
paraCamfig = figure;  % save the figure handle in a variable
surf(X,Y, errorVolumePlanePar)
hold on;
cam1.plot_camera("label", scale= 0.08)
cam2.plot_camera("label", scale= 0.08)

% Axis labels
xlabel('Point Distance in X(m)')
ylabel('Point Distance in Y(m)')
zlabel('Error Volume $(mm^3)$')
fname = 'conCamfig';
title('Error Volume vs 3D plane of Points (Parallel Cameras)')
grid on

% Boiler plate code
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(paraCamfig,'-property','FontSize'),'FontSize',25) % adjust fontsize to your document

set(findall(paraCamfig,'-property','Box'),'Box','off') % optional
set(findall(paraCamfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(paraCamfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(paraCamfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(paraCamfig,'Position');
set(paraCamfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(paraCamfig,fname,'-dpdf','-painters','-fillpage')

print(paraCamfig,fname,'-dpng','-painters')



%%
num_cams = 2;
positions = [[-1, 0, 0]; [1, 0, 0]];
% 2.2 Converging Configuration 
orientations = [[0, deg2rad(45), 0]; [0, deg2rad(-45), 0]];

errorVolumePlaneConv = zeros(length(xPoints), length(yPoints));
for i= 1:length(xPoints)
    for j = 1:length(yPoints)
        pX = xPoints(i);
        pY = yPoints(j);
        P = [pX, pY, 1];
        errorVolumePlaneConv(j,i)= resUncertainty(positions, orientations, num_cams, P, res, pixsize, f,pp);
    end
end
errorVolumePlaneConv;

cameras = cell(2,1);

for i = 1: 2
    T = se3(eul2rotm(orientations(i,:), "XYZ"), positions(i,:)); %camera to world cTw
    cameras{i} = CentralCamera(name="Camera "+i, pose=T);
end
cam1 = cameras{1};
cam2 = cameras{2};


% Boiler plate code
convCamfig = figure;  % save the figure handle in a variable
surf(X,Y, errorVolumePlaneConv)
hold on;
cam1.plot_camera("label")
cam2.plot_camera("label")

% Axis labels
xlabel('Point Distance in X(m)')
ylabel('Point Distance in Y(m)')
zlabel('Error Volume $(mm^3)$')
fname = 'conCamfig';
title('Error Volume vs 3D plane of Points (Converging Cameras)')
grid on

% Boiler plate code
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(convCamfig,'-property','FontSize'),'FontSize',25) % adjust fontsize to your document

set(findall(convCamfig,'-property','Box'),'Box','off') % optional
set(findall(convCamfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(convCamfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(convCamfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(convCamfig,'Position');
set(convCamfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(convCamfig,fname,'-dpdf','-painters','-fillpage')
print(convCamfig,fname,'-dpng','-painters')

%% 
%%
num_cams = 2;
positions = [[4, 0, 0]; [0, -4, 0]];
% 3. Perpendicular Configuration 
orientations = [[0, deg2rad(-90), 0]; [deg2rad(-90), deg2rad(0), 0]];
xPoints = -2:0.05:1;
yPoints = -1:0.05:2;
xPoints_cm= xPoints*100;
yPoints_cm = yPoints*100;
[X,Y] = meshgrid(xPoints_cm, yPoints_cm);
errorVolumePlanePerp = zeros(length(xPoints), length(yPoints));
for i= 1:length(xPoints)
    for j = 1:length(yPoints)
        pX = xPoints(i);
        pY = yPoints(j);
        P = [pX, pY, 0];
        errorVolumePlanePerp(j,i)= resUncertainty(positions, orientations, num_cams, P, res, pixsize, f,pp);
    end
end
errorVolumePlanePerp;

% Boiler plate code
convCamfig = figure;  % save the figure handle in a variable
surf(X,Y, errorVolumePlanePerp)

% Axis labels
xlabel('Point Distance in X(cm)')
ylabel('Point Distance in Y(cm)')
zlabel('Error Volume $(mm^3)$')
fname = 'conCamfig';
title('Error Volume vs 3D plane of Points (Perpendicular Cameras)')
grid on

% Boiler plate code
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(convCamfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(convCamfig,'-property','Box'),'Box','off') % optional
set(findall(convCamfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(convCamfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(convCamfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(convCamfig,'Position');
set(convCamfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(convCamfig,fname,'-dpdf','-painters','-fillpage')
print(convCamfig,fname,'-dpng','-painters')

%% 
% Multiple Camera Systems 
% # Place a 3rd Camera at a varying angle to a stereo pair
% # Place the 4 cameras at the vertices of a Square 

% System 1: 3 Cameras, where third faces and lies on a circumference of a cirle
% around the point and the other two cameras are a stereo pair 
num_cams = 3;
% theta = linspace(0,360);
theta = deg2rad(90);
P = [0,0,1]
Px = cos(theta)
Py = sin(theta)
positions = [[1, 0, 0]; [-1, 0, 0]; [Px, Py,0]];
zc= [-cos(theta), -sin(theta), 0] % tangent 
zc = zc/norm(zc)
xc =cross(zc, P)
xc = xc/norm(xc)
yc = cross(zc,xc)
yc = yc/norm(yc)
R_cam3 = [xc(:), yc(:), zc(:)]
first2camsOrient = eul2rotm([0, 0, 0])

orientations = [eye(3); eye(3); R_cam3];
cameras = cell(num_cams, 1);
for i = 1: num_cams
    T = se3(orientations(3*(i-1)+1:3*i,:), positions(i,:)); %camera to world cTw
    cameras{i} = CentralCamera(name="cam"+i,resolution= res, pixel= 1.4e-6, focal= f, pose=T, center = pp);
end
cam1 = cameras{1};
cam2 = cameras{2};
cam3 = cameras{3}

figure;
plotsphere(P, 0.15, "b"); 
hold on;
axis([-2 3 -3 3 0 4]);
cam1.plot_camera("label")
cam2.plot_camera("label")
cam3.plot_camera("label")
hold off;
title('Marker Grid and Camera Pose');
xlabel('X (m)'); 
ylabel('Y (m)'); 
zlabel('Z (m)');
grid on; 
view(45,30);