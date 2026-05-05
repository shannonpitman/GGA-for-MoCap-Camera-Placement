clear, clc
flight_envelope = [-4 4; -4 4; 0 4];
numCams =13;

numSections = floor(numCams/2);

z_mid = mean(flight_envelope(3,:));
spaceCentre = [mean(flight_envelope(1,:)), mean(flight_envelope(2,:)), z_mid];

if mod(numSections, 2) == 1
    % Odd: place one at centre, distribute rest on grid
    gridSections = numSections - 1;
    centrePoint = spaceCentre;
else
    gridSections = numSections;
    centrePoint = [];
end

if gridSections > 0
    % Find most balanced grid for gridSections
    bestNx = gridSections; bestNy = 1;
    for tryNx = 1:ceil(sqrt(gridSections))
        tryNy = ceil(gridSections / tryNx);
        if tryNx * tryNy >= gridSections
            if abs(tryNx - tryNy) < abs(bestNx - bestNy)
                bestNx = tryNx;
                bestNy = tryNy;
            end
        end
    end
    nx = bestNx;
    ny = bestNy;
    
    % Evenly spaced centres across the full grid
    x_centres = linspace(flight_envelope(1,1), flight_envelope(1,2), 2*nx+1);
    x_centres = x_centres(2:2:end); % midpoints
    y_centres = linspace(flight_envelope(2,1), flight_envelope(2,2), 2*ny+1);
    y_centres = y_centres(2:2:end);
    
    [Xc, Yc] = meshgrid(x_centres, y_centres);
    gridCentres = [Xc(:), Yc(:), repmat(z_mid, nx*ny, 1)];
    
    % If more grid cells than needed, keep most central ones
    if size(gridCentres,1) > gridSections
        dists = vecnorm(gridCentres - spaceCentre, 2, 2);
        [~, sortIdx] = sort(dists);
        gridCentres = gridCentres(sortIdx(1:gridSections), :);
    end
else
    gridCentres = [];
end

% Combine centre point (if odd) with grid centres
specs.SectionCentres = [centrePoint; gridCentres];

% --- Verification scatter plot ---
figure('Name', 'Section Centres Verification');
scatter3(specs.SectionCentres(:,1), specs.SectionCentres(:,2), specs.SectionCentres(:,3), ...
    100, 'r', 'filled');
hold on;

% Draw flight envelope boundaries
xe = flight_envelope(1,:);
ye = flight_envelope(2,:);
ze = flight_envelope(3,:);
corners = [
    xe(1) ye(1) ze(1); xe(2) ye(1) ze(1); xe(2) ye(2) ze(1); xe(1) ye(2) ze(1);  % bottom
    xe(1) ye(1) ze(2); xe(2) ye(1) ze(2); xe(2) ye(2) ze(2); xe(1) ye(2) ze(2)]; % top
edges = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
for e = 1:size(edges,1)
    plot3(corners(edges(e,:),1), corners(edges(e,:),2), corners(edges(e,:),3), 'b--', 'LineWidth', 1);
end

% Label each centre
for i = 1:size(specs.SectionCentres,1)
    text(specs.SectionCentres(i,1)+0.2, specs.SectionCentres(i,2)+0.2, specs.SectionCentres(i,3), ...
        sprintf('S%d', i), 'FontSize', 10, 'FontWeight', 'bold');
end

axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Section Centres (%d sections for %d cameras)', numSections, numCams));
view(3);
hold off;