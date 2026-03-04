function sectionCentres = generateSectionCentres(numCams, volume)
% subdivide workspace into sections to orient cameras towards closest
% section centre

numSections = floor(numCams/2);
z_mid = mean(volume(3,:)); %mid-height
spaceCentre = [mean(volume(1,:)), mean(volume(2,:)), z_mid]; %centre of the flightspace
if mod(numSections, 2)== 1
    %Odd = add the centre as a viewpoint, distribute the rest on grid
    gridSections = numSections -1;
    centrePoint = spaceCentre;
else
    gridSections = numSections;
    centrePoint = [];
end

if gridSections > 0
    % Find most balanced grid
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
    x_centres = linspace(volume(1,1), volume(1,2), 2*nx+1);
    x_centres = x_centres(2:2:end); % midpoints
    y_centres = linspace(volume(2,1), volume(2,2), 2*ny+1);
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
sectionCentres = [centrePoint; gridCentres];