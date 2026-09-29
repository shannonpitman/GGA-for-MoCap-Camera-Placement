function Chromosome = initialPopulation(VarMin, VarMax,  SectionCentres, numCams, mountRegions)
% This function generates a guided inital population where the randomised
% camera locations are placed on random boundary faces oriented towards the (proportional to amount of
% cameras) subdivided workspace 
    Chromosome = zeros(1, 6*numCams);
    camPositions = zeros(numCams, 3);
    
    % Positions: uniform within a random mountable region (mountRegions).
    % Without regions, the legacy rule: a random face of the search box.
    useMount = nargin >= 5 && ~isempty(mountRegions);

    %Pregenerate random values 
    faceIDs = randi(5,1, numCams);
    
    for c = 1:numCams
        if useMount
            camPositions(c,:) = sampleMount(mountRegions);
            continue;
        end
        chosenFace = faceIDs(c);
        
        switch chosenFace
            case 1 % +X
                camPos = [VarMax(1), unifrnd(VarMin(2), VarMax(2)), unifrnd(VarMin(3), VarMax(3))];
            case 2 % -X
                camPos = [VarMin(1), unifrnd(VarMin(2), VarMax(2)), unifrnd(VarMin(3), VarMax(3))];
            case 3 % +Y
                camPos = [unifrnd(VarMin(1), VarMax(1)),VarMax(2), unifrnd(VarMin(3), VarMax(3))];
            case 4 % -Y
                camPos = [unifrnd(VarMin(1), VarMax(1)), VarMin(2), unifrnd(VarMin(3), VarMax(3))];
            case 5 % +Z
                camPos = [unifrnd(VarMin(1), VarMax(1)), unifrnd(VarMin(2), VarMax(2)), VarMax(3)];
        end
        camPositions(c,:) = camPos;
    end

    for c =1:numCams
        camPos = camPositions(c, :); % Retrieve the camera position for the current camera
    
        distances = vecnorm(SectionCentres - camPos, 2, 2); %calcs 2-norm of each row 
        [~, closestIdx] = min(distances); %find closest camera section
        nearest_centre = SectionCentres(closestIdx,:);
    
        % Orient optical (z-) axis toward the interest point
        directionVec = nearest_centre-camPos;
        directionUnit = directionVec/norm(directionVec);
    
        % Axis-angle aim stored directly as the rotation-vector genes
        gene = [camPos, aimRotvec(directionUnit)];
        Chromosome((c-1)*6+1:c*6) = gene;    
    end
end