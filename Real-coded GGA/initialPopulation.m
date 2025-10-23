function Chromosome = initialPopulation(VarMin, VarMax,  SectionCentres, numCams)
% This function generates a guided inital population where the randomised
% camera locations are placed on random boundary faces oriented towards the (proportional to amount of
% cameras) subdivided workspace 
Chromosome = zeros(1, 6*numCams);
faceID = zeros(1, numCams);

for c = 1:numCams
    chosenFace = randi(5);
    faceID(c) = chosenFace; %store which face camera is placed on
    
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
    distances = vecnorm(SectionCentres - camPos, 2, 2); %calcs 2-norm of each row 
    [~, closestIdx] = min(distances); %find closest camera section
    nearest_centre = SectionCentres(closestIdx,:);

    % Orient optical (z-) axis toward the interest point
    directionVec = nearest_centre-camPos;
    directionUnit = directionVec/norm(directionVec);

    Z_axis = [0,0,1];
    rotAng = acos(dot(Z_axis, directionUnit));
    rotAxis = cross(Z_axis, directionUnit);
    rotAxis = rotAxis/norm(rotAxis);

    q = quaternion([cos(rotAng/2),rotAxis*sin(rotAng/2)]);
    q = normalize(q);
    q2Eul = euler(q, 'XYZ', 'point'); %euler angles in radians
    alpha = q2Eul(1);
    beta = q2Eul(2);
    gamma = q2Eul(3);
    gene = [camPos, [alpha,beta, gamma]];

    Chromosome((c-1)*6+1:c*6) = gene;    

end
end