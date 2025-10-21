function Chromosome = initialPopulation(VarMin, VarMax,  SectionCentres, numCams)
% This function generates a guided inital population where the randomised
% camera locations are oriented towards the (proportional to amount of
% cameras) subdivided workspace 
Chromosome = zeros(1, 6*numCams);
for c = 1:numCams
    camPositions = unifrnd(VarMin(1:3),VarMax(1:3), [1,3]);
    distanceCam2Centre = vecnorm(SectionCentres -camPositions,2,2); %calcs 2-norm of each row 
    [~, closestIdx] = min(distanceCam2Centre); %find closest camera section
    nearest_centre = SectionCentres(closestIdx,:);
    %Orient optical (z-) axis towards Interest Point
    distance = nearest_centre-camPositions;
    unitDist = distance/norm(distance);
    Z_axis = [0,0,1];
    rotAng = acos(dot(Z_axis, unitDist));
    rotAxis = cross(Z_axis,unitDist);
    q = quaternion([cos(rotAng/2),rotAxis*sin(rotAng/2)]);
    q = normalize(q);
    q2Eul = euler(q, 'XYZ', 'point'); %euler angles in radians
    alpha = q2Eul(1);
    beta = q2Eul(2);
    gamma = q2Eul(3);
    gene = [camPositions, [alpha,beta, gamma]];

    chromStartIdx = (c-1)*6+1;
    chromEndIdx = c*6;
    Chromosome(chromStartIdx:chromEndIdx) = gene;    

end
end