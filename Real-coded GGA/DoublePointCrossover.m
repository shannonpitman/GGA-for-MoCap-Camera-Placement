function [y1, y2] = DoublePointCrossover(x1,x2, numCams)
    q = randperm(numCams);
    j1 = min(q(1), q(2)); %Crossover point is the camera index
    j2 = max(q(1), q(2));
    % Convert camera indices to gene indices
    geneEnd1 = j1*6;            % Last gene of camera j1
    geneEnd2 = j2*6;            % Last gene of camera j2
    
    % Swap camera blocks between crossover points
    y1 = [x1(1:geneEnd1), x2(geneEnd1+1:geneEnd2), x1(geneEnd2+1:end)];
    y2 = [x2(1:geneEnd1), x1(geneEnd1+1:geneEnd2), x2(geneEnd2+1:end)];
end
