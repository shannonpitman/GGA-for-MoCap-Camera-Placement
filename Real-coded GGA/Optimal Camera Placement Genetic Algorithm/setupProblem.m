function problem = setupProblem(numCams, costFunctionType, camUpperBounds, camLowerBounds)
% Select cost function
switch costFunctionType
    case 1
        problem.CostFunction = @resUncertaintyCost;
    case 2
        problem.CostFunction = @dynamicOcclusionCost;
    case 3
        problem.CostFunction = @combinedCostFunction;
end

problem.VarMin = repmat(camLowerBounds,1,numCams);
problem.VarMax = repmat(camUpperBounds,1,numCams);
problem.nVar = 6*numCams;
end