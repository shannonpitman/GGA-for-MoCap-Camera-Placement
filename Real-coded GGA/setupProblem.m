function problem = setupProblem(numCams, costFunctionType, camUpperBounds, camLowerBounds)
% Select cost function
switch costFunctionType
    case 1
        problem.CostFunction = @resUncertainty; %moved camera chromosome to combined cost function 
    case 2
        problem.CostFunction = @dynamicOcclusion;
    case 3
        problem.CostFunction = @combinedCostFunction;
end

problem.VarMin = repmat(camLowerBounds,1,numCams);
problem.VarMax = repmat(camUpperBounds,1,numCams);
problem.nVar = 6*numCams;