function out = RunGA(problem, params, specs)
    % Specifications
    numCams = specs.Cams;
    % res = specs.Resolution;
    % pixSize = specs.PixelSize;
    % focalLength = specs.Focal;
    % PrincipalPoint = specs.PrinsipalPoint;
    SectionCentres = specs.SectionCentres;
    % TargetSpace = specs.Target;
    
    % Problem 
    CostFunction = problem.CostFunction;
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;
    
    % Params 
    MaxIt = params.MaxIt;
    nPop = params.nPop;
    beta = params.beta; % selection pressure
    pC = params.pC; % percentage children (crossover) 
    nC = round(pC*nPop/2)*2; % number of offspring (needs to be even) -> number of elements in popc
    gamma = params.gamma;
    mu = params.mu; % percentage mutated
    sigma = params.sigma;

    % Template for Empty Individuals
    empty_individual.Chromosome = []; %chromosome 
    empty_individual.Cost = [];

    % Best Solution Ever Found 
    bestsol.Cost = inf;

    % Initialization 
    pop = repmat(empty_individual, nPop, 1);
    for i = 1:nPop

        % Generate Guided Random Solution 
        pop(i).Chromosome = initialPopulation(VarMin, VarMax, SectionCentres, numCams);

        % Evaluate Solution
        pop(i).Cost = CostFunction(pop(i).Chromosome, specs);

        % Compare Solution to Best Solution Ever found 
        if pop(i).Cost < bestsol.Cost
            bestsol = pop(i);
        end 
    end

    % Best Cost of Iterations -> record of best cost after each generation 
    bestcost = nan(MaxIt,1);

    % Main Loop
    for it = 1:MaxIt

        % Selection Probabilities
        c = [pop.Cost];
        avgc = mean(c);
        if avgc ~=0
            c = c/avgc;
        end
        probs = exp(-beta*c);

        % Initialise Offsprings Population
        popc = repmat(empty_individual, nC/2, 2);  %Population children -> offspring 1 (column 1 ) offspring 2 (column 2)

        %Crossover 
        for k = 1:nC/2

            % Select Parents
            p1 = pop(RouletteWheelSelection(probs));
            p2 = pop(RouletteWheelSelection(probs));

            % Perform Crossover 
            [popc(k,1).Chromosome, popc(k,2).Chromosome] = UniformCrossover(p1.Chromosome, p2.Chromosome, gamma); %[offspring1, offspring2]

        end 

        % Convert popc to Single Column (vertical) Matrix 
        popc = popc(:);

        % Mutation 
        for l = 1:nC 

            % Perform Mutation
            popc(l).Chromosome = Mutate(popc(1).Chromosome, mu, sigma);

            % Check for variable bounds -> all variables must be greater than varMin
            % and less than varMax
            popc(l).Chromosome = max(popc(l).Chromosome, VarMin); % if greater than VarMin -> unchanged
            popc(l).Chromosome = min(popc(l).Chromosome, VarMax);

            % Evaluation 
            popc(l).Cost = CostFunction(popc(l).Chromosome, specs);

            % Compare Solution to Best Solution Ever Found 
            if popc(l).Cost < bestsol.Cost
                bestsol = pop(l);
            end
        end

        % Merge and Sort Populations
        pop = SortPopulation([pop;popc]);

        % Remove Extra Individuals 
        pop = pop(1:nPop);

        % Update Best Cost of Iteration
        bestcost(it) = bestsol.Cost;
        
        % Display Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost = ', num2str(bestcost(it))]);

    end

    % Results
    out.pop = pop;
    out.bestsol = bestsol;
    out.bestcost = bestcost;
end 