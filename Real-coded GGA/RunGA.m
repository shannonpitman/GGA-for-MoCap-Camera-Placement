function out = RunGA(problem, params, specs)
    % Specifications
    numCams = specs.Cams;
    SectionCentres = specs.SectionCentres;
    warmStartUsed = specs.warmStart;
    warmStartChromosome = specs.warmChromosomes;

    % Problem 
    CostFunction = problem.CostFunction;
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;
    nVar = problem.nVar;
    
    % Params 
    MaxIt = params.MaxIt;
    nPop = params.nPop;
    % beta = params.beta; % selection pressure for roulette wheel 
    pC = params.pC; % percentage children (crossover) 
    nC = round(pC*nPop/2)*2; % number of offspring (needs to be even) -> number of elements in popc
    mu = params.mu; % percentage mutated
    sigma = params.sigma;
    TournamentSize = params.Tournamentsize;
    elitismDelay = params.elitismDelay;

    % Template for Empty Individuals
    empty_individual.Chromosome = []; %chromosome 
    empty_individual.Cost = [];

    % Best Solution Ever Found 
    bestsol.Cost = inf;

    % Initialization 
    pop = repmat(empty_individual, nPop, 1);
    costs = zeros(nPop,1);
    
    chromosomes = zeros(nPop, nVar); 

    for i = 1:nPop
        % Generate Guided Random Solution 
        chromosomes(i,:) = initialPopulation(VarMin, VarMax, SectionCentres, numCams);
    end
    if warmStartUsed
        for i = 1:size(warmStartChromosome,1)
            chromosomes(i,:)= warmStartChromosome(i,:);
        end
    end
    %Evaluate in Parallel
    parfor i = 1:nPop
        % Evaluate Solution
        costs(i) = CostFunction(chromosomes(i,:), specs);
    end
    
    %Assign Results 
    for i = 1:nPop
        pop(i).Chromosome = chromosomes(i,:);
        pop(i).Cost = costs(i);
    end

    % Find initial best Solution
    [minCost, minIdx] = min(costs);
    if minCost < bestsol.Cost
        bestsol = pop(minIdx);
    end 

    % Best Cost of Iterations -> record of best cost after each generation 
    bestcost = nan(MaxIt,1);
    bestChromosomes = nan(MaxIt, nVar);
    avgcost = nan(MaxIt,1); % average cost of the entire population 
    topTenAvgCost = nan(MaxIt,1); % average cost of the top 10 solutions 

    % Main Loop
    for it = 1:MaxIt
        % Uncomment if Roulette Wheel Selection is preferable
        % % Selection Probabilities
        % c = [pop.Cost];
        % avgc = mean(c);
        % if avgc ~=0
        %     c = c/avgc;
        % end
        % probs = exp(-beta*c); % Boltzmann distribution 
        % 
        % % Initialise Offsprings Population
        % popc = repmat(empty_individual, nC/2, 2); %Population children -> offspring 1 (column 1 ) offspring 2 (column 2)
        % 
        % % Parent chromosomes
        % parent_indices = RouletteWheelSelection(probs, nC*2);
        % parent_indices = reshape(parent_indices, nC, 2);
        % parent_chromosomes = zeros(nC, nVar);
        % 
        % for k = 1:nC
        %     parent_chromosomes(k,:) = pop(parent_indices(k,1)).Chromosome;
        % end
        % 

        %Tournament Selection 
        parent_chromosomes = Tournament(pop, nVar, nC, TournamentSize);

        popc = repmat(empty_individual, nC/2, 2);


        % Perform crossover with pre-extracted data
        for k = 1:nC/2
            idx1 = 2*k-1;
            idx2 = 2*k;
            [popc(k,1).Chromosome, popc(k,2).Chromosome] = DoublePointCrossover(parent_chromosomes(idx1,:), parent_chromosomes(idx2,:), numCams);
        end

        % Convert popc to Single Column (vertical) Matrix 
        popc = popc(:);

        % Mutation 
        for l = 1:nC
            % adaptive_m
            % u = mu* exp(-it/MaxIt); % Decrease mutation rate as convergence improves
            popc(l).Chromosome = Mutate(popc(l).Chromosome, mu, sigma);
    
            if mod(it, 10) == 0 || rand < 0.2  % or 20% chance each generation
                    popc(l).Chromosome = fixPoorCameras(popc(l).Chromosome, specs, 0.05);
            end
   
            popc(l).Chromosome = max(popc(l).Chromosome, VarMin); % if greater than VarMin -> unchanged
            popc(l).Chromosome = min(popc(l).Chromosome, VarMax);
            % Evaluation 
        end

        parfor l = 1:nC
            popc(l).Cost= CostFunction(popc(l).Chromosome, specs);
        end

        % Find best in offspring
        [minOffspringCost, minOffspringIdx] = min([popc.Cost]);
        if minOffspringCost< bestsol.Cost
            bestsol = popc(minOffspringIdx);
        end

        % Merge populations efficiently
        if it <= elitismDelay
            pop = SortPopulation(popc);
            % Remove Extra Individuals 
            pop = pop(1:nPop);
        else
            pop = SortPopulation([pop;popc]);
            % Remove Extra Individuals 
            pop = pop(1:nPop);
        end

        % Update Best Cost of Iteration
        bestcost(it) = bestsol.Cost;
        bestChromosomes(it,:) = bestsol.Chromosome;

        % Calculate Average Costs - vectorized
        allCosts = [pop.Cost];
        avgcost(it) = mean(allCosts);
        topTenAvgCost(it) = mean(allCosts(1:10));
        
        % Display Iteration Information
        if mod(it, 5) == 0 || it == 1
            fprintf('Iteration %3d: Best = %.6f, Avg = %.6f, Top10 = %.6f\n', it, bestcost(it), avgcost(it), topTenAvgCost(it));
        end

    end

    % Results
    out.pop = pop;
    out.bestsol = bestsol;
    out.bestcost = bestcost;
    out.bestChromosomes = bestChromosomes;
    out.avgcost = avgcost;
    out.topTenAvgCost = topTenAvgCost;
end 