function out = RunGA(problem, params, specs)
    % Specifications
    numCams = specs.Cams;
    SectionCentres = specs.SectionCentres;
    
    % Problem 
    CostFunction = problem.CostFunction;
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;
    nVar = problem.nVar;
    
    % Params 
    MaxIt = params.MaxIt;
    nPop = params.nPop;
    beta = params.beta; % selection pressure
    pC = params.pC; % percentage children (crossover) 
    nC = round(pC*nPop/2)*2; % number of offspring (needs to be even) -> number of elements in popc
    mu = params.mu; % percentage mutated
    sigma = params.sigma;

    % Template for Empty Individuals
    empty_individual.Chromosome = []; %chromosome 
    empty_individual.Cost = [];

    % Best Solution Ever Found 
    bestsol.Cost = inf;

    % Initialization 
    pop = repmat(empty_individual, nPop, 1);
    costs = zeros(nPop,1);
    pop(1).Chromosome = [ 2.442268657559323	4.500000000000000	2.053445556884173	2.031104950848467	0.001993426075398	0.298328077466267	4.406159222656118	-4.500000000000000	0.628774334251690	-1.504336780033701	-0.630674976872758	0.674503542729945	-2.055481818393845	-4.500000000000000	1.146039253354012	-1.572002378191020	-0.030396267365364	2.526621294823446	5.000000000000000	-1.101699569116084	4.687619266380412	2.932764681180698	-0.920198698149629	-1.709547586795046	-4.812818580459252	-3.009084393410572	4.800000000000000	-3.005420406819754	-0.017997508967026	-1.207042480548926	5.000000000000000	-4.323585646279842	1.173364976902663	-1.601278406431096	-0.913779712919404	1.344598154530524	-4.999493535876553	-4.077142830967917	0.825726417434871	-1.476973064671660	0.923198048282691	-0.867005773230155];
    pop(1).Cost = CostFunction(pop(1).Chromosome, specs);
    costs(1) = pop(1).Cost;    
    
    parfor i = 2:nPop
        % Generate Guided Random Solution 
        pop(i).Chromosome = initialPopulation(VarMin, VarMax, SectionCentres, numCams);
        % Evaluate Solution
        pop(i).Cost = CostFunction(pop(i).Chromosome, specs);
        costs(i) = pop(i).Cost;
    end
    % Find initial best Solution
    [minCost, minIdx] = min(costs);
    if minCost < bestsol.Cost
        bestsol = pop(minIdx);
    end 

    % Best Cost of Iterations -> record of best cost after each generation 
    bestcost = nan(MaxIt,1);
    bestChromosomes = nan(MaxIt, nVar);
    avgcost = nan(MaxIt,1); %average cost of the entire population 
    topTenAvgCost = nan(MaxIt,1); %average cost of the top 10 solutions 

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
            [popc(k,1).Chromosome, popc(k,2).Chromosome] = DoublePointCrossover(p1.Chromosome, p2.Chromosome,nVar); %[offspring1, offspring2]

        end 

        % Convert popc to Single Column (vertical) Matrix 
        popc = popc(:);

        % Mutation 
        for l = 1:nC 

            % Perform Mutation
            % adaptive_m
            % u = mu* exp(-it/MaxIt); % Decrease mutation rate as convergence improves
            popc(l).Chromosome = Mutate(popc(1).Chromosome, mu, sigma);

            % Check for variable bounds -> all variables must be greater than varMin
            % and less than varMax
            popc(l).Chromosome = max(popc(l).Chromosome, VarMin); % if greater than VarMin -> unchanged
            popc(l).Chromosome = min(popc(l).Chromosome, VarMax);

            % Evaluation 
            popc(l).Cost= CostFunction(popc(l).Chromosome, specs);
        end
        
        % Find best in offspring
        [minOffspringCost, minOffspringIdx] = min([popc.Cost]);
        if minOffspringCost< bestsol.Cost
            bestsol = popc(minOffspringIdx);
        end

        % Merge and Sort Populations
        pop = SortPopulation([pop;popc]);

        % Remove Extra Individuals 
        pop = pop(1:nPop);

        % Update Best Cost of Iteration
        bestcost(it) = bestsol.Cost;
        bestChromosomes(it,:) = bestsol.Chromosome;

        % Calculate Average Costs
        allCosts = [pop.Cost];
        avgcost(it) = mean(allCosts);
        topTenAvgCost(it) = mean(allCosts(1:10));
        
        % Display Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost = ', num2str(bestcost(it)), '  Avg Cost = ', num2str(avgcost(it))]);

    end

    % Results
    out.pop = pop;
    out.bestsol = bestsol;
    out.bestcost = bestcost;
    out.bestChromosomes = bestChromosomes;
    out.avgcost = avgcost;
    out.topTenAvgCost = topTenAvgCost;
end 