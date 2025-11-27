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
    
    chromosomes = zeros(nPop, nVar);

    % pop(1).Chromosome = [ 2.442268657559323	4.500000000000000	2.053445556884173	2.031104950848467	0.001993426075398	0.298328077466267	4.406159222656118	-4.500000000000000	0.628774334251690	-1.504336780033701	-0.630674976872758	0.674503542729945	-2.055481818393845	-4.500000000000000	1.146039253354012	-1.572002378191020	-0.030396267365364	2.526621294823446	5.000000000000000	-1.101699569116084	4.687619266380412	2.932764681180698	-0.920198698149629	-1.709547586795046	-4.812818580459252	-3.009084393410572	4.800000000000000	-3.005420406819754	-0.017997508967026	-1.207042480548926	5.000000000000000	-4.323585646279842	1.173364976902663	-1.601278406431096	-0.913779712919404	1.344598154530524	-4.999493535876553	-4.077142830967917	0.825726417434871	-1.476973064671660	0.923198048282691	-0.867005773230155];
    % pop(1).Cost = CostFunction(pop(1).Chromosome, specs);
    % costs(1) = pop(1).Cost;    
    % pop(2).Chromosome = [5.000000000000000	3.135293348324546	1.795549488573229	1.648898066818321	-0.809489616591172	-0.732073310686882	5.000000000000000	-4.320533111305447	0.507980219777034	-1.476013469425586	-0.774576311879347	-0.080573556444475	4.896745385074944	-2.240068458379641	1.424049127293288	-2.215676212732957	-1.282317601389465	1.967342403343143	5.000000000000000	-3.237069360343214	0.458565590845076	-1.170764661173703	-1.068140036684545	1.023287923603299	5.000000000000000	1.822852881123080	1.606594661711947	1.736046137426712	0.763478556750995	1.157088648554378	-4.981946829345784	-3.754979006399330	1.061232114819257	-1.608409025747987	0.877630914268858	-0.901046250819087	-4.974169068633028	-0.838988812216320	0.760472457820022	1.260067040873139	1.329552685313434	1.265891584248071];
    % pop(2).Cost = CostFunction(pop(2).Chromosome, specs);
    % costs(2) = pop(2).Cost;   

    for i = 1:nPop
        % Generate Guided Random Solution 
        chromosomes(i,:) = initialPopulation(VarMin, VarMax, SectionCentres, numCams);
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
        popc = repmat(empty_individual, nC/2, 2); %Population children -> offspring 1 (column 1 ) offspring 2 (column 2)
        
        % Pre-extract all parent chromosomes
        parent_indices = RouletteWheelSelection(probs, nC*2);
        parent_indices = reshape(parent_indices, nC, 2);
        parent_chromosomes = zeros(nC, nVar);

        for k = 1:nC
            parent_chromosomes(k,:) = pop(parent_indices(k,1)).Chromosome;
        end
        
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
        pop = SortPopulation([pop;popc]);

        % Remove Extra Individuals 
        pop = pop(1:nPop);

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