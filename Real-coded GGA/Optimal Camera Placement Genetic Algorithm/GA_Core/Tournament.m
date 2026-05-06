function selectedChromosomes = Tournament(pop, nVar, n, k)
% TOURNAMENTSELECTION Vectorized tournament selection
%   costs: vector of fitness values (lower is better)
%   n: number of variables
%   k: tournament size 
%
    selectedChromosomes = zeros(n, nVar);
    N = length(pop);
    costs = [pop.Cost];
    
    % Each row is one tournament, columns are the k participants
    tournaments = randi(N, n, k); 
    
    % Get costs for all participants: n x k matrix
    tournamentCosts = costs(tournaments);
    
    % Find winner (minimum cost) in each tournament
    [~, winnerCol] = min(tournamentCosts, [], 2);
    
    % Convert column indices to actual population indices
    % Use linear indexing: (row-1)*k + col gives linear index into tournaments
    linearIdx = (1:n)' + (winnerCol - 1) * n;
    selectedIndices = tournaments(linearIdx);
    
    for i = 1:n
        selectedChromosomes(i, :) = pop(selectedIndices(i)).Chromosome;
    end
end