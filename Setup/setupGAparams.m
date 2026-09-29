function params = setupGAparams(MaxIt, nPop, cfg)
% GA parameters. Mutation, crossover and tournament settings come from a
% runConfig struct when one is given; otherwise the historical defaults.
if nargin < 3 || isempty(cfg)
    cfg = runConfig();
end
params.MaxIt = MaxIt;
params.nPop = nPop;
params.beta = 1;
params.pC = cfg.CrossoverFraction; %probability of crossover 
params.gamma = 0.1;
params.mu = cfg.MutationRate; %probability of mutation
params.sigma = cfg.MutationSigma;
params.Tournamentsize = cfg.TournamentSize;
params.elitismDelay = 0; %MaxIt/3; % prioritise diversity initially 
