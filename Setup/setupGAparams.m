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
params.sigma = [repmat(cfg.MutationSigmaPos, 1, 3), repmat(cfg.MutationSigmaRot, 1, 3)]; % per camera: [m m m rad rad rad]
params.Tournamentsize = cfg.TournamentSize;
params.elitismDelay = 0; %MaxIt/3; % prioritise diversity initially 
