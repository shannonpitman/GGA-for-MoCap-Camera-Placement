function [seeds, info] = warmSeedsForRun(cfg, condKey, problem, params)
%WARMSEEDSFORRUN  Warm-start rows for one run, from the condition's archive.
%
%   [seeds, info] = warmSeedsForRun(cfg, condKey, problem, params)
%
%   Seed count: cfg.WarmSeedCount if set, otherwise
%   round(cfg.WarmSeedFraction * population size). Method and diversity
%   threshold come from cfg.WarmSeedMode / cfg.WarmSeedMinDist.
%   Empty seeds (and a warning) if the archive has nothing for condKey yet.
    nSeed = numSeeds(cfg, params.nPop);
    pool = loadSeedArchive(cfg.Preset, condKey);
    if isempty(pool)
        warning('warmSeedsForRun:EmptyArchive', ...
            'No archived solutions for %s / %s: this warm run starts cold.', cfg.Preset, condKey);
    end
    [seeds, info] = selectWarmSeeds(pool, nSeed, cfg.WarmSeedMode, cfg.WarmSeedMinDist, problem);
    info.PoolSize = numel(pool);
    fprintf('    Warm start: %d seeds (%s; %d distinct, %d filled, %d mutants) from %d archived\n', ...
        size(seeds, 1), info.Method, info.nDistinct, info.nFilled, info.nMutant, info.PoolSize);
end

function n = numSeeds(cfg, nPop)
    if ~isempty(cfg.WarmSeedCount)
        n = cfg.WarmSeedCount;
    else
        n = round(cfg.WarmSeedFraction * nPop);
    end
    n = max(1, min(n, nPop));
end
