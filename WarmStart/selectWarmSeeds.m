function [seeds, info] = selectWarmSeeds(pool, nSeed, method, minDist, problem)
%SELECTWARMSEEDS  Chromosomes to inject into a warm-started GA population.
%
%   [seeds, info] = selectWarmSeeds(pool, nSeed, method, minDist, problem)
%
%   pool     struct array from loadSeedArchive (Chromosome, Cost), any order
%   nSeed    number of rows to inject (always returned exactly, if the pool
%            is not empty)
%   method   'top_n'       the nSeed best archived solutions, kept greedily
%                          by cost subject to a diversity constraint: a
%                          candidate is taken only if its layoutDistance to
%                          every seed already taken is >= minDist. If that
%                          leaves fewer than nSeed, the remaining slots are
%                          filled with the best rejected candidates, then
%                          (if the pool itself is smaller than nSeed) with
%                          mutated copies of the best.
%            'best_mutant' the old scheme, generalised: the best archived
%                          solution plus nSeed-1 copies of it, every gene
%                          mutated with sigma 0.5 (Mutate(x, 1, 0.5)).
%                          No distance check.
%   minDist  [m] diversity threshold for 'top_n' (0 turns it off)
%   problem  from buildRunSpecs (VarMin/VarMax for clamping mutants)
%
%   info.nFromPool / nDistinct / nFilled / nMutant / Costs (archive costs of
%   the injected pool rows, NaN for mutants) record exactly what was used.
%   RunGA wraps the rotation vectors and projects positions onto the mount
%   regions for every warm row, so nothing here needs to.

    if nargin < 4 || isempty(minDist), minDist = 0; end
    info = struct('Method', method, 'nRequested', nSeed, 'nFromPool', 0, ...
                  'nDistinct', 0, 'nFilled', 0, 'nMutant', 0, 'Costs', []);
    seeds = zeros(0, problem.nVar);
    if isempty(pool) || nSeed < 1
        return;
    end

    [~, order] = sort([pool.Cost]);
    pool  = pool(order);
    chrom = vertcat(pool.Chromosome);
    cost  = [pool.Cost].';
    if size(chrom, 2) ~= problem.nVar
        error('selectWarmSeeds:WrongLength', ...
            'Archive chromosomes have %d genes, this run needs %d.', size(chrom, 2), problem.nVar);
    end

    switch lower(method)
        case 'best_mutant'
            best = chrom(1, :);
            seeds = best;
            info.Costs = cost(1);
            for k = 2:nSeed
                seeds(k, :) = clampToBounds(Mutate(best, 1, 0.5), problem); %#ok<AGROW>
            end
            info.nFromPool = 1;
            info.nDistinct = 1;
            info.nMutant   = nSeed - 1;

        case 'top_n'
            taken = false(size(chrom, 1), 1);
            picked = [];
            for i = 1:size(chrom, 1)
                if numel(picked) >= nSeed, break; end
                far = true;
                if minDist > 0
                    for j = picked
                        if layoutDistance(chrom(i, :), chrom(j, :)) < minDist
                            far = false;
                            break;
                        end
                    end
                end
                if far
                    picked(end+1) = i; %#ok<AGROW>
                    taken(i) = true;
                end
            end
            info.nDistinct = numel(picked);

            % Top up with the best rejected candidates, by cost
            rest = find(~taken).';
            nFill = min(nSeed - numel(picked), numel(rest));
            picked = [picked, rest(1:nFill)];
            info.nFilled = nFill;

            seeds = chrom(picked, :);
            info.Costs = cost(picked);
            info.nFromPool = numel(picked);

            % Pool smaller than nSeed: mutated copies of the best
            for k = size(seeds, 1)+1:nSeed
                seeds(k, :) = clampToBounds(Mutate(chrom(1, :), 1, 0.5), problem);
                info.Costs(end+1) = NaN;
                info.nMutant = info.nMutant + 1;
            end

        otherwise
            error('selectWarmSeeds:UnknownMethod', ...
                'Unknown warm-seed method "%s" (use ''top_n'' or ''best_mutant'').', method);
    end
    info.Costs = info.Costs(:);
end

function x = clampToBounds(x, problem)
    x = min(max(x, problem.VarMin), problem.VarMax);
end
