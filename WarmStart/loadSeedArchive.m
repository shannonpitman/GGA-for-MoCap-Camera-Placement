function pool = loadSeedArchive(preset, condKey)
%LOADSEEDARCHIVE  Every archived solution for one preset and condition.
%
%   pool = loadSeedArchive(preset, condKey)
%
%   pool is a struct array (empty if no archive yet) with fields
%     Chromosome  1 x 6N row (rotation-vector genes)
%     Cost        cost under the condition's cost function
%     Source      'best' (a run's best solution) or 'finalpop' (another
%                 member of that run's final population)
%     RunFile     run .mat basename it came from
%     WarmStart   whether that run was warm-started
%     Added       datetime it was archived
%   sorted by Cost, best first.
    f = seedArchiveFile(preset, condKey);
    if ~isfile(f)
        pool = emptyPool();
        return;
    end
    S = load(f, 'archive');
    pool = S.archive;
    if isempty(pool)
        pool = emptyPool();
        return;
    end
    [~, order] = sort([pool.Cost]);
    pool = pool(order);
end

function pool = emptyPool()
    pool = struct('Chromosome', {}, 'Cost', {}, 'Source', {}, ...
                  'RunFile', {}, 'WarmStart', {}, 'Added', {});
end
