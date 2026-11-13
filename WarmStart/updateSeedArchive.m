function nAdded = updateSeedArchive(preset, condKey, chroms, costs, sources, runFile, warmStart)
%UPDATESEEDARCHIVE  Add solutions to the warm-start archive of one condition.
%
%   nAdded = updateSeedArchive(preset, condKey, chroms, costs, sources, runFile, warmStart)
%
%   chroms   K x 6N chromosomes      costs  K x 1
%   sources  K x 1 cellstr ('best' / 'finalpop'), or one char for all
%
%   Chromosomes already in the archive (identical to 1e-9) are skipped, so
%   re-running a backfill is harmless. Called by saveResults after every run
%   (best + top fraction of the final population) and by backfillSeedArchive.
    if isempty(chroms), nAdded = 0; return; end
    if ischar(sources), sources = repmat({sources}, size(chroms, 1), 1); end

    f = seedArchiveFile(preset, condKey);
    if ~isfolder(fileparts(f)), mkdir(fileparts(f)); end
    archive = loadSeedArchive(preset, condKey);
    if isempty(archive)
        existing = zeros(0, size(chroms, 2));
    else
        existing = vertcat(archive.Chromosome);
    end

    nAdded = 0;
    now_ = datetime('now');
    for k = 1:size(chroms, 1)
        c = chroms(k, :);
        if ~isempty(existing) && size(existing, 2) == numel(c) && ...
                any(all(abs(existing - c) < 1e-9, 2))
            continue;
        end
        archive(end+1) = struct('Chromosome', c, 'Cost', costs(k), ...
            'Source', sources{k}, 'RunFile', runFile, ...
            'WarmStart', logical(warmStart), 'Added', now_); %#ok<AGROW>
        existing(end+1, :) = c; %#ok<AGROW>
        nAdded = nAdded + 1;
    end

    if nAdded > 0
        [~, order] = sort([archive.Cost]);
        archive = archive(order); %#ok<NASGU>
        save(f, 'archive');
    end
end
