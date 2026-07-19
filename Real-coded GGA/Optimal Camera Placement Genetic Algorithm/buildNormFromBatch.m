function normTable = buildNormFromBatch(schedule, cfg, varargin)
%BUILDNORMFROMBATCH  Derive CF3 normalisation constants from batch CF1/CF2 runs.
%
% Instead of a separate calibration pass (buildNormalisationSchedule), this
% reads the CF1 (resolution) and CF2 (occlusion) runs already completed in a
% batchRunGA schedule and builds the per-instance normalisation table from
% them. Because the batch's CF1/CF2 runs use the FULL production budget
% (repeats + warm start), the utopia they define is a genuine lower bound, so
% CF3 = (raw - utopia)/norm can no longer go negative on runs from the same
% batch.
%
% For each instance (TargetType x GridMode x NumCameras x Spacing) that has at
% least one completed CF1 AND one completed CF2 run:
%   utopiaUnc = min BestCost over CF1 runs;  chromUnc = its chromosome
%   utopiaOcc = min BestCost over CF2 runs;  chromOcc = its chromosome
%   nadirOcc  = dynamicOcclusion(chromUnc)   (occlusion at resolution optimum)
%   nadirUnc  = resUncertainty(chromOcc)     (resolution at occlusion optimum)
%   uncertNorm = max(nadirUnc - utopiaUnc, eps)
%   occlNorm   = max(nadirOcc - utopiaOcc, eps)
%
% Writes Results/normTable.mat in the same schema as buildNormalisationSchedule
% so getNormConstants / setupCostParams read it transparently. Any prior table
% is archived alongside first (normTable_backup_<timestamp>.mat).
%
% USAGE:
%   buildNormFromBatch(schedule, cfg);                 % called inside batchRunGA
%   S = load('Results/Logs/BatchLog_XXXX.mat');
%   buildNormFromBatch(S.schedule, S.batchConfig);     % standalone from a log
%
% Name-Value:
%   'OutFile' - normTable.mat path. Default Results/normTable.mat
%   'Verbose' - print per-instance summary. Default true
%   'Archive' - back up any existing table before overwrite. Default true

    addProjectPaths();

    p = inputParser;
    addParameter(p, 'OutFile', '', @ischar);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'Archive', true, @islogical);
    parse(p, varargin{:});
    o = p.Results;

    projectRoot = fileparts(mfilename('fullpath'));
    if isempty(o.OutFile)
        o.OutFile = fullfile(projectRoot, 'Results', 'normTable.mat');
    end
    if ~isfolder(fileparts(o.OutFile))
        mkdir(fileparts(o.OutFile));
    end

    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform', 'Normal'};

    % Keep only completed runs (failed / pending carry no usable cost).
    done = schedule(strcmp({schedule.Status}, 'done'));
    if isempty(done)
        error('buildNormFromBatch:noRuns', 'No completed runs in schedule.');
    end

    cf = [done.CostFunc];
    tt = [done.TargetType];
    gm = [done.GridMode];
    nc = [done.NumCams];
    sp = [done.Spacing];

    % Suppress any incidental figures during cross-evaluation.
    prevVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanupObj = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis)); %#ok<NASGU>

    % Unique instances, independent of cost function.
    keyMat = [tt(:), gm(:), nc(:), sp(:)];
    [instances, ~, ic] = unique(keyMat, 'rows', 'stable');
    ic = ic(:)';
    nInst = size(instances, 1);

    normTable = struct([]);
    nBuilt = 0;

    for k = 1:nInst
        ittt = instances(k, 1);  igm = instances(k, 2);
        inc  = instances(k, 3);  isp = instances(k, 4);

        inThis = (ic == k);
        idxCF1 = find(inThis & cf == 1);
        idxCF2 = find(inThis & cf == 2);

        if isempty(idxCF1) || isempty(idxCF2)
            if isempty(idxCF1), miss = 'CF1'; else, miss = 'CF2'; end
            if o.Verbose
                fprintf('  [skip] %s/%s %dC sp%.2f: no completed %s run\n', ...
                    ttNames{ittt}, gmNames{igm}, inc, isp, miss);
            end
            continue;
        end

        cf1Costs = [done(idxCF1).BestCost];
        cf2Costs = [done(idxCF2).BestCost];
        % min ignores NaN by default; completed runs carry a real BestCost.
        [utopiaUnc, a1] = min(cf1Costs);
        [utopiaOcc, a2] = min(cf2Costs);
        chromUnc = done(idxCF1(a1)).BestChromosome;
        chromOcc = done(idxCF2(a2)).BestChromosome;

        % Reconstruct the instance specs exactly as batchRunGA does, then
        % cross-evaluate: the opposite objective at each single-objective
        % optimum defines the nadir.
        specs = buildInstanceSpecs(inc, ittt, igm, isp, cfg);
        nadirOcc = evalComponent(chromUnc, specs, 'occ');
        nadirUnc = evalComponent(chromOcc, specs, 'unc');

        uncertNorm = max(nadirUnc - utopiaUnc, eps);
        occlNorm   = max(nadirOcc - utopiaOcc, eps);

        row = struct( ...
            'TargetType',   ittt, ...
            'GridMode',     igm, ...
            'NumCameras',   inc, ...
            'Spacing',      isp, ...
            'utopiaUnc',    utopiaUnc, ...
            'nadirUnc',     nadirUnc, ...
            'uncertNorm',   uncertNorm, ...
            'utopiaOcc',    utopiaOcc, ...
            'nadirOcc',     nadirOcc, ...
            'occlNorm',     occlNorm, ...
            'chromUnc',     chromUnc, ...
            'chromOcc',     chromOcc, ...
            'replicateUnc', cf1Costs(:), ...
            'replicateOcc', cf2Costs(:), ...
            'elapsedSec',   sum([done(idxCF1).ElapsedTime, done(idxCF2).ElapsedTime]));

        normTable = appendRow(normTable, row);
        nBuilt = nBuilt + 1;

        if o.Verbose
            fprintf(['  %s/%s %dC sp%.2f | utopiaUnc=%.4g nadirUnc=%.4g -> uncertNorm=%.4g' ...
                     ' | utopiaOcc=%.4g nadirOcc=%.4g -> occlNorm=%.4g\n'], ...
                ttNames{ittt}, gmNames{igm}, inc, isp, ...
                utopiaUnc, nadirUnc, uncertNorm, utopiaOcc, nadirOcc, occlNorm);
        end
    end

    if isempty(normTable)
        error('buildNormFromBatch:noInstances', ...
            'No instance had both a completed CF1 and CF2 run; nothing written.');
    end

    % Archive any prior table before overwriting.
    if o.Archive && isfile(o.OutFile)
        [d, n, e] = fileparts(o.OutFile);
        bak = fullfile(d, sprintf('%s_backup_%s%s', n, datestr(now, 'yyyymmdd_HHMMSS'), e));
        copyfile(o.OutFile, bak);
        if o.Verbose, fprintf('  archived previous table -> %s\n', bak); end
    end

    saveTable(o.OutFile, normTable, cfg, ttNames, gmNames);
    if o.Verbose
        fprintf('buildNormFromBatch: wrote %d instance(s) to %s\n', nBuilt, o.OutFile);
    end
end

%% Helpers
function specs = buildInstanceSpecs(numCams, tt, gm, sp, cfg)
% Mirror the specs construction in batchRunGA's main loop so the cross-eval
% is on an identical target space / geometry.
    volume = getfielddef(cfg, 'Volume', [-4 4; -4 4; 0 4]);
    maxH   = getfielddef(cfg, 'UGV_MaxHeight', 0.5);
    zSpc   = getfielddef(cfg, 'UGV_ZSpacing', 0.25);

    if tt == 2
        volume(3, :)  = [0, maxH];
        zSp           = min(zSpc, maxH);
        targetSpacing = [sp, sp, zSp];
    else
        targetSpacing = sp;
    end

    specs = setupHardwareSpecs(numCams);
    specs.WeightUncertainty = 0.5;
    specs.WeightOcclusion   = 0.5;
    specs.TargetType = tt;
    specs.TargetMode = gm;
    specs.Target = generateTargetSpace(volume, gm, targetSpacing);
    specs.NumPoints = size(specs.Target, 1);
    specs.spacing = sp;
    if tt == 2
        specs.spacingZ = zSp;
    end
    specs.SectionCentres = generateSectionCentres(numCams, volume);
    specs.UseNormTable = false;   % must not depend on the file we are writing
    specs = setupCostParams(specs);
end

function val = evalComponent(chrom, specs, which)
    [cameras, camCenters] = setupCameras(chrom, specs.Cams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    switch which
        case 'unc'
            val = resUncertainty(specs, cameras, camCenters);
        case 'occ'
            val = dynamicOcclusion(specs, cameras, camCenters);
        otherwise
            error('evalComponent:bad', 'which must be ''unc'' or ''occ''.');
    end
end

function tbl = appendRow(tbl, row)
    if isempty(tbl)
        tbl = row;
    else
        tbl(end + 1) = orderfields(row, tbl);
    end
end

function v = getfielddef(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end

function saveTable(outFile, normTable, cfg, ttNames, gmNames)
% Persist the struct array, a readable table view, and provenance metadata,
% matching buildNormalisationSchedule's on-disk schema.
    ps = getfielddef(cfg, 'PopulationSize', []);
    if isempty(ps)
        popMeta = 'auto: numCams*numParams*10';
    else
        popMeta = ps;
    end

    meta = struct( ...
        'CreatedOn',      datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        'Method',         'batch-derived (utopia = min CF1/CF2 BestCost, nadir = cross-eval)', ...
        'Source',         'buildNormFromBatch', ...
        'Repeats',        getfielddef(cfg, 'NumRepeats', NaN), ...
        'Generations',    getfielddef(cfg, 'MaxGenerations', NaN), ...
        'PopulationSize', popMeta, ...
        'WarmStart',      ~getfielddef(cfg, 'SkipWarmStart', false));

    tt = [normTable.TargetType]';
    gm = [normTable.GridMode]';
    normSummary = table( ...
        tt, gm, [normTable.NumCameras]', [normTable.Spacing]', ...
        [normTable.utopiaUnc]', [normTable.nadirUnc]', [normTable.uncertNorm]', ...
        [normTable.utopiaOcc]', [normTable.nadirOcc]', [normTable.occlNorm]', ...
        'VariableNames', {'TargetType','GridMode','NumCameras','Spacing', ...
        'utopiaUnc','nadirUnc','uncertNorm','utopiaOcc','nadirOcc','occlNorm'});
    normSummary.TargetTypeName = ttNames(tt)';
    normSummary.GridModeName   = gmNames(gm)';

    save(outFile, 'normTable', 'normSummary', 'meta', '-v7.3');
end
