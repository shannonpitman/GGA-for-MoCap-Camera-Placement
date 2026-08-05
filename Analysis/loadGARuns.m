function [filteredLog, filterDesc] = loadGARuns(varargin)
% loadGARuns  Load and filter GA runs from a master log file.
%   [filteredLog, filterDesc] = loadGARuns('Name', Value, ...)
%
%   Returns the filtered runLog struct array and a human-readable string
%   describing the active filters (useful for filenames and annotations).
%
%   Name-Value Parameters:
%     'LogFile'      - Path to .mat log file
%                      (default: <projectRoot>/Results/Logs/GGA_RunsLog.mat)
%     'NumCameras'   - Scalar camera count to keep
%     'CostFunction' - Scalar cost function type (1, 2, or 3)
%     'WarmStart'    - Logical true/false
%     'TargetType'   - 1 (UAV) or 2 (UGV)
%     'GridMode'     - 1 (Uniform) or 2 (Normal)
%     'Spacing'      - Grid spacing in metres (tolerance 0.001)
%     'MatchParams'  - true to keep only runs whose GAParams match the
%                      first passing run (default: false)

    % Default log lives under <projectRoot>/Results/Logs/.
    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'NumCameras',   [],  @isnumeric);
    addParameter(p, 'CostFunction', [],  @isnumeric);
    addParameter(p, 'WarmStart',    [],  @islogical);
    addParameter(p, 'TargetType',   [],  @isnumeric);
    addParameter(p, 'GridMode',     [],  @isnumeric);
    addParameter(p, 'Spacing',      [],  @isnumeric);
    addParameter(p, 'MatchParams',  false, @islogical);
    parse(p, varargin{:});

    opts = p.Results;

    %% Load — fall back to the legacy root-level location if needed.
    if ~isfile(opts.LogFile)
        legacy = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'GGA_RunsLog.mat');
        if isfile(legacy)
            opts.LogFile = legacy;
        else
            error('loadGARuns:fileNotFound', ...
                  'Log file not found: %s', opts.LogFile);
        end
    end
    load(opts.LogFile, 'runLog');
    N = length(runLog);
    fprintf('Loaded %d total GA runs from %s\n', N, opts.LogFile);

    %% Backwards compatibility — ensure newer fields exist
    for fn = {'TargetType', 'GridMode', 'Spacing', 'NumTargetPoints'}
        if ~isfield(runLog, fn{1})
            [runLog.(fn{1})] = deal(NaN);
        end
    end

    %% Build filter mask
    keep = true(N, 1);
    descParts = {};

    if ~isempty(opts.NumCameras)
        keep = keep & ([runLog.NumCameras]' == opts.NumCameras);
        descParts{end+1} = sprintf('%dC', opts.NumCameras);
    end
    if ~isempty(opts.CostFunction)
        keep = keep & ([runLog.CostFunctionType]' == opts.CostFunction);
        descParts{end+1} = sprintf('CF%d', opts.CostFunction);
    end
    if ~isempty(opts.WarmStart)
        keep = keep & ([runLog.WarmStart]' == opts.WarmStart);
        if opts.WarmStart, descParts{end+1} = 'Warm';
        else,              descParts{end+1} = 'Cold'; end
    end
    if ~isempty(opts.TargetType)
        keep = keep & ([runLog.TargetType]' == opts.TargetType);
        descParts{end+1} = sprintf('TT%d', opts.TargetType);
    end
    if ~isempty(opts.GridMode)
        keep = keep & ([runLog.GridMode]' == opts.GridMode);
        descParts{end+1} = sprintf('GM%d', opts.GridMode);
    end
    if ~isempty(opts.Spacing)
        keep = keep & (abs([runLog.Spacing]' - opts.Spacing) < 0.001);
        descParts{end+1} = sprintf('sp%.0fcm', opts.Spacing*100);
    end

    %% Match GA parameters (optional)
    if opts.MatchParams && sum(keep) > 1
        refIdx = find(keep, 1);
        ref = runLog(refIdx).GAParams;
        for i = 1:N
            if keep(i)
                g = runLog(i).GAParams;
                if g.MaxIt ~= ref.MaxIt || g.nPop ~= ref.nPop || ...
                   g.beta  ~= ref.beta  || g.mu   ~= ref.mu
                    keep(i) = false;
                end
            end
        end
        descParts{end+1} = 'matchedParams';
    end

    %% Apply
    filteredLog = runLog(keep);
    nFilt = length(filteredLog);

    if nFilt == 0
        error('loadGARuns:noMatch', 'No runs match the specified filters.');
    end

    if isempty(descParts)
        filterDesc = 'all';
    else
        filterDesc = strjoin(descParts, '_');
    end

    fprintf('  -> %d runs after filtering (%s)\n\n', nFilt, filterDesc);
end