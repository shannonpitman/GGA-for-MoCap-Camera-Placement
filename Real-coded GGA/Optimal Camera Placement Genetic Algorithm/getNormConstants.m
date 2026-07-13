function nc = getNormConstants(varargin)
%GETNORMCONSTANTS  Look up per-instance normalisation constants.
%
% Retrieves the normalisation constants built by buildNormalisationSchedule
% for a given instance (TargetType x GridMode x NumCameras x Spacing) from
% the stored normTable.mat.
%
% USAGE:
%   nc = getNormConstants('TargetType',1,'GridMode',1,'NumCameras',7);
%   nc = getNormConstants(specs);  % pull the key straight off a specs struct
%   [uNorm,oNorm] = deal(nc.uncertNorm, nc.occlNorm);
%
% Name-Value parameters (or pass a specs struct as the first argument):
%   'TargetType' - 1 (UAV) or 2 (UGV).           Default 1
%   'GridMode'   - 1 (Uniform) or 2 (Normal).    Default 1
%   'NumCameras' - camera count.                 Default 7
%   'Spacing'    - grid spacing [m].             Default 1.0
%   'File'       - normTable.mat path. Default Results/normTable.mat
%
% OUTPUT struct nc with fields:
%   uncertNorm, occlNorm, utopiaUnc, nadirUnc, utopiaOcc, nadirOcc,
%   and the matched TargetType/GridMode/NumCameras/Spacing.
%
% Errors if no matching instance is stored.

    % Allow getNormConstants(specs) where specs carries the instance key.
    leading = {};
    if nargin >= 1 && isstruct(varargin{1})
        s = varargin{1};
        if isfield(s,'TargetType'), leading = [leading, {'TargetType', s.TargetType}]; end
        if isfield(s,'TargetMode'), leading = [leading, {'GridMode',   s.TargetMode}]; end
        if isfield(s,'Cams'),       leading = [leading, {'NumCameras', s.Cams}];       end
        if isfield(s,'spacing'),    leading = [leading, {'Spacing',    s.spacing}];    end
        varargin = [leading, varargin(2:end)];
    end

    projectRoot = fileparts(mfilename('fullpath'));

    p = inputParser;
    addParameter(p, 'TargetType', 1, @isnumeric);
    addParameter(p, 'GridMode', 1, @isnumeric);
    addParameter(p, 'NumCameras', 7, @isnumeric);
    addParameter(p, 'Spacing', 1.0, @isnumeric);
    addParameter(p, 'File', fullfile(projectRoot, 'Results', 'normTable.mat'), @ischar);
    parse(p, varargin{:});
    o = p.Results;

    if ~isfile(o.File)
        error('getNormConstants:noFile', 'Normalisation table not found: %s\nRun buildNormalisationSchedule first.', o.File);
    end

    S = load(o.File, 'normTable');
    if ~isfield(S, 'normTable') || isempty(S.normTable)
        error('getNormConstants:empty', 'normTable in %s is empty.', o.File);
    end
    T = S.normTable;

    idx = find([T.TargetType] == o.TargetType & [T.GridMode] == o.GridMode & ...
        [T.NumCameras] == o.NumCameras & abs([T.Spacing] - o.Spacing) < 1e-9, 1);

    if isempty(idx)
        error('getNormConstants:noMatch', ...
            'No stored instance for TargetType=%d GridMode=%d NumCameras=%d Spacing=%.3g.', ...
            o.TargetType, o.GridMode, o.NumCameras, o.Spacing);
    end

    r = T(idx);
    nc = struct( ...
        'uncertNorm', r.uncertNorm, ...
        'occlNorm',   r.occlNorm, ...
        'utopiaUnc',  r.utopiaUnc, ...
        'nadirUnc',   r.nadirUnc, ...
        'utopiaOcc',  r.utopiaOcc, ...
        'nadirOcc',   r.nadirOcc, ...
        'TargetType', r.TargetType, ...
        'GridMode',   r.GridMode, ...
        'NumCameras', r.NumCameras, ...
        'Spacing',    r.Spacing);
end
