function fullPath = resolveRunPath(filename, numCams, runDir)
% resolveRunPath  Locate a per-run .mat (or .txt) file under the Results tree.
%
%   fullPath = resolveRunPath(filename) searches in this order:
%       1. As given (already absolute or relative-to-pwd that exists).
%       2. <projectRoot>/Results/<N>Cams/<filename> for N in {6,7,8},
%          where the cam count is parsed from the filename when possible.
%       3. <projectRoot>/Results/Logs/<filename> (for log files).
%       4. <projectRoot>/<filename> (legacy: file still at root).
%
%   fullPath = resolveRunPath(filename, numCams) restricts the cam-count
%   search to the supplied numCams when known (faster, avoids ambiguity).
%
%   fullPath = resolveRunPath(filename, numCams, runDir) checks runDir
%   and runDir/<N>Cams/ before falling back to the project root.
%
%   If no candidate exists, returns the most likely intended path so the
%   caller can produce a meaningful error.
%
%   This helper exists so logged RunFilename entries (stored as bare
%   basenames like '7Cams_Run_20260331_001929.mat') keep working after the
%   results were moved into Results/<N>Cams/.

    if nargin < 2, numCams = []; end
    if nargin < 3 || isempty(runDir), runDir = ''; end

    % 1. Use the path as-is if it already resolves.
    if isfile(filename)
        fullPath = filename;
        return;
    end

    % Project root is the parent folder of Analysis/ (where this file lives).
    projectRoot = fileparts(fileparts(mfilename('fullpath')));

    % If numCams was not supplied, try to parse it from the filename
    % (matches the "%dCams_Run_..." convention used by saveResults.m).
    if isempty(numCams)
        tok = regexp(filename, '^(\d+)Cams_', 'tokens', 'once');
        if ~isempty(tok)
            numCams = str2double(tok{1});
        end
    end

    candidates = {};

    if ~isempty(runDir)
        candidates{end+1} = fullfile(runDir, filename);
        if ~isempty(numCams)
            candidates{end+1} = fullfile(runDir, sprintf('%dCams', numCams), filename);
        end
    end

    if ~isempty(numCams)
        candidates{end+1} = fullfile(projectRoot, 'Results', sprintf('%dCams', numCams), filename);
    else
        for n = [6 7 8]
            candidates{end+1} = fullfile(projectRoot, 'Results', sprintf('%dCams', n), filename); %#ok<AGROW>
        end
    end

    % Logs / curated files
    candidates{end+1} = fullfile(projectRoot, 'Results', 'Logs', filename);
    candidates{end+1} = fullfile(projectRoot, 'Results', filename);
    candidates{end+1} = fullfile(projectRoot, filename);

    for i = 1:numel(candidates)
        if isfile(candidates{i})
            fullPath = candidates{i};
            return;
        end
    end

    % No match — return the most likely intended path so error messages
    % from the caller are informative.
    if ~isempty(numCams)
        fullPath = fullfile(projectRoot, 'Results', sprintf('%dCams', numCams), filename);
    else
        fullPath = fullfile(projectRoot, 'Results', filename);
    end
end
