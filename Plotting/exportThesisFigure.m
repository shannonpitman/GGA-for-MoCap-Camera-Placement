function paths = exportThesisFigure(fig, outName, varargin)
%EXPORTTHESISFIGURE  Write one figure out in every format the thesis needs.
%
%   exportThesisFigure(fig, outName) writes
%       <outName>.pdf   vector, for LaTeX \includegraphics
%       <outName>.png   raster at 300 dpi, for Word / slides / quick review
%
%   outName carries no extension. Both files land beside each other so a
%   figure can be swapped between vector and raster without renaming.
%
%   Name-Value parameters:
%     'Formats'    cellstr subset of {'pdf','png'}.  Default both.
%     'Resolution' PNG dpi.                          Default 300.
%     'Background' exportgraphics BackgroundColor.   Default 'white'.
%     'Quiet'      true to suppress the "Saved:" line. Default false.
%
%   Why a helper: every plotGA_* / plot*_GAvsOptiTrack function used to
%   inline its own exportgraphics(fig, [outName '.pdf'], ...) call. Adding
%   a second output format meant editing thirteen files, and the export
%   options had already begun to drift between them. Centralising means
%   the whole figure set stays uniform — same background, same dpi, same
%   naming — which is exactly what a dissertation figure set needs.
%
%   See also gaPlotStyle, applyThesisStyle.

    p = inputParser;
    addParameter(p, 'Formats',    {'pdf', 'png'}, @(x) iscellstr(x) || ischar(x));
    addParameter(p, 'Resolution', 300,            @isnumeric);
    addParameter(p, 'Background', 'white');
    addParameter(p, 'Quiet',      false,          @islogical);
    parse(p, varargin{:});
    o = p.Results;

    formats = o.Formats;
    if ischar(formats), formats = {formats}; end

    % Strip an extension if the caller passed one by habit.
    [d, b, e] = fileparts(outName);
    if any(strcmpi(e, {'.pdf', '.png'}))
        outName = fullfile(d, b);
    end

    if ~isempty(d) && ~isfolder(d)
        mkdir(d);
    end

    paths = cell(1, numel(formats));

    for k = 1:numel(formats)
        switch lower(formats{k})
            case 'pdf'
                paths{k} = [outName '.pdf'];
                exportgraphics(fig, paths{k}, ...
                    'ContentType',     'vector', ...
                    'BackgroundColor', o.Background);
            case 'png'
                paths{k} = [outName '.png'];
                exportgraphics(fig, paths{k}, ...
                    'Resolution',      o.Resolution, ...
                    'BackgroundColor', o.Background);
            otherwise
                error('exportThesisFigure:badFormat', ...
                      'Unsupported format "%s".', formats{k});
        end
    end

    if ~o.Quiet
        [~, base] = fileparts(outName);
        fprintf('Saved: %s.{%s}\n', base, strjoin(formats, ','));
    end
end
