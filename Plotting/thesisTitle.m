function h = thesisTitle(target, str, sty, varargin)
%THESISTITLE  A figure-spanning title that fits, and stays black.
%
%   h = thesisTitle(tl, str, sty)        % tiledlayout title
%   h = thesisTitle(fig, str, sty)       % sgtitle on a subplot figure
%   h = thesisTitle(..., 'MaxChars', 60) % wrap width, default 58
%
%   Two problems this fixes, both of which showed up across the whole
%   figure set:
%
%     1. Long one-line headings ran off both edges of the exported page.
%        MATLAB does not wrap or shrink them; it just clips. This wraps
%        the string onto as many lines as it needs, breaking at spaces
%        and preferring a break after ':' or before '(' or 'vs' so the
%        split lands somewhere meaningful.
%
%     2. Layout titles and sgtitles ignore axes-level styling and export
%        as a faded grey. This forces pure black, matching
%        applyThesisStyle's contract for everything else.
%
%   Name-Value parameters:
%     'MaxChars'   soft wrap width in characters.  Default 58
%     'FontSize'   default sty.FontSizeTitle
%     'FontWeight' default 'bold'
%
%   See also applyThesisStyle, gaPlotStyle.

    p = inputParser;
    addParameter(p, 'MaxChars',   58, @isnumeric);
    addParameter(p, 'FontSize',   [], @isnumeric);
    addParameter(p, 'FontWeight', 'bold', @ischar);
    parse(p, varargin{:});
    o = p.Results;

    if isempty(o.FontSize), o.FontSize = sty.FontSizeTitle; end

    lines = wrapTitle(char(str), o.MaxChars);

    args = {'FontSize',   o.FontSize, ...
            'FontWeight', o.FontWeight, ...
            'FontName',   sty.FontName, ...
            'Color',      'k'};

    if isa(target, 'matlab.graphics.layout.TiledChartLayout')
        h = title(target, lines, args{:});
    elseif isa(target, 'matlab.ui.Figure')
        h = sgtitle(target, lines, args{:});
        set(h, 'Color', 'k');
    else
        h = title(target, lines, args{:});
    end
end


function lines = wrapTitle(s, maxChars)
%WRAPTITLE  Greedy word wrap with preferred break points.
    s = strtrim(s);
    if numel(s) <= maxChars
        lines = {s};
        return;
    end

    words = strsplit(s, ' ');
    lines = {};
    cur   = '';
    for i = 1:numel(words)
        w = words{i};
        if isempty(cur)
            cand = w;
        else
            cand = [cur ' ' w];
        end

        if numel(cand) <= maxChars
            cur = cand;
            % A colon is a natural heading break — take it if we are
            % already past half the line, rather than dragging the next
            % clause up.
            if endsWith(cur, ':') && numel(cur) > maxChars/2
                lines{end+1} = cur; %#ok<AGROW>
                cur = '';
            end
        else
            if isempty(cur)
                lines{end+1} = w; %#ok<AGROW>
            else
                lines{end+1} = cur; %#ok<AGROW>
                cur = w;
            end
        end
    end
    if ~isempty(cur)
        lines{end+1} = cur;
    end
end
