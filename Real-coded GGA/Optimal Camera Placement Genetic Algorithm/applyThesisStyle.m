function applyThesisStyle(target, varargin)
% APPLYTHESISSTYLE  Force white background + black text on a figure or axes.
%   applyThesisStyle(fig)       — applies to a figure and all its axes
%   applyThesisStyle(ax)        — applies to a single axes object
%   applyThesisStyle(...,'IncludeColorbars', false)  — skip colorbar restyling
%
% This helper enforces the thesis style requirement that every figure has
% a pure white background and pure black axis lines, ticks, labels,
% titles, legends, and (where present) colorbar text. It is called near
% the end of every plotGA_* function so that no axes-level styling can
% override the global look.
%
% Why a helper instead of inline code: MATLAB's default
% AxesColor/GridColor varies subtly by release, and inline styling
% creates drift between plot files. Centralising guarantees consistency
% across every figure produced for the thesis.

    p = inputParser;
    addParameter(p, 'IncludeColorbars', true, @islogical);
    parse(p, varargin{:});
    incCB = p.Results.IncludeColorbars;

    % --- Resolve figure / axes handle list ----------------------------
    if isempty(target) || ~isgraphics(target)
        return;
    end

    if isa(target, 'matlab.ui.Figure')
        fig = target;
        % Recurse: also style every axes in the figure.
        set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
        axList = findall(fig, 'Type', 'axes');
        for k = 1:numel(axList)
            styleSingleAxes(axList(k));
        end
        % sgtitle — colour the supertitle text black if present
        st = findall(fig, 'Type', 'text', '-and', 'Tag', 'suptitle');
        for k = 1:numel(st)
            set(st(k), 'Color', 'k');
        end
        if incCB
            cbList = findall(fig, 'Type', 'colorbar');
            for k = 1:numel(cbList)
                styleColorbar(cbList(k));
            end
        end
    else
        styleSingleAxes(target);
    end
end

% =====================================================================

function styleSingleAxes(ax)
% Apply pure-black axes/text/grid on a single axes handle.
    if ~isgraphics(ax), return; end

    % Axis line and tick colour
    set(ax, ...
        'Color',           'w', ...
        'XColor',          'k', ...
        'YColor',          'k', ...
        'ZColor',          'k', ...
        'GridColor',       'k', ...
        'MinorGridColor',  'k', ...
        'GridAlpha',       0.15, ...
        'MinorGridAlpha',  0.10);

    % Title and labels — explicit black so themes / dark backgrounds
    % cannot bleed through.
    if ~isempty(ax.Title) && isgraphics(ax.Title)
        set(ax.Title,  'Color', 'k');
    end
    if ~isempty(ax.XLabel) && isgraphics(ax.XLabel)
        set(ax.XLabel, 'Color', 'k');
    end
    if ~isempty(ax.YLabel) && isgraphics(ax.YLabel)
        set(ax.YLabel, 'Color', 'k');
    end
    if ~isempty(ax.ZLabel) && isgraphics(ax.ZLabel)
        set(ax.ZLabel, 'Color', 'k');
    end

    % Legend (if attached)
    if ~isempty(ax.Legend) && isgraphics(ax.Legend)
        set(ax.Legend, ...
            'Color',     'w', ...
            'TextColor', 'k', ...
            'EdgeColor', 'k');
    end

    % Free-floating text annotations parented to this axes
    txt = findall(ax, 'Type', 'text');
    for k = 1:numel(txt)
        % Only recolour text whose colour is currently the MATLAB
        % default (dark grey). Preserve user-set colours such as
        % red 'overall best' annotations.
        c = get(txt(k), 'Color');
        if isnumeric(c) && all(abs(c - [0.15 0.15 0.15]) < 0.05)
            set(txt(k), 'Color', 'k');
        end
        % Force any background boxes to white if they are off-white
        bg = get(txt(k), 'BackgroundColor');
        if (ischar(bg) && strcmpi(bg, 'none')) || isempty(bg)
            % leave as-is
        end
    end
end

% =====================================================================

function styleColorbar(cb)
% Apply black axis/text on a colorbar handle.
    if ~isgraphics(cb), return; end
    set(cb, 'Color', 'k');
    if ~isempty(cb.Label) && isgraphics(cb.Label)
        set(cb.Label, 'Color', 'k');
    end
end
