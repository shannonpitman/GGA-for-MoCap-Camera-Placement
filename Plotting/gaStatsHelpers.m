function helpers = gaStatsHelpers()
% gaStatsHelpers  Bundle of statistical / annotation helpers used by
% plotGA_* functions.
%
%   h = gaStatsHelpers();
%   p     = h.mannWhitney(x, y);            % two-sided U-test p-value
%   p     = h.wilcoxonSigned(x, y);         % paired signed-rank p-value
%   stars = h.sigStars(p);                  % '' | 'n.s.' | '*' | '**' | '***'
%   h.drawSigBracket(ax, x1, x2, y, label); % bracket + significance label
%
% These wrappers prefer Statistics Toolbox (ranksum / signrank) but fall
% back to a normal-approximation implementation so the plotting code does
% not hard-fail on a minimal MATLAB install.

    helpers.mannWhitney     = @mannWhitneyU;
    helpers.wilcoxonSigned  = @wilcoxonSignedRank;
    helpers.sigStars        = @significanceStars;
    helpers.drawSigBracket  = @drawSigBracket;
end


% ---------------------------------------------------------------------
function p = mannWhitneyU(x, y)
% Two-sided Mann-Whitney U-test p-value (a.k.a. Wilcoxon rank-sum).
% Uses Statistics Toolbox `ranksum` if available, else a normal
% approximation (with tie correction) suitable for n >= 4.

    x = x(:); y = y(:);
    x = x(~isnan(x));
    y = y(~isnan(y));
    if isempty(x) || isempty(y)
        p = NaN; return;
    end

    if exist('ranksum', 'file') == 2
        try
            p = ranksum(x, y);
            return;
        catch
            % fall through to manual implementation
        end
    end

    n1 = numel(x);
    n2 = numel(y);
    combined = [x; y];
    [~, ord] = sort(combined);
    ranks = zeros(size(combined));
    ranks(ord) = 1:numel(combined);
    % Tie correction: average ranks for equal values
    [uVals, ~, gIdx] = unique(combined);
    for k = 1:numel(uVals)
        m = gIdx == k;
        if sum(m) > 1
            ranks(m) = mean(ranks(m));
        end
    end

    R1 = sum(ranks(1:n1));
    U1 = R1 - n1*(n1+1)/2;
    U  = min(U1, n1*n2 - U1);

    % Normal approximation with tie correction
    nT = n1 + n2;
    [~, ~, gIdx] = unique(combined);
    tieCounts = accumarray(gIdx, 1);
    tieTerm   = sum(tieCounts.^3 - tieCounts);

    muU    = n1*n2/2;
    sigmaU = sqrt(n1*n2/12 * ((nT+1) - tieTerm/(nT*(nT-1))));
    if sigmaU == 0
        p = 1; return;
    end
    z = (U - muU) / sigmaU;
    p = 2 * (1 - 0.5*(1 + erf(abs(z)/sqrt(2))));
    p = max(min(p, 1), 0);
end


% ---------------------------------------------------------------------
function p = wilcoxonSignedRank(x, y)
% Two-sided Wilcoxon signed-rank p-value for paired samples.
% Uses Statistics Toolbox `signrank` if available, else normal
% approximation.

    x = x(:); y = y(:);
    n = min(numel(x), numel(y));
    if n < 2, p = NaN; return; end
    d = x(1:n) - y(1:n);
    d = d(d ~= 0);
    if isempty(d), p = 1; return; end

    if exist('signrank', 'file') == 2
        try
            p = signrank(x(1:n), y(1:n));
            return;
        catch
            % fall through
        end
    end

    absD = abs(d);
    [~, ord] = sort(absD);
    ranks = zeros(size(absD));
    ranks(ord) = 1:numel(absD);
    [uVals, ~, gIdx] = unique(absD);
    for k = 1:numel(uVals)
        m = gIdx == k;
        if sum(m) > 1
            ranks(m) = mean(ranks(m));
        end
    end
    Wpos = sum(ranks(d > 0));
    Wneg = sum(ranks(d < 0));
    W    = min(Wpos, Wneg);

    nW = numel(d);
    muW    = nW*(nW+1)/4;
    sigmaW = sqrt(nW*(nW+1)*(2*nW+1)/24);
    if sigmaW == 0
        p = 1; return;
    end
    z = (W - muW) / sigmaW;
    p = 2 * (1 - 0.5*(1 + erf(abs(z)/sqrt(2))));
    p = max(min(p, 1), 0);
end


% ---------------------------------------------------------------------
function s = significanceStars(p)
% Map p-value to APA-style significance label.
    if isnan(p)
        s = '';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end


% ---------------------------------------------------------------------
function drawSigBracket(ax, x1, x2, y, label, fontSize)
% Draw a horizontal bracket between (x1,y) and (x2,y) with a centred
% label slightly above. Used to display significance labels above
% paired box plots.
    if nargin < 6 || isempty(fontSize), fontSize = 8; end
    tickH = 0.02 * range(ylim(ax));
    if tickH <= 0, tickH = 0.01; end

    line(ax, [x1, x1, x2, x2], [y - tickH, y, y, y - tickH], ...
        'Color', 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    text(ax, (x1+x2)/2, y + 0.4*tickH, label, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', fontSize, 'FontName', 'Helvetica', ...
        'Color', 'k');
end
