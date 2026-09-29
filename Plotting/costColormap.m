function cmap = costColormap(n)
%COSTCOLORMAP  Sequential green -> light blue -> blue -> purple -> red map.
%
%   Low cost reads as green (good), high cost as red (bad). Five distinct
%   hues give more visible steps than a two-hue ramp, so small cost
%   differences stand out. Deliberately avoids yellow and near-white:
%   both wash out when projected onto a bright screen.
%
%   cmap = costColormap()     % 256 levels
%   cmap = costColormap(64)

    if nargin < 1 || isempty(n), n = 256; end

    anchors = [0.10 0.60 0.28;   % green          (low)
               0.30 0.70 0.90;   % light blue
               0.12 0.30 0.78;   % blue
               0.52 0.20 0.64;   % purple
               0.84 0.11 0.13];  % red            (high)

    x  = linspace(0, 1, size(anchors,1));
    xq = linspace(0, 1, n);
    cmap = [interp1(x, anchors(:,1), xq, 'linear')', ...
            interp1(x, anchors(:,2), xq, 'linear')', ...
            interp1(x, anchors(:,3), xq, 'linear')'];
end
