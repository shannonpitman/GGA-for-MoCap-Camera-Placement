function s = gaPlotStyle()
% gaPlotStyle  Publication-quality defaults for GA results figures.
%   s = gaPlotStyle() returns a struct of styling parameters used by all
%   plotGA_* functions. Centralised here so every figure is consistent.
%
%   Designed for a single-column thesis layout (~160 mm text width).
%   Figures exported as vector PDF via exportgraphics.
%
%   THESIS STYLE LOCK: every plotGA_* function ends with a call to
%   applyThesisStyle(fig) which forces white background + black text,
%   axes, ticks, legends, titles, labels, and colorbars. Do not rely on
%   per-axes colour settings being preserved; the helper overrides them.

%% Figure dimensions (inches — MATLAB default unit)
s.FigWidthFull   = 6.3;    % full text width (~160 mm)
s.FigWidthHalf   = 3.1;    % half text width (~80 mm)
s.FigWidthDouble = 13.0;   % side-by-side (cold|warm subplot, full-width landscape)
s.FigHeight      = 3.5;    % default height
s.FigHeightTall  = 5.0;    % for 2-row subfigures
s.FigHeightWide  = 5.0;    % for FigWidthDouble layouts

%% Fonts (bumped 2026-05-13 so legends/callouts are legible at \textwidth)
s.FontName       = 'Helvetica';  % clean sans-serif; pairs well with LaTeX body
s.FontSizeAxis   = 13;           % axis labels
s.FontSizeTick   = 12;           % tick numbers
s.FontSizeLegend = 11;           % legend entries
s.FontSizeAnnot  = 11;           % annotations / text boxes
s.FontSizeTitle  = 14;           % panel titles (subplot headers)

%% Lines & markers
s.LineWidth      = 1.6;
s.LineWidthThin  = 0.9;          % individual run traces
s.MarkerSize     = 6;
s.MarkerSizeLg   = 8;

%% Background and text colours (thesis lock)
s.BackgroundColor = 'w';   % pure white
s.TextColor       = 'k';   % pure black
s.GridColor       = 'k';   % pure black grid (low alpha)
s.GridAlpha       = 0.15;

%% Colour palettes (RGB, 0-1)
% Camera counts: 6=blue, 7=amber, 8=teal
s.CameraColors   = [0.12 0.47 0.71;   % 6C
    0.85 0.53 0.10;   % 7C
    0.17 0.63 0.53];  % 8C

% Cost functions: CF1=blue, CF2=orange, CF3=green
s.CostFuncColors = [0.20 0.40 0.75;   % Resolution Uncertainty
    0.85 0.35 0.15;   % Dynamic Occlusion
    0.25 0.65 0.30];  % Combined

% Warm/Cold
s.WarmColor      = [0.80 0.25 0.20];
s.ColdColor      = [0.20 0.45 0.75];

% Target type: UAV=blue, UGV=orange
s.TargetColors   = [0.20 0.50 0.80;   % UAV
    0.90 0.55 0.20];  % UGV

% Grid mode: Uniform=green, Normal=purple
s.GridColors     = [0.30 0.65 0.35;   % Uniform
    0.55 0.30 0.65];  % Normal

% Shading alpha for envelopes
s.EnvelopeAlpha  = 0.20;

%% Label maps (index-matched to integer codes)
s.CostFuncNames  = {'Resolution Uncertainty', 'Dynamic Occlusion', 'Combined'};
s.CostFuncShort  = {'Res. Unc.', 'Dyn. Occ.', 'Combined'};
s.TargetNames    = {'UAV', 'UGV'};
s.GridNames      = {'Uniform', 'Normal'};

%% Export settings
s.ExportFormat   = 'pdf';        % vector output for LaTeX
s.ExportDPI      = 300;          % fallback if exporting raster
s.ExportBgColor  = 'white';      % exportgraphics BackgroundColor
end
