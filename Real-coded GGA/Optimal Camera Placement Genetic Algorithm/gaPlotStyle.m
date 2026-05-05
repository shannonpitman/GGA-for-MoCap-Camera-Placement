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
s.FigWidthFull   = 6.3;   % full text width (~160 mm)
s.FigWidthHalf   = 3.1;   % half text width (~80 mm)
s.FigHeight      = 3.5;   % default height
s.FigHeightTall  = 5.0;   % for 2-row subfigures

%% Fonts
s.FontName       = 'Helvetica';  % clean sans-serif; pairs well with LaTeX body
s.FontSizeAxis   = 9;
s.FontSizeTick   = 8;
s.FontSizeLegend = 8;
s.FontSizeAnnot  = 7;           % annotations / text boxes

%% Lines & markers
s.LineWidth      = 1.4;
s.LineWidthThin  = 0.8;          % individual run traces
s.MarkerSize     = 5;
s.MarkerSizeLg   = 7;

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
