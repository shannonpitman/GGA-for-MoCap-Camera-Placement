function b = costBreakdown(chrom, specs)
%COSTBREAKDOWN  Raw and weighted CF3 components of one chromosome.
%
%   b = costBreakdown(chrom, specs)
%
%   b.RawUnc, b.RawOcc  raw resolution-uncertainty / dynamic-occlusion costs
%   b.Junc, b.Jocc, b.J weighted, utopia-shifted CF3 terms (cf3Terms), on the
%                       specs' normalisation table; J = Junc + Jocc
    [cams, cc] = setupCameras(chrom, specs.Cams, specs.Resolution, specs.Focal, ...
        specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);
    b.RawUnc = resUncertainty(specs, cams, cc);
    b.RawOcc = dynamicOcclusion(specs, cams, cc);
    [b.J, b.Junc, b.Jocc] = cf3Terms(b.RawUnc, b.RawOcc, specs);
end
