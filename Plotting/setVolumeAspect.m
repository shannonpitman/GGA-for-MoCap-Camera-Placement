function setVolumeAspect(ax, TargetSpace)
%SETVOLUMEASPECT  A 3D aspect ratio that works for both target modalities.
%
%   setVolumeAspect(ax, specs.Target)
%
%   The UAV target space is an 8 x 8 x 4 m volume; the UGV target space is
%   an 8 x 8 x 0.5 m floor slab. `axis equal` is right for the first and
%   unusable for the second — a true-scale slab renders as a line, leaves a
%   band of white above and below it in its tile, and pushes the panel
%   title up into the layout heading.
%
%   So: equal aspect when the volume is reasonably cubic, and a fixed plot
%   box with the z axis exaggerated when it is a slab. x and y stay equal
%   to each other either way, which is what matters for reading the
%   footprint; the z exaggeration is obvious from the tick labels.
%
%   See also plotHeatmap_GAvsOptiTrack, plotCostField_GAvsOptiTrack.

    r = range(TargetSpace, 1);          % [rx ry rz]
    rxy = max(r(1:2));

    if rxy > 0 && r(3) / rxy < 0.25
        % Slab: fix the plot box instead of the data scale.
        pbaspect(ax, [1 1 0.42]);
        daspect(ax, 'auto');
    else
        axis(ax, 'equal');
    end
end
