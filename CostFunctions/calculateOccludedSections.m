function Q = calculateOccludedSections(viewAngles, onAxis, camViewVectors, minAngle, maxAngle)
% Sum of occluder orientations [deg] with no triangulable front-side pair.
%
% The occluder is a vertical plane through the target point, described by
% its horizontal direction phi; cameras whose direction lies in
% (phi, phi + 180) are on its front side. Visibility only changes when the
% plane passes a camera's view vector or its extension, so the view lines
% split the circle into 2n sections of constant visibility (Rahimian &
% Kearney 2017, Sec. 3.4). Each section is tested at its middle orientation.

    lineAngles = viewAngles(~onAxis);
    if isempty(lineAngles)
        % Every camera is on the occluder axis: visibility never changes
        lineAngles = 0;
    end

    boundaries = unique(mod([lineAngles, lineAngles + 180], 360));
    boundaries = [boundaries, boundaries(1) + 360];

    Q = 0;
    for s = 1:numel(boundaries) - 1
        sectionSize = boundaries(s + 1) - boundaries(s);
        occluderAngle = mod(boundaries(s) + sectionSize/2, 360);

        frontCams = findFrontCameras(viewAngles, onAxis, occluderAngle);

        if ~checkSectionTriangulability(camViewVectors, frontCams, minAngle, maxAngle)
            Q = Q + sectionSize;
        end
    end
end
