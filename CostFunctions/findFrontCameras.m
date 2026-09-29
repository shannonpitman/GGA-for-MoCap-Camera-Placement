function frontCams = findFrontCameras(viewAngles, onAxis, occluderAngle)
%FINDFRONTCAMERAS Cameras on the front side of a vertical occluder
%   The occluder plane runs along occluderAngle [deg]; its front side holds
%   directions in (occluderAngle, occluderAngle + 180). Cameras on the
%   occluder axis (onAxis) are never occluded. Indices refer to the columns
%   of viewAngles, i.e. the same order as the view vectors.
    angleDiffs = mod(viewAngles - occluderAngle, 360);
    frontCams = find(onAxis | (angleDiffs > 0 & angleDiffs < 180));
end
