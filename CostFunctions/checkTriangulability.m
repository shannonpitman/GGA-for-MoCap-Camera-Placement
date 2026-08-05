function isTriangulable = checkTriangulability(viewVectors, idx1, idx2, minAngle, maxAngle)
%Check if two cameras satisfy triangulability constraint
%isTriangulable: Boolean indicating if pair is triangulable

    vec1 = viewVectors(:, idx1);
    vec2 = viewVectors(:, idx2);

    dotProduct = dot(vec1, vec2); %3D convergence angle
    
    % Valid range for acos
    dotProduct = max(-1, min(1, dotProduct));
    convergenceAngle = acosd(dotProduct); % degrees
    isTriangulable = (convergenceAngle >= minAngle) && (convergenceAngle <= maxAngle);
end