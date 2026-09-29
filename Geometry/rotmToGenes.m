function r = rotmToGenes(R)
%ROTMTOGENES  Rotation vector (exponential map, |r| <= pi) from a rotation matrix.
%   Inverse of genesToRotm. Handles the theta ~ pi case, where the axis
%   comes from the symmetric part of R.
    c = max(-1, min(1, (trace(R) - 1) / 2));
    theta = acos(c);
    if theta < 1e-9
        r = [R(3,2) - R(2,3), R(1,3) - R(3,1), R(2,1) - R(1,2)] / 2;
        return;
    end
    if pi - theta > 1e-6
        k = [R(3,2) - R(2,3), R(1,3) - R(3,1), R(2,1) - R(1,2)] / (2*sin(theta));
    else
        % theta ~ pi: R = 2kk' - I, take the best-conditioned column
        B = (R + eye(3)) / 2;
        [~, j] = max(diag(B));
        k = B(:, j).' / sqrt(max(B(j, j), eps));
        % fix sign from the (small) antisymmetric part when available
        s = [R(3,2) - R(2,3), R(1,3) - R(3,1), R(2,1) - R(1,2)];
        if dot(k, s) < 0, k = -k; end
    end
    r = rotvecWrap(k / norm(k) * theta);
end
