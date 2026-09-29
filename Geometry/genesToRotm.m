function R = genesToRotm(r)
%GENESTOROTM  Camera-to-world rotation from a camera's orientation genes.
%
%   R = genesToRotm(r) with r = theta*k (1x3 rotation vector / exponential
%   map: axis k, angle theta in rad). Rodrigues' formula:
%       R = I + sin(theta) K + (1 - cos(theta)) K^2,   K = skew(k)
%   The single place orientation genes are turned into a rotation matrix,
%   so the parameterisation can never disagree between files.
    r = r(:);
    theta = norm(r);
    if theta < 1e-12
        R = eye(3) + skew(r);           % first-order, exact as theta -> 0
        return;
    end
    K = skew(r / theta);
    R = eye(3) + sin(theta)*K + (1 - cos(theta))*(K*K);
end

function S = skew(v)
    S = [    0, -v(3),  v(2);
          v(3),     0, -v(1);
         -v(2),  v(1),     0];
end
