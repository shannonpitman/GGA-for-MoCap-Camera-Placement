function y = Mutate(x, mu, sigma)
% Gaussian mutation. Each gene mutates with probability mu.
% sigma is a scalar, or a 1x6 per-camera vector [pos pos pos rot rot rot]
% (metres and radians) that is repeated for every camera.
    if numel(sigma) == 6 && numel(x) > 6
        sigma = repmat(sigma(:).', 1, numel(x)/6);
    end
    if isscalar(sigma)
        sigma = repmat(sigma, size(x));
    end
    flag = (rand(size(x))< mu); %which indices mutate (logical operator) 
    y = x;
    r = randn(size(x));
    y(flag) = x(flag)+ sigma(flag).*r(flag); %guassian steps  
    %here it is possible that the values go beyond the feasible range of the variables
end 
