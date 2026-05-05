function [y1,y2] = UniformCrossover(x1,x2, xFace1, xFace2, gamma)
    
    % alpha = rand(size(x1));
    alpha = unifrnd(-gamma, 1+gamma, size(x1)); %improves exploration
    y1 = alpha.*x1+ (1-alpha).*x2;
    y2 = alpha.*x2 + (1-alpha).*x1; %here it is possible that the values go beyond the feasible range of the variables
end