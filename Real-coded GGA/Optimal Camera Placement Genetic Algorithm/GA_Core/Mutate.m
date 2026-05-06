function y = Mutate(x, mu, sigma)
    flag = (rand(size(x))< mu); %which indices mutate (logical operator) 
    y = x;
    r = randn(size(x));
    y(flag) = x(flag)+ sigma*r(flag); %guassian steps  
    %here it is possible that the values go beyond the feasible range of the variables
end 