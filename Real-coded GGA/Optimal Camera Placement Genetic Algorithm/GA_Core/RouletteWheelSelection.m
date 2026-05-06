function indices = RouletteWheelSelection(p, n)
    c = [0; cumsum(p(:))];  
    r = rand(n, 1) * c(end);
    indices = discretize(r, c);
end