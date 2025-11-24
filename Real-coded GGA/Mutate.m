function y = Mutate(x, flags, sigma)
    % Optimized mutation that only modifies flagged genes
    y = x;
    if any(flags)
        r = randn(1, length(x));
        y(flags) = x(flags) + sigma * r(flags);
    end
end