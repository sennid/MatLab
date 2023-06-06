function f = func6(t, x, alpha)
    f = zeros(4,1);
    f(1) = x(2);
    f(2) = -alpha .* x(1);
    f(3) = x(4);
    f(4) = -alpha .* x(3);
end