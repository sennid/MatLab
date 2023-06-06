function dx = func10(t,x)
    
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -x(1) + 2*t - pi;
end