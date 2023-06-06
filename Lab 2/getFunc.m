function T = getFunc(n)
    T = zeros(n+1);
    T(1, 1) = 1;
    T(2, 2) = 1;
    for i = 3:(n+1)
        T(i, :) = [0, 2*T(i-1, 1:end-1)] - T(i-2, :);
    end
end