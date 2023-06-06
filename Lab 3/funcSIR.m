function f = funcSIR(t, x, beta1, beta2, alpha, gamma, delta, k)
    f = zeros(4,1);
    f(1) = -beta1.*x(3).*x(1) + k.*x(2);
    f(2) = -beta2.*x(3).*x(2) + alpha.*x(4) - k.*x(2);
    f(3) = beta1.*x(3).*x(1) + beta2.*x(3).*x(2) - gamma.*x(3);
    f(4) = gamma.*x(3) - delta.*x(4) - alpha.*x(4);
end