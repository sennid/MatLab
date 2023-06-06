function integ = simpson(X, f)
h = X(2) - X(1);
a = X(1);
b = X(end);
X1 = a:(h/2):b;
Y = f(X1);
Y(isnan(Y)) = 1;

integ = h/6*(f(a)+f(b) + 2*sum(Y(2:2:end)) + 4*sum(Y(1:2:end-1)));

end