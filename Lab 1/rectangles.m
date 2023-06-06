function integ = rectangles(X, f)

h = X(2) - X(1);
Y = f(X);
Y(isnan(Y)) = 1;

integ = h*sum(Y);

end
