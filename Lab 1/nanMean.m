function m = nanMean(X)
B = sum(~isnan(X));
A = sum(X, 'omitnan');
m = A./B;
