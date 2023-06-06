function timevec = time(f, nVec)
timevec = zeros(length(nVec),1);
vec = zeros(10,1);
for i = 1:length(nVec)
    n = nVec(i);
    for j = 1:10
        tic
        f(X, Y)
        vec(j) = toc*1e6;
    end
    timevec(i) = median(vec);
end