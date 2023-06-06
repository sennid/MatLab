function timevec = time1(nVec)
timevec = zeros(length(nVec),1);
vec = zeros(10,1);
for i = 1:length(nVec)
    n = nVec(i);
    for j = 1:10
        A = ceil(10*rand(n));
        B = ceil(10*rand(n));
        tic
        C = A*B;
        vec(j) = toc*1e6;
    end
    timevec(i) = median(vec);
end