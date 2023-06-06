function timevec = time2(nVec)
timevec = zeros(length(nVec),1);
vec = zeros(5,1);
for k = 1:length(nVec)
    fprintf('progress = %.2f%% \n', k/length(nVec)*100);
    n = nVec(k);
    C = zeros(n);
    for l = 1:5
        A = ceil(10*rand(n));
        B = ceil(10*rand(n));
        tic
        for i = 1:n
            for j = 1:n
                C(i,j) = A(i,:)*B(:,j);
            end
        end
        vec(l) = toc*1e6;
    end
    timevec(k) = median(vec);
end