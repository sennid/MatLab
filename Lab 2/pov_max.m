function pos = pov_max(fMat)
    fMat1 = zeros(size(fMat));
    fMat1(1,:) = fMat(1,:);
    fMat1(2:end,:) = fMat(1:end-1,:);

    fMat2 = zeros(size(fMat));
    fMat2(end,:) = fMat(end,:);
    fMat2(1:end-1,:) = fMat(2:end,:);

    pos1 = find((fMat>=fMat1)&(fMat>=fMat2));

    fMat1 = zeros(size(fMat));
    fMat1(:,1) = fMat(:,1);
    fMat1(:,2:end) = fMat(:,1:end-1);

    fMat2 = zeros(size(fMat));
    fMat2(:,end) = fMat(:,end);
    fMat2(:,1:end-1) = fMat(:,2:end);

    pos2 = find((fMat>=fMat1)&(fMat>=fMat2));

    pos = intersect(pos1,pos2);
    
end