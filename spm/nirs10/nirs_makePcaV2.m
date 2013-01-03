function [d_PCA] = nirs_makePcaV2(d,channels,nSV)
d=d-repmat(mean(d,1),[size(d,1) 1]) ;     % substract the mean to have zero mean
[w s v] = svd(d);
for i0 = 1:nSV
    %Remove nSV eigenvalues -- svd will put the larger ones first
    s(i0,i0) = 0;
end
d_PCA = w * s * v';




