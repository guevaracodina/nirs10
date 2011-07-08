function Q = interpolation_kernel_short(nch,W)
s1 = W.s1;
s2 = W.s2;
rchn = W.rchn;
cchn = W.cchn;
%identity over remaining channels
mtx_eye = eye(nch);
%mesh of topographically projected brain size; note s1 and s2 inverted
[x, y] = meshgrid(1:s2, 1:s1);

B = zeros(s1, s2, nch);
for kk = 1:nch
    %a grid for channel positions?
    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic'); %could try nearest?
    if  kk == 1
        %mask the NaN
        mask = isnan(grid_eye);
        mask = 1- mask;
        [rmask{1}, cmask{1}] = find(mask == 1);
        index_mask0 = mask == 0;
    end
    grid_eye(index_mask0) = 0;
    B(:,:, kk) = grid_eye;
end
Q.B = B;
Q.rmask = rmask;
Q.cmask = cmask;
end