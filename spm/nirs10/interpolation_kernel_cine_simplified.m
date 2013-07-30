function Q = interpolation_kernel_cine_simplified(W)
nch = length(W.ch_HbO);
s1 = W.s1;
s2 = W.s2;
rchn = W.rchn;
cchn = W.cchn;
%identity over remaining channels
mtx_eye = eye(nch);
[x, y] = meshgrid(1:s2, 1:s1);
for kk = 1:nch
    warning('off') %no big deal if some channels are superposed: the values will be averaged,
    %Unless one channel is bad, and one is good, and then the
    %average is useless...
    if W.no_interpolation
        grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'nearest');
    else
        if W.AllowExtrapolation
            grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'v4');
        else
            grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
        end
    end
    warning('on')
    if kk == 1
        mask = 1 - isnan(grid_eye);
        index_mask0 = find(mask == 1);
        B = zeros(nch,length(index_mask0));
    end
    B(kk, :) = grid_eye(index_mask0)';
    Q.index_mask = index_mask0;
end
Q.B = B;
