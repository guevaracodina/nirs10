function Q = interpolation_kernel(h1,W,LKC,UseCorrelRes)
try
    %Here we do 2 things
    %1- generate the interpolation kernels, B, Bx, By, B_volume
    %2- generate the interpolated beta and
    switch h1
        case 1
            ch = W.ch_HbO;
        case 2
            ch = W.ch_HbR;
        case 3
            ch = W.ch_HbT;
    end
    nch = length(ch);
    if ~W.Avg
        res = W.res(ch,:)';
        corr_beta = W.corr_beta;
        mtx_var = diag(W.var(ch)); %old NIRS_SPM version
    end
    s1 = W.s1;
    s2 = W.s2;
    %When using LKC, B, Bx, By are defined differently, and B_volume is
    %like the old B
    if nch > 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %1- generate the interpolation kernels, B, Bx, By, B_volume
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rchn = W.rchn;
        cchn = W.cchn;
        %identity over remaining channels
        mtx_eye = eye(nch);
        [x, y] = meshgrid(1:s2, 1:s1);
        if LKC || UseCorrelRes || W.Avg
            B_volume = zeros(s1, s2, nch);
        else
            B = zeros(s1, s2, nch);
            Bx = zeros(s1, s2, nch);
            By = zeros(s1, s2, nch);
        end
        for kk = 1:nch
            %try
            %    grid_eye = TriScatteredInterp(X, Y, V);
            %catch
            grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
            %end
            if kk == 1
                mask = 1 - isnan(grid_eye);
                if LKC || UseCorrelRes || W.Avg
                    index_mask0 = find(mask == 1);
                    L2 = zeros(1,nch);
                else
                    [rmask{1}, cmask{1}] = find(mask == 1);
                    index_mask0 = find(mask == 0);
                end
            end
            if LKC || UseCorrelRes || W.Avg
                B(kk, :) = grid_eye(index_mask0)';
                grid_eye(mask == 0) = 0;
                B_volume(:,:,kk) = grid_eye;
                [Bx_ch By_ch] = gradient(grid_eye);
                Bx(kk,:) = Bx_ch(index_mask0)';
                By(kk,:) = By_ch(index_mask0)';
                if ~W.Avg
                if kk == nch && LKC
                    L2 = calc_LKC(B_volume, mask, res,'individual');
                end
                end
            else
                grid_eye(index_mask0) = 0;
                [temp_Bx, temp_By] = gradient(grid_eye);
                temp_Bx(index_mask0) = 0;
                temp_By(index_mask0) = 0;
                B(:,:, kk) = grid_eye;
                Bx(:,:,kk) = temp_Bx;
                By(:,:,kk) = temp_By;
            end
        end
        Q.index_mask = index_mask0;
        if LKC || UseCorrelRes || W.Avg
            Q.B_volume = B_volume;
        else
            Q.rmask = rmask;
            Q.cmask = cmask;
        end
        if LKC
            L2 = sum(L2);
            r = sqrt(L2./pi);
            L1 = pi * r;
            L0 = 1;
            Q.L2 = L2;
            Q.L1 = L1;
            Q.L0 = L0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %2- generate the interpolated beta and
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculation of interp_beta and interp_var
        if LKC || UseCorrelRes
            Q.ibeta = W.beta(:,ch) * B;
            mtx_var = W.varch(ch, ch);
            %ensure mtx_var is symmetrical -- at machine-precision,
            %it is not due to round-off errors -- but that's not enough to
            %ensure that the eigenvalues are real
            mtx_var = (mtx_var+mtx_var')/2;
            [V_X D_X] = eig(mtx_var);
            tmp = D_X.^(1/2) * V_X' * B;
            Q.ivar = sum(tmp.^2,1);
            %         if ~isreal(Q.ivar)
            %             disp('Some values of ivar are imaginary -- check. Removing them for now');
            %             disp(['Min imaginary: ' num2str(min(imag(Q.ivar(:)))) ', max imaginary: ' num2str(max(imag(Q.ivar(:)))) ...
            %                 ', min real: ' num2str(min(real(Q.ivar(:)))) ', max real ' num2str(max(real(Q.ivar(:))))]);
            Q.ivar = real(Q.ivar); %Imaginary values are always very small
            %         end
        else %for tube only, for kappa calculation
            %Calculate covariance by SVD
            if ~W.Avg
                cov_beta = kron(mtx_var, corr_beta);
                [U, S, V] = svd(cov_beta);
                cov_beta_r = U*(S.^(0.5))*V';
                Q.cov_beta_r =cov_beta_r;
            else %for Avg
                Q.ibeta = W.beta(:,ch) * B;               
                Q.ivar =  W.covbeta(:,ch) * B; %Imaginary values are always very small
            end
        end
    else %not enough channels for interpolation
        disp('not enough channels for interpolation');
    end
    Q.B = B;
    Q.Bx = Bx;
    Q.By = By;
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end