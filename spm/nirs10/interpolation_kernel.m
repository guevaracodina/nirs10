function Q = interpolation_kernel(h1,W,LKC,UseCorrelRes,DoStats)
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
    if ~W.Avg && UseCorrelRes
        res = W.res(ch,:)';
        corr_beta = W.corr_beta;
        mtx_var = diag(W.var(ch)); %old NIRS_SPM version
    end
    
    %**********************************************************************
    %To have the pseudo-residuals in case of using the Avg method
    %Ke Peng, 2012-07-18
    %**********************************************************************
    if LKC && DoStats
        res = W.res(ch,:)';
    end
    %**********************************************************************
    
    s1 = W.s1;
    s2 = W.s2;
    %When using LKC, B, Bx, By are defined differently, and B_volume is
    %like the old B
    
    %warning('off') %this is to turn off the
    %Warning: Duplicate x-y data points detected: using average values for duplicate points
    %This is a problem that should be investigated, eventually...
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
            warning('off') %no big deal if some channels are superposed: the values will be averaged,
            %Unless one channel is bad, and one is good, and then the
            %average is useless...
            
            if W.no_interpolation
                grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'nearest');
                grid_eye_for_LKC = grid_eye;
            else
                grid_eye_for_LKC = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
                if W.AllowExtrapolation
                    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'v4');
                else
                    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
                end
            end
            warning('on')
            
            %end
            if kk == 1
                mask = 1 - isnan(grid_eye);
                mask_for_LKC = 1 - isnan(grid_eye_for_LKC);
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
                %**********************************************************
                %try to calculate LKC while using avg method
                %Ke Peng, 2012-07-18
                %**********************************************************
                if ~W.Avg || LKC
                    
                    %if ~W.Avg
                    if kk == nch && LKC && DoStats
                        L2 = calc_LKC(B_volume, mask_for_LKC, res,'individual');
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
        if LKC && DoStats
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
        
        %******************************************************************
        %Strict the condition
        %Ke Peng, 2012-07-18
        %******************************************************************
        %if LKC || UseCorrelRes
        if (LKC || UseCorrelRes) && ~W.Avg
            %******************************************************************
            %    if LKC || UseCorrelRes
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
                if DoStats
                    %Q.ivar =  W.covbeta(:,ch) * B; %Imaginary values are always very small
                    %Q.ivar = W.covbeta(:,ch) * reshape(Q.B_volume,[],nch)'; This formula is identical to the one above
                    for i0 = 1:size(W.covbeta,1)
                        Q.ivar(i0,:) = sum(B.* (diag(W.covbeta(i0,ch)) * B),1);
                    end
                    %Q.t = reshape(Q.ibeta(1,:) ./ (Q.ivar(1,:).^(0.5)),[362 434]) ;
                    %figure; imagesc(Q.t); colorbar
                end
            end
        end
    else %not enough channels for interpolation
        disp('not enough channels for interpolation');
    end
    %warning('on')
    Q.B = B;
    Q.Bx = Bx;
    Q.By = By;
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end