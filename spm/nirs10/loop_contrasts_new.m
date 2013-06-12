function C = loop_contrasts_new(h1,Q,W,xCon,Z,f1,cSigma0,xX)
try
    switch h1
        case 1
            ch = W.ch_HbO;
        case 2
            ch = W.ch_HbR;
        case 3
            ch = W.ch_HbT;
    end
    s1 = W.s1;
    s2 = W.s2;
    B = Q.B;
    if ~W.Avg
        nC = size(xCon,2);
    else
        nC = size(W.beta,1);
    end
    if Z.DoStats
        stat_map = zeros(nC,s1,s2); %T-stat
        beta_map = zeros(nC,s1,s2);
        cov_map = zeros(nC,s1,s2);
        for c1 = 1:nC
            stat_tmp = zeros(s1,s2);
            beta_tmp = zeros(s1,s2);
            cov_tmp = zeros(s1,s2);
            if ~W.Avg
                c = xCon(c1).c;
            else
                c = zeros(nC+1,1);
                c(c1,1) = 1;
                xCon(c1).STAT = 'T';
            end
            if xCon(c1).STAT == 'T'
                % covariance of interpolated beta
                if ~W.Avg
                    cxCor = c' * W.corr_beta * c;
                    stat = (c' * Q.ibeta)./sqrt(cxCor .* Q.ivar);
                    beta_tmp(Q.index_mask) = c' * Q.ibeta;
                    cov_tmp(Q.index_mask) = cxCor .* Q.ivar;
                else
                    c = c(1:end-1);
                    Q.ivar(isnan(Q.ivar)) = max(Q.ivar(~isnan(Q.ivar)));
                    Q.ibeta(isnan(Q.ibeta)) = 0;
                    Q.ivar(Q.ivar < 0) = -Q.ivar(Q.ivar < 0);
                    stat = (c' * Q.ibeta(1:nC,:))./sqrt(c' * Q.ivar(1:nC,:));
                    beta_tmp(Q.index_mask) = c' * Q.ibeta(1:nC,:);
                    cov_tmp(Q.index_mask) = c' * Q.ivar(1:nC,:);
                end
                beta_map(c1,:,:) = beta_tmp;
                cov_map(c1,:,:) = cov_tmp;
                if ~isreal(stat)
                    stat = real(stat);
                    disp('Some imaginary values removed in t-stat calculation');
                end
                stat_tmp(Q.index_mask) = stat;
                stat_map(c1,:,:) = stat_tmp;
            else
                ResSS = xX.ResSSch(ch, ch);
                cSigma = cSigma0{f1,c1}(ch,ch);
                
                [V_X D_X] = eig(ResSS);
                [V_X0 D_X0] = eig(cSigma);
                tmp = D_X.^(1/2) * V_X' * B;
                ip_ResSS = sum(tmp.^2,1);
                tmp = D_X0.^(1/2) * V_X0' * B;
                ip_ResSS_X0 = sum(tmp.^2,1);
                F_stat = ((ip_ResSS_X0 - ip_ResSS)./xCon(c1).trRV)./(ip_ResSS./xX.trRV);
                if ~isreal(F_stat)
                    F_stat = real(F_stat);
                    disp('Some imaginary values removed for F_stat');
                end
                stat_tmp(Q.index_mask) = F_stat;
                stat_map(c1,:,:) = stat_tmp;
                beta_map(c1,:,:) = stat_tmp;
            end
        end
        if Z.spatial_LPF %does not work properly -- do not use
            K.k1 = s1;
            K.k2 = s2;
            K.radius = Z.radius;
            K = spatial_LPF('set',K);
            stat_map = spatial_LPF('lpf',K,stat_map);
        end
        C.stat_map = stat_map;
        C.beta_map = beta_map; %Used at the group level
        C.cov_map = cov_map;
    else %simple interpolation of the betas
        if nC == 1 && size(Q.ibeta,1) == 1
            beta_map = zeros(1,s1,s2);
            beta_tmp = zeros(s1,s2);
            Q.ibeta(isnan(Q.ibeta)) = 0;
            beta_tmp(Q.index_mask) = Q.ibeta;
            beta_map(1,:,:) = beta_tmp;
            C.beta_map = beta_map; %Used at the group level
        else
            nC = size(Q.ibeta,1);
            %beta_map = zeros(nC,s1,s2);
            beta_tmp = zeros(nC,s1,s2);
            Q.ibeta(isnan(Q.ibeta)) = 0;
            beta_tmp(:,Q.index_mask) = Q.ibeta(1:nC,:);
            %beta_map = beta_tmp;
            C.beta_map = beta_tmp; %Used at the group level
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end