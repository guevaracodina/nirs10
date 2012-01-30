function C = loop_contrasts(h1,Q,W,xCon,Z)
C = [];
switch h1
    case 1
        ch = W.ch_HbO;
    case 2
        ch = W.ch_HbR;
    case 3
        ch = W.ch_HbT;
end
mtx_var = diag(W.var(ch));
beta = W.beta(:,ch);
beta = beta(:);
corr_beta = W.corr_beta;
cov_beta_r = Q.cov_beta_r;
rmask = Q.rmask;
cmask = Q.cmask;
s1 = W.s1;
s2 = W.s2;
rmv = rmask{1};
cmv = cmask{1};
nC = size(xCon,2);
nm = length(rmv);
%preallocate
sum_kappa = zeros(nC,1);
kappa = zeros(nC,s1,s2);
c_interp_beta = zeros(nC,s1,s2); %analog of con0001.nii in SPM
c_cov_interp_beta = zeros(nC,s1,s2);
c_interp_T = zeros(nC,s1,s2); %T-stat
c_interp_ess = zeros(nC,s1,s2); %analog of ess0001.nii in SPM
c_interp_ess0 = zeros(s1,s2); %extra-sum of squares of full model (ResMS)
c_interp_F = zeros(nC,s1,s2); %F-stat
sz_xCon  = size(xCon(1).c,1);
%identity matrix of size number of regressors
tmp = eye(sz_xCon);
for kk = 1:nm
    %this is different for HbO and HbR
    B2(:,1) = Q.B(rmv(kk), cmv(kk), :);
    B2x(:,1) = Q.Bx(rmv(kk), cmv(kk), :);
    B2y(:,1) = Q.By(rmv(kk), cmv(kk), :);
    B3 = kron(B2, tmp);
    B3x = kron(B2x, tmp);
    B3y = kron(B2y, tmp);
    B3t = kron(B2', tmp);
    d = (B2'*mtx_var*B2); %ResSS/TrRV
    c_interp_ess0(rmv(kk), cmv(kk)) = d; %ResSS/TrRV -- recall normalization of ResSS by TrRV is included here (while it is a scale factor in the SPM nifti)
    for c1 = 1:nC
        c = xCon(c1).c;
        if xCon(c1).STAT == 'T'
            %this is the same for HbO and HbR
            c_corr_beta = c' * corr_beta * c;
            kappa(c1,rmv(kk),cmv(kk)) = calc_kappa(B3,B3x,B3y,cov_beta_r,c);
            c_interp_beta(c1,rmv(kk), cmv(kk)) = (c' * B3t) * beta;
            c_cov_interp_beta(c1,rmv(kk), cmv(kk)) = (B2'*mtx_var*B2) * c_corr_beta;
            c_interp_T(c1,rmv(kk), cmv(kk)) = c_interp_beta(c1,rmv(kk), cmv(kk))/(c_cov_interp_beta(c1,rmv(kk), cmv(kk)))^0.5;           
        else
            %This is an F-stat; if xCon(c1).STAT == 'F'
            %Do a GLM for the reduced model for each F contrast
            hsqr = xCon(c1).h * B3t;
            %Numerator of F-test
            c_ResSS =  (beta'*hsqr')*(hsqr*beta);
            %Interpolate
            c_interp_F(c1,rmv(kk), cmv(kk)) = c_ResSS/(d*xCon(c1).trRV);
            c_interp_ess(c1,rmv(kk), cmv(kk)) = c_ResSS/xCon(c1).trRV; %note that we also normalize here by eidf
        end
    end
end
for c1 = 1:nC
    if xCon(c1).STAT == 'T'
        tm = kappa(c1,:,:);
        sum_kappa(c1) = sum(tm(:));
    end
end

if Z.spatial_LPF %does not work properly -- do not use
    K.k1 = s1;
    K.k2 = s2;
    K.radius = W.radius;
    K = spatial_LPF('set',K);
    c_interp_beta = spatial_LPF('lpf',K,c_interp_beta);
    c_cov_interp_beta = spatial_LPF('lpf',K,c_cov_interp_beta);
    %should interpolate all the rest too...
end
C.sum_kappa = sum_kappa;
C.c_interp_beta = c_interp_beta;
C.c_cov_interp_beta = c_cov_interp_beta;
C.c_interp_T = c_interp_T;
C.c_interp_F = c_interp_F;
C.c_interp_ess0 = c_interp_ess0;
C.c_interp_ess = c_interp_ess;