function F_stat = nirs_calculate_F_contrasts(SP,Y,c)
SPM = SP.SPM;

% load parameters from SPM_nirs structures
X = SPM.xX.xKXs.X; % design matrix
ResSS = SPM.xX.ResSS;
nScan = size(SPM.xX.X,1);
% generation of S matrix (filtering) for calculation of trMV & trMVMV
switch SPM.xX.K.LParam.type
    case {'hrf', 'Gaussian'}
        S = SPM.xX.K.KL;
    case 'none'
        S = speye(nScan);
end

nCont = size(X, 2);
c0 = eye(nCont) - c * pinv(c);
X0 = X * c0;
cSigma = (KY' * KY) - (KY' * X0) * pinv(X0) * KY;

% calculation of trace(MV) & trace(MVMV)
[trMV trMVMV ] = approx_trRV(X,pKX,S,c);
% end of calculation of trace MV and trace MVMV
df(1) = trMVMV ./ (trMV.^2);
df(2) = SPM.xX.erdf;


F_stat = [];
% interpolating matrix dependent-term calculation
nch = length(chs{kk});
ResSS_set = zeros(nch);
cSigma_set = zeros(nch);
for aa = 1:nch
    for bb = 1:nch
        ResSS_set(aa,bb) = ResSS(chs{kk}(aa), chs{kk}(bb));
        cSigma_set(aa,bb) = cSigma(chs{kk}(aa), chs{kk}(bb));
    end
end
[V_X D_X] = eig(ResSS_set);
[V_X0 D_X0] = eig(cSigma_set);
tmp = D_X.^(1/2) * V_X' * B;
ip_ResSS = sum(tmp.^2,1);
tmp = D_X0.^(1/2) * V_X0' * B;
ip_ResSS_X0 = sum(tmp.^2,1);
F_stat = ((ip_ResSS_X0 - ip_ResSS)./trMV)./(ip_ResSS./SPM.xX.trRV);
  

