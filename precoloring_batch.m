function [SPM] = precoloring_batch(varargin)
% NIRS_SPM: this function is used for estimation of GLM parameters.
SPM = varargin{1};
Y = varargin{2};
overrideJ0 = 0;
if size(varargin) == 3
    J0 = varargin{3};
    overrideJ0 = 1;
end
[nScan nBeta] = size(SPM.xX.X);

% Set up of the filters / Design space and projector matrix [pseudoinverse] for WLS
%==================================================================
switch SPM.xX.K.HParam.type
    case 'Wavelet-MDL'
        tmp_K = SPM.xX.K;
        tmp_K.HParam.type = '';
        SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(tmp_K, SPM.xX.X)); % KX
        SPM.xX.K.X = SPM.xX.X;
        clear tmp_K;
    case 'DCT'
        SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.xX.X)); % KX
    case 'none'
        SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.xX.X)); % KX ?
end

SPM.xX.xKXs.X = full(SPM.xX.xKXs.X);
SPM.xX.pKX = spm_sp('x-', SPM.xX.xKXs); % projector

%filtering of the data
KY = spm_filter_HPF_LPF_WMDL(SPM.xX.K, Y);
%For testing
%figure; nn = 8; plot(Y(:,nn)); hold on; plot(KY(:,nn),'r'); hold off

%GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
%SPM.VResMS if we had full images as with fMRI
SPM.xX.beta = SPM.xX.pKX * KY; % beta : least square estimate
res = spm_sp('r', SPM.xX.xKXs, KY); % Residuals 
SPM.xX.ResSS = sum(res.^2); % Residual SSQ
%data no longer required
clear KY Y res

%more filtering operations
switch SPM.xX.K.LParam.type
    case {'hrf', 'Gaussian'}
        S = SPM.xX.K.KL;
    case 'none'
        S = speye(nScan);
end

switch SPM.xX.K.HParam.type 
    case 'DCT'
         S = S - SPM.xX.K.X0 * (SPM.xX.K.X0' * S);                                    
        %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
end

%Difficulty with calculation of TrRV and especially TrRVRV is that the
%nScan values that contribute to the traces are very badly distributed,
%with a few very large values. Strategy is to calculate the values that
%contribute to TrRV exactly, and since the indices where they occur will
%be the same as for trRVRV, approximate trRVRV by those values.

%computation of var1 and var2 is fairly quick
var1 = S * S';
var2 = var1 - SPM.xX.xKXs.X * (SPM.xX.pKX * var1); % RV * e(kk)
%var3 = var1*var2 is the long calculation - instead...
%find extreme values (k2) of var2
RVd = diag(var2);
trRV = sum(RVd); %exact
RVs = std(RVd);
RVm = median(RVd);
k2 = abs(RVd-RVm)>RVs; %exceeding one standard deviation 
%largest contributions to var3
var3e = var1 * var2(:,k2); %e for extreme
var4e = var3e - SPM.xX.xKXs.X * (SPM.xX.pKX * var3e);
RVRVe = diag(var4e(k2,:));
%find typical value of RVRV contribution:
sample_size = 100;
var3t = var1 * var2(:,1:sample_size); %t for typical
var4t = var3t - SPM.xX.xKXs.X * (SPM.xX.pKX * var3t);
RVRVt = median(diag(var4t(1:sample_size,:))); %median is insensitive to extreme values
%trRVRV approximated by sum of extreme values and contribution from typical
%values
trRVRV = sum(RVRVe) + RVRVt*(nScan-length(find(k2)));

%In one example, "exact" calculation of trRV, trRVRV gave:
% 341.25, 238.84
%while truncation to 100 random samples gave
% 334.19 , 202.19
%and trunction to 1000 random samples gave
% 339.53, 206.81
%so the error is as much as 15%

SPM.xX.trRV = trRV; % <R'*y'*y*R>
SPM.xX.trRVRV = trRVRV; %- Satterthwaite
SPM.xX.erdf = trRV^2/trRVRV;
% SPM.xX.Bcov = SPM.xX.pKX*V*SPM.xX.pKX';
SPM.xX.Bcov = (SPM.xX.pKX * S);
SPM.xX.Bcov = SPM.xX.Bcov * SPM.xX.Bcov';
%SPM.nirs.step = 'estimation';

try
    K = SPM.xX.K;
    K = rmfield(K, 'X');
    K = rmfield(K, 'KL');
    SPM.xX.K = K;
    %clear K;
end
end