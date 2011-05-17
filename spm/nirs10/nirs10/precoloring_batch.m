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
SPM.KY = KY;
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

if SPM.generate_trRV
    [trRV trRVRV ] = approx_trRV(SPM.xX.xKXs.X,SPM.xX.pKX,S);
else
    trRV = 0;
    trRVRV = 0;
end
SPM.xX.trRV = trRV; % <R'*y'*y*R>
SPM.xX.trRVRV = trRVRV; %- Satterthwaite
try SPM.xX.erdf = trRV^2/trRVRV; end
% SPM.xX.Bcov = SPM.xX.pKX*V*SPM.xX.pKX';
SPM.xX.Bcov = (SPM.xX.pKX * S);
SPM.xX.Bcov = SPM.xX.Bcov * SPM.xX.Bcov';
%SPM.nirs.step = 'estimation';
end