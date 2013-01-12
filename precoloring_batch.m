function [SPM res] = precoloring_batch(varargin)
% NIRS_SPM: this function is used for estimation of GLM parameters.
SPM = varargin{1};
Y = varargin{2};
nScan = size(SPM.xX.X,1);

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
%GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
%SPM.VResMS if we had full images as with fMRI
SPM.xX.beta = SPM.xX.pKX * KY; % beta : least square estimate
%Residuals are grouped as [HbO; HbR; HbT] and saved for each session
res = spm_sp('r', SPM.xX.xKXs, KY); % Residuals
% update for calculating the channel-wise least-square residual correlation
% date: Aug 10, 2011
SPM.xX.ResSS = sum(res.^2); %old residuals
ResSS = (KY' * KY) - (KY' * SPM.xX.xKXs.X) * SPM.xX.pKX * KY; %channel by channel
SPM.xX.ResSSch = ResSS; % Residual sum of squares
%data no longer required
SPM.KY = KY;
res = res'; %residuals now as channels by time, ready to save as .nii
clear KY Y 

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
if ~isfield(SPM,'TrRVRVexact')
    SPM.TrRVRVexact = 0; %approximation
end

if SPM.generate_trRV
    [trRV trRVRV ] = approx_trRV(SPM.xX.xKXs.X,SPM.xX.pKX,S,[],SPM.TrRVRVexact,SPM);
else
    trRV = 0;
    trRVRV = 0;
end
SPM.xX.trRV = trRV; % <R'*y'*y*R>
SPM.xX.trRVRV = trRVRV; %- Satterthwaite
try 
    SPM.xX.erdf = trRV^2/trRVRV; 
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem calculating degrees of freedom');    
end
% SPM.xX.Bcov = SPM.xX.pKX*V*SPM.xX.pKX';
SPM.xX.Bcov = (SPM.xX.pKX * S);
SPM.xX.Bcov = SPM.xX.Bcov * SPM.xX.Bcov';