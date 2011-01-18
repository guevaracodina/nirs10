function [pE,pC] = nirs_hdm_priors_Huppert1(m,h)
% returns priors for a hemodynamic dynamic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m,[h])
% m   - number of inputs
% h   - number of hemodynamic modes (default = 3)
%
% pE  - prior expectations
% pC  - prior covariances
%
% (5) biophysical parameters
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(3) - transit time                                   (t0)
%   P(4) - exponent for Fout(v)                           (alpha)
%   P(5) - resting oxygen extraction                      (E0)
%   P(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   %not used here but used in measure
%          model spm_gx_hdm.m -- do not move this constant to other than
%          P(6)!
%According to Huppert et al, K is estimated from baseline values of SO2,
%flow, PtO2 and PvO2
%   P(7) - CMRO2 signal decay                             d(ds2/dt)/ds2)
%   P(8) - CMRO2 autoregulation                           d(ds2/dt)/dm)
%   P(9) - V0 Baseline volume fraction                    V0
%   P(10) - oxygen diffusion constant                      (K)
%
%   P(10 + 1:m)   - input efficacies                       d(ds/dt)/du)
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_priors.m 2050 2008-09-05 19:15:50Z klaas $



% default: 3 hemodynamic [eigen]modes
%---------------------------------------------------------------------------
if nargin < 2
    h = 3;
end

% biophysical parameters with prior expectation
%---------------------------------------------------------------------------
pE    = [   0.65      0.41      0.98      0.32      0.34];

% covariance restricted to h modes (v) scaled by eigenvalues (e) {see below)
%---------------------------------------------------------------------------
v     = [
        -0.0535    0.0095   -0.1117   -0.0040   -0.0026
        -0.0604   -0.0319    0.0430   -0.0077    0.0026
         0.1116   -0.0347   -0.2539   -0.0169   -0.0115
         0.1985    0.1698    0.4984   -0.4493    0.4434
         0.0029    0.2081    1.9582   -0.5209   -1.1634]';

e     = [2.1225    1.2006    0.3519    0.0039    0.0012];

% set variance of minor modes to zero
%---------------------------------------------------------------------------
i     = (h + 1):5;  %MODIFIER
e(i)  = 0;
pC    = v*diag(e)*v';

% append scaling parameter for epsilon: prior variance = 1/32,
%for this reason, this parameter must be the last one before the efficacies
%---------------------------------------------------------------------------
pE    = [pE 0];
pC    = blkdiag(pC,1/32);  %MODIFIER?

%append P(7) to P(10): 0.65  0.41  0.052 0.1
pE    = [pE  0.65  0.41  0.052 0.1];
%add arbitrary prior covariances for the extra parameters
pC    = blkdiag(pC,0.05, 0.03, 0.1, 0.1);

% append 2*m efficacy priors
%---------------------------------------------------------------------------
pE    = [pE(:); zeros(2*m,1)];
pC    = blkdiag(pC,eye(2*m,2*m)*32);

return