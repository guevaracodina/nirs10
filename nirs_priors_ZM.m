function [pE,pC] = nirs_priors_ZM(m,h)
% returns priors for a hemodynamic dynamic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m,[h])
% m   - number of inputs
% h   - number of hemodynamic modes (default = 3)
%
% pE  - prior expectations
% pC  - prior covariances
%
% (5) biophysical parameters
%    P(1) - signal decay                  d(ds/dt)/ds)
%    P(2) - autoregulation                d(ds/dt)/df)
%    P(3) - transit time                  (t0)
%    P(4) - exponent for Fout(v)          (alpha)
%    P(5) - resting oxygen extraction     (E0)
%   P(6) - time constant of vascular tone w               (tau_w)
%   P(7) - gain parameter b = b0 V0                       (b)
%   P(8) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   
%
%
% plus (m) efficacy priors
%    P(9) - ....
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

% biophysical parameters with prior expectation and
%note that vascular tone and gain parameters have been scaled down by 10
%---------------------------------------------------------------------------
pE    = [   0.65      0.41      0.98      0.32      0.34   0.668   0.760];

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
%i     = (h + 1):5;  %MODIFIER
%e(i)  = 0;
pC    = v*diag(e)*v';

%Append for vascular tone and gain parameters
pC = blkdiag(pC,diag([0.01 0.1]));

% append m efficacy priors
%---------------------------------------------------------------------------
pE    = [pE(:); zeros(m,1)];
pC    = blkdiag(pC,eye(m,m)*32);

