function [pE,pC] = spm_hdm_priors_YO(m)
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
%    P(6) - ratio of intra- to extra-     (epsilon)
%           vascular components of the
%           gradient echo signal   
%
% plus (m) efficacy priors
%    P(7) - ....
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_priors.m 4052 2010-08-27 19:22:44Z karl $



% default: 3 hemodynamic [eigen]modes
% %---------------------------------------------------------------------------
% if nargin < 2
%     h = 3;
% end
% h=5;
% biophysical parameters with prior expectation and
%---------------------------------------------------------------------------
%pE    = [   0.65      0.41      0.98      0.32      0.34  ];
%pE = [ 0.6564    0.4148    1.4404    0.3188    0.3378]; %    0.0268]; %estimated with BOLD+ASL
%pE = [ 0.6308    0.3685    1.5050    0.2898    0.3062]; %    0.0121]; %Session 2
%pE = [ 0.6436    0.3917    1.4727    0.3043    0.3220]; %    0.0195]; %Sessions Combined
% covariance restricted to h modes (v) scaled by eigenvales (e) {see below)
%---------------------------------------------------------------------------
% v     = [
%         -0.0535    0.0095   -0.1117   -0.0040   -0.0026
%         -0.0604   -0.0319    0.0430   -0.0077    0.0026
%          0.1116   -0.0347   -0.2539   -0.0169   -0.0115
%          0.1985    0.1698    0.4984   -0.4493    0.4434
%          0.0029    0.2081    1.9582   -0.5209   -1.1634]';
% 
% e     = 100*[2.1225    1.2006    0.3519    0.0039    0.0012];
% 
% % set variance of minor modes to zero
% %---------------------------------------------------------------------------
% i     = (h + 1):5;
% e(i)  = 0;
% pC    = v*diag(e)*v'/32;

% pC = [  0.0449    0.0005   -0.0017   -0.0000   -0.0000
%     0.0005    0.0060    0.0012    0.0000    0.0000
%    -0.0017    0.0012    0.1360   -0.0003   -0.0004
%    -0.0000    0.0000   -0.0003    0.0041    0.0000
%    -0.0000    0.0000   -0.0004    0.0000    0.0077];
% 
% % append scaling parameter for epsilon: prior variance = 1/32,
% %---------------------------------------------------------------------------
% %pE    = [pE 0];
% pE = [pE 0.0195];
% pC    = blkdiag(pC,1/32*10);
%pE = [ 0.6564    0.4148    1.4404    0.3188    0.3378    0.0268]; %BOLD+ASL %estimated with BOLD+ASL
pE = [ 0.6066    0.3876    1.5570    0.3295    0.3459   -0.1527]; %from BOLD
% pC = [   0.0449    0.0005   -0.0017   -0.0000   -0.0000    0.0005
%     0.0005    0.0060    0.0012    0.0000    0.0000   -0.0007
%    -0.0017    0.0012    0.1360   -0.0003   -0.0004    0.0044
%    -0.0000    0.0000   -0.0003    0.0041    0.0000    0.0002
%    -0.0000    0.0000   -0.0004    0.0000    0.0077    0.0001
%     0.0005   -0.0007    0.0044    0.0002    0.0001    0.3088]; %From
%     BOLD+ASL

  pC = [ 0.0380    0.0016   -0.0043   -0.0000    0.0000    0.0003
    0.0016    0.0052    0.0040    0.0002    0.0001   -0.0026
   -0.0043    0.0040    0.1224   -0.0011   -0.0009    0.0165
   -0.0000    0.0002   -0.0011    0.0040   -0.0000    0.0011
    0.0000    0.0001   -0.0009   -0.0000    0.0077    0.0007
    0.0003   -0.0026    0.0165    0.0011    0.0007    0.2947]; %From BOLD
% if modal == 2 %add a parameter for the relative strength of flow to BOLD
%      
% end

% append m efficacy priors
%---------------------------------------------------------------------------
%pE    = [pE(:); zeros(m,1)];
pE    = [pE(:); ones(m,1)]; %PP
%pC    = blkdiag(pC,eye(m,m)*2e-7);
pC    = blkdiag(pC,eye(m,m));

return


% NOTES: sample covariances from Friston et al (2000)
%---------------------------------------------------------------------------
qC    = [0.0150    0.0052    0.0283    0.0002   -0.0027
         0.0052    0.0020    0.0104    0.0004   -0.0013
         0.0283    0.0104    0.0568    0.0010   -0.0069
         0.0002    0.0004    0.0010    0.0013   -0.0010
        -0.0027   -0.0013   -0.0069   -0.0010    0.0024];


% NOTES: Reduce rank of prior covariances for computational expediancy
%---------------------------------------------------------------------------

% assume independent priors in parameter space
%---------------------------------------------------------------------------
qC    = diag(diag(qC));


% model specification (single node DCM)
%---------------------------------------------------------------------------
M.f   = 'spm_fx_hdm';
M.g   = 'spm_gx_hdm';
M.x   = [0 0 0 0]';
M.pE  = [pE(1:6) 1];
M.m   = 1;
M.n   = 4;
M.l   = 1;
M.N   = 32;
M.dt  = 1/2;

% compute partial derivatives w.r.t. hemodynamic parameters [J] dy(t)/dp
%---------------------------------------------------------------------------
P     = M.pE;
p     = length(P);
dp    = 1e-6;
[k J] = nirs_nlsi(M);
for i = 1:5
    M.pE    = P;
    M.pE(i) = M.pE(i) + dp;
    [k q]   = nirs_nlsi(M);
    Jq(:,i) = (q - J)/dp;
end

% implied covariance of impulse response
%---------------------------------------------------------------------------
Cq    = Jq*qC*Jq';

% reduce to h hemodynamic modes in measurement space
%---------------------------------------------------------------------------
[v e] = spm_svd(Cq);
e     = diag(e);
v     = pinv(Jq)*v;
qC    = v*diag(e)*v';


% NOTES: graphics - eigenvalues of qC
%---------------------------------------------------------------------------
subplot(2,2,1)
bar(e)
xlabel('eigen mode')
title('eigenvalue')
set(gca,'XLim',[0 6])
axis square
grid on

% graphics - response differentials
%---------------------------------------------------------------------------
subplot(2,2,2)
plot([1:M.N]*M.dt,Jq*v(:,1),[1:M.N]*M.dt,Jq*v(:,2),'-.')
xlabel('PST {secs}')
title('hemodynamic modes')
axis square
grid on
