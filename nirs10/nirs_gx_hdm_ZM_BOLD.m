function [y] = nirs_gx_hdm_ZM_BOLD(x,u,P,M)
% Simulated BOLD response to input.
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% y(1)    - BOLD response (%)
% y(2)    - ASL flow
% x    - state vector     (see spm_fx_dcm)
% P    - Parameter vector (see spm_fx_dcm)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_hdm.m 3666 2010-01-10 17:54:31Z klaas $


% biophysical constants for 1.5 T: 
%==========================================================================

% echo time (seconds)
%--------------------------------------------------------------------------
try
    TE = M(1).TE;
catch
    TE = 0.03; %MODIFIER
end

%resting venous volume
%--------------------------------------------------------------------------
V0    = 100*0.08;                                

% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
r0 = 100; %3T r0    = 25; %r0 = 100; %Huppert 2008 -- Mildner 2001

% frequency offset at the outer surface of magnetized vessels
%--------------------------------------------------------------------------
nu0 = 80.6; %3T nu0   = 40.3; %nu0 = 80.6; %Huppert Mildner

% region-specific resting oxygen extraction fractions
%-------------------------------------------------------------------------- 
E0    = P(5); 

% region-specific ratios of intra- to extravascular components of
% the gradient echo signal (prior mean = 1, log-normally distributed 
% scaling factor)
%--------------------------------------------------------------------------
epsi  = exp(P(8));    %MODIFIER
 
% coefficients in BOLD signal model
%--------------------------------------------------------------------------
k1    = 4.3.*nu0.*E0.*TE;
k2    = epsi.*r0.*E0.*TE;
k3    = 1 - epsi;
 
% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x     = exp(x); 

% BOLD signal
%--------------------------------------------------------------------------
v     = x(3);
q     = x(4);
switch M.modal
    case 1
        y     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v));
    case 2
        cf = 1.;
        y(2) = cf*(x(2)-1); %flow 
        y(1)     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v)); 
    case 3
        cf = 1.;
        y = cf*(x(2)-1); %flow
    case 4 %HbO+HbR
        cf1 = 1.; cf2 = 1.;
        y(1) = cf1*(x(3)-1); %Assume v = HbT = HbO+ HbR
        y(2) = cf2*(x(4)-1); %Assume q = HbR
end
y = y(:);


