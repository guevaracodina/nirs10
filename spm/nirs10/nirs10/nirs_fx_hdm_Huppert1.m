function [f] = nirs_fx_hdm_Huppert1(x,u,P,M)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1) - vascular signal                                    s
%   x(2) - rCBF                                           log(f)
%   x(3) - venous volume                                  log(v)
%   x(4) - dHb                                            log(q)
%   x(5) - CMRO2 inducing signal                              s2
%   x(6) - m=CMRO2/CMRO20                                 log(m) 
%   x(7) - CtO2 (Tissue oxygen concentration)             log(CtO2)
%   x(8) - CvO2 (Venous oxygen concentration)             log(CvO2) 

% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
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
%   P(10 + 1:2*m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Huppert, Allen, Diamond, Boas 2008 Estimating CMRO2 from fMRI with a
% dynamic multicompartment Windkessel model
% Ref2: Huppert, Allen, Benav, Jones, Boas 2007 A multicompartment vascular
% model for inferring baseline and functional changes in CMRO2 and arterial
% dilation
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% LIOM - adapted from:

% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hdm.m 2495 2008-11-27 12:18:33Z karl $

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x([2 3 4 6 7 8]) = exp(x([2 3 4 6 7 8])); %excluding flow inducing signals

%define m, number of stimuli
m = length(u); %check!

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(3)^(1/P(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - P(5))^(1/x(2)))/P(5);

% implement differential state equations
%--------------------------------------------------------------------------
f(1)     = P(11:10+m)'*u(:) - P(1)*x(1) - P(2)*(x(2) - 1);
f(2)     = x(1)/x(2); 
f(3)     = (x(2) - fv/x(5))/(P(3)*x(3)); 
f(4)     = (ff*x(2) - fv*x(4)/x(3))/(P(3)*x(4));
f(5)     = P(11+m:10+2*m)'*u(:) - P(7)*x(5)-P(8)*(x(6)-1);
%x(5) and x(6): CMRO2 inducing signal modeled on equations for x(1) and x(2)
f(6)     = x(5)/x(6);
%dCtO2/dt = (1/V0) * (K * (CvO2 - CtO2) - CMRO2))
f(7)     = (P(10)*(x(8)-x(7))-x(6))/(P(9)*x(7));
%dCvO2/dt = (1/v) * ( f_in * CaO2 - f_out * CvO2 - CvO2 dv/dt)
%            -K * (CvO2-CtO2)
% take CaO2 = 95% as a constant
f(8)     = ( (x(2)*0.95-fv*x(8))/x(3) -x(8)*f(3)  - P(10)*(x(8)-x(7)))/x(8); %???
f        = f(:);

% adjust motion for DEM (that uses time-bins as units of time)
%--------------------------------------------------------------------------
try, global dt, f  = f*dt; end

return