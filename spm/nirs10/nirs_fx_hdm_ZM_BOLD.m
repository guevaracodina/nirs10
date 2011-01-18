function [f] = nirs_fx_hdm_ZM_BOLD(x,u,P,M)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1) - vascular signal                                    s
%   x(2) - rCBF                                           log(f)
%   x(3) - venous volume                                  log(v)
%   x(4) - dHb                                            log(q)
%   x(5) - "w=A/A_SS": normalized vascular tone           log(w)
% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(3) - transit time                                   (t0)
%   P(4) - exponent for Fout(v)                           (alpha)
%   P(5) - resting oxygen extraction                      (E0)
%   P(6) - time constant of vascular tone w (divided by 10)(tau_w)
%   P(7) - gain parameter b = b0 V0  (divided by 10)      (b)
%   P(8) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   
%
%   P(8 + 1:m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Zheng Mayhew, Time Invariant visco-elastic windkessel model relating
% blood flow and blood volume, NeuroImage 47 (2009) 1371.
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
x(2:end) = exp(x(2:end)); 

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(3)^(1/P(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - P(5))^(1/x(2)))/P(5);

% implement differential state equations
%--------------------------------------------------------------------------
f(1)     = P(9:end)'*u(:) - P(1)*x(1) - P(2)*(x(2) - 1);
f(2)     = x(1)/x(2); 
f(3)     = (x(2) - fv/x(5))/(P(3)*x(3)); 
f(4)     = (ff*x(2) - fv*x(4)/x(3))/(P(3)*x(4));
f(5)     = (-x(5)+exp(-10*P(7)*f(3)))/(10*P(6)*x(5));
f        = f(:);

% adjust motion for DEM (that uses time-bins as units of time)
%--------------------------------------------------------------------------
try, global dt, f  = f*dt; end

return