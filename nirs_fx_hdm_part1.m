function [f] = nirs_fx_hdm_part1(x,u,P,M)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1) - vascular signal                                    s
%   x(2) - rCBF                                           log(f)
% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(2 + 1:m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________

% LIOM - adapted from:

% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hdm.m 2495 2008-11-27 12:18:33Z karl $

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(2:end) = exp(x(2:end)); 


% implement differential state equations
%--------------------------------------------------------------------------
f(1)     = P(3:end)'*u(:) - P(1)*x(1) - P(2)*(x(2) - 1);
f(2)     = x(1)/x(2); 
f        = f(:);

% adjust motion for DEM (that uses time-bins as units of time)
%--------------------------------------------------------------------------
try, global dt, f  = f*dt; end

return