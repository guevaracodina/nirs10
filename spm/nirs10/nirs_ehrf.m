function [hrf,p] = nirs_ehrf(RT,P)
% Returns a hemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,[p])
% RT   - scan repeat time
% p    - parameters of the response function (two gamma functions)
%
%                                                     defaults
%                                                    (seconds)
%   p(1) - delay of response (relative to onset)         8
%   p(2) - dispersion of response                        4
%   p(3) - onset (seconds)                               0
%   p(4) - length of kernel (seconds)                   32
%
% hrf  - hemodynamic response function
% p    - parameters of the response function
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hrf.m 3716 2010-02-08 13:58:09Z karl $


% global parameter
%--------------------------------------------------------------------------
try
    fMRI_T = spm_get_defaults('stats.fmri.t');
catch
    fMRI_T = 16;
end

% default parameters
%--------------------------------------------------------------------------
p   = [8 4 0 32]; %gives a peak between 4 and 5 s, still substantial after 10 to 15 s
if nargin > 1
    p(1:length(P)) = P;
end

% modelled hemodynamic response function 
%--------------------------------------------------------------------------
dt  = RT/fMRI_T;
u   = [0:(p(4)/dt)] - p(3)/dt;
hrf = spm_Gpdf(u,p(1)/p(2),dt/p(2));
hrf = hrf([0:(p(4)/RT)]*fMRI_T + 1);
hrf = hrf'/sum(hrf);
