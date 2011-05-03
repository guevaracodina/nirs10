% SaJ() - fminsearch on CGamma. Called by weighted_estimate.m
%
% Usage: = >> See weighted_estimate.m
%
% Inputs:
%   cGamma
%   spectral_exponent
%   J
%   cutoff
%
% Outputs:
%   var_approx
%
% Note:


function var_approx = SaJ(cGamma,spectral_exponent,J,cutoff)

var_approx = (cGamma.*((2.^(J+1)).^spectral_exponent))./...
    ((2.*pi).^spectral_exponent.*(1-spectral_exponent)) *...
    1/cutoff.^spectral_exponent*(cutoff.^spectral_exponent- 2^((J+1)-spectral_exponent*(J+1))*cutoff);