% SdJ() - fminsearch on CGamma. Called by weighted_estimate.m
%
% Usage: = >> See weighted_estimate.m
%
% Inputs:
%   cGamma
%   spectral_exponent
%   index_j
%
% Outputs:
%   var_detail
%
% Note:


function var_detail = SdJ(cGamma, spectral_exponent, index_j)

var_detail = cGamma.*(2.^index_j).^spectral_exponent.*(2-2.^spectral_exponent)./( (2.*pi).^spectral_exponent.*(1-spectral_exponent) );