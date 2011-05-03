% loglike() - Called by weighted_estimate.m
%
% Usage: = >> See weighted_estimate.m
%
% Inputs:
%   cGamma
%   spectral_exponent
%   N
%   J
%   nbr_level
%   L
%   approx_coeff
%   c
%   cutoff
%
% Outputs:
%   LL
%
% Note:


function LL = loglike(cGamma, spectral_exponent, N, J,nbr_level, L,approx_coeff,c,cutoff)

temp =0;

%on recalcule sigma carré
var_approx = SaJ(cGamma,spectral_exponent,J,cutoff); % RaJ
for index_j = 1:nbr_level
    var_detail(index_j)=SdJ(cGamma,spectral_exponent, index_j); % les Rdj
end

Sigma2 = sigma_ml(approx_coeff, c, var_approx, var_detail,nbr_level, J); % equation 21



for index_j = 1:nbr_level
    temp = temp + N/2^index_j * log (SdJ(cGamma,spectral_exponent, index_j)); % 2eme partie de equation 20
end

var1 = -1/2*(  N*log(2*pi*Sigma2) +  log(SaJ(cGamma, spectral_exponent,J,cutoff))  + temp); %ajout de la premiere partie

temp = (approx_coeff).^2./SaJ(cGamma,spectral_exponent,J,cutoff); % 3eme partie de l'equation

temp2 = 0;
for index_j = 1:nbr_level
    temp2 = temp2+sum(cell2mat(c{1,index_j}).^2)./SdJ(cGamma,spectral_exponent,index_j);% 4eme partie
end

var2 = -1/(2*Sigma2)*(temp+temp2); % multiplication par -1/(2*Sigma2) de la 3eme et 4eme partie

LL = -(var1+var2); % finaly LL(gamma)
