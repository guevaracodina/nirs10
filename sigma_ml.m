% sigma_ml() - Called by weighted_estimate.m
%
% Usage: = >> See weighted_estimate.m
%
% Inputs:
%   approx_residual
%   c
%   SaJ
%   SdJ
%   J
%   nbr_level
%
% Outputs:
%   var_out
%
% Note:

function var_out = sigma_ml(approx_residual,c,SaJ,SdJ,J,nbr_level)

temp = 0;
for index_j = 1:nbr_level
      
    temp = temp+sum(cell2mat(c{1,index_j}).^2)./SdJ(index_j);

end

var_out = 1/(2.^J)*(approx_residual.^2/SaJ+temp);