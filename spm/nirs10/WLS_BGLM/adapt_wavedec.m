% adapt_wavedec() - Create L matrix used in Matlab's wavelet toolbox. Used
% by weighted_estimate.m
%
% Usage: = >>
%
% Inputs:
%   S: coeeficient du bruit ( residu)
%   J0: scale 
%
% Outputs:
%   L 1 by nb_Total_coeff+2
%
% Note:


function [L] = adapt_wavedec(S,J0)

n = length(S) ; %2^J coeff
J = ceil(log(n)/log(2)); %on retouve le J
if 2^J ~= n %control
    error('Dyadlength: n != 2^J')
end

L = zeros(1,J0+2); % Taille 1 by scale max+2
L(end)= n; % on met le nombre total de coeff dans la derniere case

temp = n;
for index=1:J0
    temp = temp/2;
    L(end-index) = temp; % dans l'avant derniere case on a le nombre de coeff de la scale 1 et ainsi de suite jusqua JO
end

L(1) = temp; % Dans la premiere case on met le nombre de coeff de scale JO ( correspond je pense au nombre de approx coeff