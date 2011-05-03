% atom_correlation() - Perform correlation between data and paradigm.
%
% Usage: = >> See weighted_estimate.m
%
% Inputs:
%   f
%   response: paradigm
%   L0
%   N
%
% Outputs:
%   matrix_corr
%
% Note:

function [matrix_corr] = atom_correlation(f, response,L0,N)

matrix_corr = zeros(1,N);
start_correlation = 1;
temp3 = zeros(1,N);

for index_xi = start_correlation:N
    
    temp3(index_xi)=1;
    wave_ind = IWT_PO(temp3,L0,f); %individual wavelet=psi(k,j) en fonction de la place du 1 dans la matrice temp
   
    %normalisation
    coeff = sum(response(1,:).^2)/length(response(1,:)); %moyenne quadratique
    resp_norm = response(1,:)/sqrt(coeff);
    coeff_2 = sum(wave_ind(1,:).^2)/length(wave_ind(1,:));
    wave_ind_norm = wave_ind(1,:)/sqrt(coeff_2);
    
    %perform correlation between paradigm and atom
    tmp=corrcoef(resp_norm',wave_ind');
    
    %remet le coeff a 0
    temp3(index_xi)=0;
    
    corr = tmp(1,2);
    matrix_corr(index_xi) = abs(corr);
    
end

