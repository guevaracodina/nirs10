function out = GetOpt_ppts(varargin)
% To use subject's specific values, please CODE your formula to create
% appropriate factors to scale absorption and diffusion coefficients.
%
% FORMAT opt_ppts = GetOpt_ppts(choice)
% choice           - NIRS matrix
%
% FORMAT opt_ppts = GetOpt_ppts(choice, p, isubj)
% choice - NIRS matrix
% p      - path to the file where subject's mua and musprime are stored
% isubj  - subject's number (same as in the name of the subject's folder) 
%
% If mean over values for young adults are not equal to litterature values 
% (Strangman, 2003 : Factors affecting the accuracy of near-infrared 
% spectroscopy concentration calculations for focal changes in oxygenation 
% parameters), these scalling factors could be (for a study with young and 
% old adults groups) : 
% mu(isubj) = mu(litt)*mu(isubj,modality)/mean(mu(young subjects,modality))
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry
% 2011

if length(varargin)>1
    choice ='subj';
    p = varargin{2};
    isubj = varargin{3};
else
    choice ='litt';
end

% mua mus g n
%690
opt_ppts(1,1,:) = [0.0178   12.5   0.9   1.4];
opt_ppts(1,2,:) = [0.0178   12.5   0.9   1.4];
opt_ppts(1,3,:) = [0.0004   0.06   0.9   1.4];
opt_ppts(1,4,:) = [0.0101   10.0   0.9   1.4];
opt_ppts(1,5,:) = [0.0159-0.006   8.00   0.9   1.4];
opt_ppts(1,6,:) = [0.000   0.0   0   1];% air inside sinuses
opt_ppts_perturb{1} =  opt_ppts(1,1,:);%job.MC_parameters.perturbationPpties_l1+
%830
opt_ppts(2,1,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,2,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,3,:) = [0.0026   0.06   0.9   1.4];
opt_ppts(2,4,:) = [0.0136   8.60   0.9   1.4];
opt_ppts(2,5,:) = [0.0191-0.006  6.60   0.9   1.4];
opt_ppts(2,6,:) = [0.000   0.0   0   1];% air inside sinuses
opt_ppts_perturb{2} =  opt_ppts(2,1,:);%job.MC_parameters.perturbationPpties_l1+

switch choice
    case 'litt'
        %default values are handed out
    case 'subj'
        load(p);
        
        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
        %%% TO BE CODED BY USER TO GENERATE APPROPRIATE SCALLING FACTORS
        %%% FOR ABSORPTION AND DIFFUSION COEFFICIENTS
        c_a1 = 1;
        c_s1 = 1;
        c_a2 = 1;
        c_s2 = 1;
        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

        op = opt_ppts;
        op(1,[1 2 4 5],1) = c_a1*op(1,[1 2 4 5],1);
        op(1,[1 2 4 5],2) = c_s1*op(1,[1 2 4 5],2);
        op(2,[1 2 4 5],1) = c_a2*op(2,[1 2 4 5],1);
        op(2,[1 2 4 5],2) = c_s2*op(2,[1 2 4 5],2);
        
        opt_ppts = op;
end
out{1} = opt_ppts;
out{2} = opt_ppts_perturb;
end