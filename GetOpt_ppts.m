function out = GetOpt_ppts(varargin)
% mua mus g n
%830
opt_ppts{2,1} = [0.0186   11.1   0.9   1.4];
opt_ppts{2,2} = [0.0186   11.1   0.9   1.4];
opt_ppts{2,3} = [0.0026   0.10   0.9   1.4];
opt_ppts{2,4} = [0.0136   8.60   0.9   1.4];
opt_ppts{2,5} = [0.0191   6.60   0.9   1.4];
opt_ppts{2,6} = [0.0001   0.01   1   1.4];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{2} =  opt_ppts{2,1};%job.MC_parameters.perturbationPpties_l1+

%690
opt_ppts{1,1} = [0.0178   12.5   0.9   1.4];
opt_ppts{1,2} = [0.0178   12.5   0.9   1.4];
opt_ppts{1,3} = [0.0004   0.10   0.9   1.4];
opt_ppts{1,4} = [0.0101   10.0   0.9   1.4];
opt_ppts{1,5} = [0.0159   8.00   0.9   1.4];
opt_ppts{1,6} = [0.0001   0.01   1   1.4];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{1} =  opt_ppts{1,1};%job.MC_parameters.perturbationPpties_l1+

out{1} = opt_ppts;
out{2} = opt_ppts_perturb;
end