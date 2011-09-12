function out = GetOpt_ppts(varargin)

if length(varargin)>1
    choice ='subj';
    p = varargin{2};
    isubj = varargin{3};
else
    choice ='litt';
end
% mua mus g n
%830
opt_ppts(2,1,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,2,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,3,:) = [0.0026   0.10   0.9   1.4];
opt_ppts(2,4,:) = [0.0136   8.60   0.9   1.4];
opt_ppts(2,5,:) = [0.0191   6.60   0.9   1.4];
opt_ppts(2,6,:) = [0.0001   0.01   1   1.4];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{2} =  opt_ppts(2,1,:);%job.MC_parameters.perturbationPpties_l1+

%690
opt_ppts(1,1,:) = [0.0178   12.5   0.9   1.4];
opt_ppts(1,2,:) = [0.0178   12.5   0.9   1.4];
opt_ppts(1,3,:) = [0.0004   0.10   0.9   1.4];
opt_ppts(1,4,:) = [0.0101   10.0   0.9   1.4];
opt_ppts(1,5,:) = [0.0159   8.00   0.9   1.4];
opt_ppts(1,6,:) = [0.0001   0.01   1   1.4];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{1} =  opt_ppts(1,1,:);%job.MC_parameters.perturbationPpties_l1+

switch choice
    case 'litt'
        %default values are handed out       
    case 'subj'
        load(p);     
        
%         c_a1 = Opt_ppts(isubj,1,1)/mean(Opt_ppts(:,1,1));
%         c_s1 = Opt_ppts(isubj,1,2)/mean(Opt_ppts(:,1,2));
%         c_a2 = Opt_ppts(isubj,2,1)/mean(Opt_ppts(:,2,1));
%         c_s2 = Opt_ppts(isubj,2,2)/mean(Opt_ppts(:,2,2));


        c_a1 = 2;
        c_s1 = 1;
        c_a2 = 2;
        c_s2 = 1;
        % calcul de la nouvelle valeurs des mu pour toutes les couches sauf
        % le CSF qui ne contient aucune vasculature.
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