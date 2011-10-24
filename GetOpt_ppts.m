function out = GetOpt_ppts(varargin)

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
opt_ppts(1,3,:) = [0.0004   0.10   0.9   1.4];
opt_ppts(1,4,:) = [0.0101   10.0   0.9   1.4];
opt_ppts(1,5,:) = [0.0159   8.00   0.9   1.4];
opt_ppts(1,6,:) = [0.000   0.0   0   1];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{1} =  opt_ppts(1,1,:);%job.MC_parameters.perturbationPpties_l1+
%830
opt_ppts(2,1,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,2,:) = [0.0186   11.1   0.9   1.4];
opt_ppts(2,3,:) = [0.0026   0.10   0.9   1.4];
opt_ppts(2,4,:) = [0.0136   8.60   0.9   1.4];
opt_ppts(2,5,:) = [0.0191   6.60   0.9   1.4];
opt_ppts(2,6,:) = [0.000   0.0   0   1];% air a l interieur de la tete : ex sinus
opt_ppts_perturb{2} =  opt_ppts(2,1,:);%job.MC_parameters.perturbationPpties_l1+

switch choice
    case 'litt'
        %default values are handed out
    case 'subj'
        load(p);
        load('D:\Users\Clément\DPF_testDuncan\P2S_age.mat');
        
        lim = 50;
        j = find(ages(:,2)<lim);
        v = find(ages(:,2)>lim);
        
        %         if ages(isubj,2)<lim%j
        %             c_a1 = mean(Opt_ppts(j,1,1))/mean(Opt_ppts(:,1,1));
        %             c_s1 = mean(Opt_ppts(j,1,2))/mean(Opt_ppts(:,1,2));
        %             c_a2 = mean(Opt_ppts(j,2,1))/mean(Opt_ppts(:,2,1));
        %             c_s2 = mean(Opt_ppts(j,2,2))/mean(Opt_ppts(:,2,2));
        %         else%v
        %             c_a1 = mean(Opt_ppts(v,1,1))/mean(Opt_ppts(:,1,1));
        %             c_s1 = mean(Opt_ppts(v,1,2))/mean(Opt_ppts(:,1,2));
        %             c_a2 = mean(Opt_ppts(v,2,1))/mean(Opt_ppts(:,2,1));
        %             c_s2 = mean(Opt_ppts(v,2,2))/mean(Opt_ppts(:,2,2));
        %         end
        
        % 
        load('D:\Users\Clément\DPF_article\analyses_dpf11\Mua_MichP2S44sujets.mat');
        load('D:\Users\Clément\DPF_article\analyses_dpf11\Mus_MichP2S44sujets.mat');
        % W weight en exponentielle negative
        W = [0.0321 0.0871 0.2369 0.6439];
        c_a1 = Mua(isubj,[1 2 3 4])*W'/mean(Mua(:,[1 2 3 4])*W');
        c_s1 = Musprime(isubj,[1 2 3 4])*W'/mean(Musprime(:,[1 2 3 4])*W');
        c_a2 = Mua(isubj,[13 14 15 16])*W'/mean(Mua(:,[13 14 15 16])*W');
        c_s2 = Musprime(isubj,[13 14 15 16])*W'/mean(Musprime(:,[13 14 15 16])*W');
        %         % calcul de la nouvelle valeurs des mu pour toutes les couches sauf
        %         % le CSF qui ne contient aucune vasculature.
        op = opt_ppts;
        op(1,[1 2 4 5],1) = c_a1*op(1,[1 2 4 5],1);
        op(1,[1 2 4 5],2) = c_s1*op(1,[1 2 4 5],2);
        op(2,[1 2 4 5],1) = c_a2*op(2,[1 2 4 5],1);
        op(2,[1 2 4 5],2) = c_s2*op(2,[1 2 4 5],2);
        opt_ppts = op;
        %%%%%%%%%%% DPF_article : MCsims
        % % % % % %         load('D:\Users\Clément\DPF_article\MCsims\scripts\muFactors');
        % % % % % %         op = opt_ppts;
        % % % % % %         op(1,[1 2 4 5],1) = mua_vF690*op(1,[1 2 4 5],1);
        % % % % % %         op(1,[1 2 4 5],2) = mus_vF690*op(1,[1 2 4 5],2);
        % % % % % %         op(2,[1 2 4 5],1) = mua_vF830*op(2,[1 2 4 5],1);
        % % % % % %         op(2,[1 2 4 5],2) = mus_vF830*op(2,[1 2 4 5],2);
        % % % % % %         opt_ppts = op;
        end
        out{1} = opt_ppts;
        out{2} = opt_ppts_perturb;
end