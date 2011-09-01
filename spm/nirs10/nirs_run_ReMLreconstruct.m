function out = nirs_run_ReMLreconstruct(job)
% Clément Bonnéry May 2011

% Le systeme qu'on resout est le suivant :
% on est en un UNIQUE POINT TEMPOREL donc
%              ---------------------
%
% / Y(t0) = Xsens*[ext]*[DHbO(vx);DHbR(vx)] + epsilon_channel-noise
% \ [DHbO(vx);DHbR(vx)] = 0 + omega
%
% ce que l on peut ecrire sous la forme
% /Y(t0)\   /Xsens*[ext] Xsens*[ext]\  /omega\   /epsilon_channelnoise\
% |  0  | = |      1           0    |* \beta0/ + |       -omega        |
% \  0  /   \     0           1     /            \       -beta0        /
%
% Y(t0)                 = instant t0 dans le fichier .nirs
% Xsens                 = matrice sensitivite en un unique point temporel
% [ext]                 = matrice des coefficients d extinction
% epsilon_channelnoise = bruit dans les canaux
% omega (notre beta !)  = effet du paradigme

%%%% lire le BOLD du contraste correspondant et prevoir d'enregistrer sous
%%%% differents noms em fonction du contrates puisqu apres il faudra
%%%% continuer avec le GLM et le meme contraste


for Idx=1:size(job.NIRSmat,1)
    
    try
        clear NIRS
        load(job.NIRSmat{Idx,1});
        %selection of working directory
        if ~isempty(job.dir_in{1})
            f = job.dir_in;
            cs_dir = fileparts(f{1,:});
            [dummy cs_ldir] = fileparts(cs_dir); %name of simulation
            ics = 1;
            while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                ics =ics+1;
            end
            %current simulation
            cs = NIRS.Cs.mcs{ics};
        else
            %take last simulation
            cs = NIRS.Cs.mcs{end};
            cs_dir = cs.dir;
        end
        
        ctm = {};
        
        % b8i and segR have the same orientations since there is no
        % permutation (C and Matlab conventions in nirs_run_configMC)
        VsegR = spm_vol(cs.segR);
        jobR.image_in ={cs.segR};
        jobR.out_dir = cs_dir;
        jobR.out_dim = [1 1 1];
        jobR.out_dt = 'same';
        jobR.out_vxsize = job.sens_vxsize;
        jobR.out_autonaming = 0;
        temp.segRR =  nirs_resize(jobR);
        
        % segRR has the same size as the sensitivity matrix
        VsegRR = spm_vol(temp.segRR);
        YsegRR = spm_read_vols(VsegRR);
        YsegRR = reshape(YsegRR,[1 prod(VsegRR.dim)]);
        
        % positions for Qp
        Y_c1 = find(YsegRR==1);% mask for GM
        Y_c5 = find(YsegRR==5);% mask for skin
        clear YsegRR;
        
        % Pairs....
        %C_cs = [4 32];
        C_cs = cs.C; %ancien Cmc
        NC_cs = length(C_cs); %Total number of measurements
        
        %%% X %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % p330
        % on prend l'identite pour \Omegachapeau en premiere approximation
        %%% pour l'instant pas de SVD
        load(fullfile(cs_dir,'sens.mat'));
        Xmc = sens;
        clear sens
        % On veut reconstruire efficacement. Comme on a une resolution
        % proche du cm, on sous echantillonne la matrice de sensitivite...
        for Ci =1:NC_cs
            %inefficient, why write .nii
            Xmci = reshape(Xmc(Ci,:),VsegR.dim);
            Vmc = nirs_create_vol(fullfile(cs_dir,'Xmci.nii'),...
                VsegR.dim, VsegR.dt, VsegR.pinfo, VsegR.mat, Xmci);
            
            
            %this resizing has the same parameters has the one of segR to
            %segRR (lines 47 to 52)
            jobR.image_in ={Vmc.fname};
            jobR.out_autonaming = 1;
            jobR.out_prefix = 'R';
            out =  nirs_resize(jobR);
            clear Xmci Vmc
            Vmc = spm_vol(out);
            Ymc = spm_read_vols(Vmc);
            Xmc_cm(Ci,:) = reshape(Ymc,[1 prod(Vmc.dim)]);
        end
        % Xmc : matrice de sensitivite retaillee avec des plus gros voxels !
        Xmc = Xmc_cm;
        clear Xmc_cm Vmc Ymc sens
        Nvx = size(Xmc,2);
        NC2mi = NC_cs/2;
        
        Xwl{1} = sparse(Xmc(1:NC2mi,:));
        Xwl{2} = sparse(Xmc(NC2mi+1:end,:));
        sXmc = sparse(Xmc);%%%%%%%% SHBTB
        clear Xmc
        
        beta_prior = zeros(2*Nvx,1);
        
        ext1 = GetExtinctions(NIRS.Cf.dev.wl(1));
        ext2 = GetExtinctions(NIRS.Cf.dev.wl(2));
        Xsens = blkdiag(Xwl{1}, Xwl{2});
        
        E0 = [ext1(1) ext1(2); ext2(1) ext2(2)];
        sX = Xsens * kron(E0,sparse(eye(Nvx)));
        clear Xsens
        switch job.ReML_method
            case 0
            case 1
                tic%%%%%%%% SHBTB %This requires a lot of memory, takes about 1 minute
                Xbar = sparse([sXmc sX ; speye(Nvx) sparse(zeros(Nvx,2*Nvx)) ; sparse(zeros(2*Nvx,Nvx)) speye(2*Nvx)]);%%%%%%%% SHBTB
                toc%%%%%%%% SHBTB
        end
        for ifnirs=1:size(NIRS.Dt.fir.pp(end),2)
            fnirs = load(NIRS.Dt.fir.pp(end).p{1,ifnirs},'-mat');
            ctm.Y = NIRS.Dt.fir.pp(end).p{1,ifnirs};
            
            %%% covariances %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Qn : Covariance du bruit de mesure
            % huppert_direct_2008 : (section 3.6 p 12) The initial seed of R was calculated from linear regression of the data with the temporal basis.
            for idx=1:size(NIRS.Cf.dev.wl,2) % over the wavelengths
                lst=(1:NC_cs/2)+NC_cs/2*(idx-1);
                Qn{idx}=sparse(lst,lst,ones(size(lst)),NC_cs,NC_cs);
            end
            %  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
            Qp{1}=sparse(Y_c5,Y_c5,ones(length(Y_c5),1),2*Nvx,2*Nvx);         % Skin layer- HbO % Prior knowledge of location of ROI (eq 22 huppert_hierarchical_2010)
            Qp{2}=sparse(Y_c1,Y_c1,ones(length(Y_c1),1),2*Nvx,2*Nvx);         % Brain layer- HbO
            Qp{3}=sparse(Nvx+Y_c5,Nvx+Y_c5,ones(length(Y_c5),1),2*Nvx,2*Nvx); % Skin layer- HbR
            Qp{4}=sparse(Nvx+Y_c1,Nvx+Y_c1,ones(length(Y_c1),1),2*Nvx,2*Nvx); % Brain layer- HbR
            Qp{5}=sparse(Nvx+1:2*Nvx,1:Nvx,-ones(Nvx,1),2*Nvx,2*Nvx)+sparse(1:Nvx,1+Nvx:2*Nvx,-ones(Nvx,1),2*Nvx,2*Nvx); % Quantifie la covariance entre HbO et HbR dans chacun des voxels
            %             if ~isempty(job.subj(1,is).boldmask{:})                           % BOLD
            %                 %Qp{6}=sparse(1);
            %             elseif ~isfield(NIRS.Cm,'bold')
            %                 %'NIRS.Cm.bold'
            %                 %Qp{6}=sparse();
            %             end
            
            daate = strrep(datestr(now),':','-');
            tm_dir = ['Re_' daate '_a' ctm.alg '_' int2str(job.sens_vxsize) 'mm'];
            ctm.p = fullfile(cs_dir,tm_dir);
            if ~exist(ctm.p,'dir'),mkdir(ctm.p);end
            if ~isfield(NIRS,'Tm'), NIRS.Tm ={}; end
            if isfield(NIRS.Tm,'tmrs')
                itm = size(NIRS.Tm.tmrs,2)+1;
            else
                itm=1;
            end
            ctm.n = tm_dir;
            ctm.sens_vxsize = job.sens_vxsize;
            
            [dum,namRR,extRR] = fileparts(temp.segRR);
            ctm.segRR = fullfile(ctm.p,[namRR extRR]);
            copyfile(temp.segRR,ctm.segRR);
            delete(temp.segRR);
            NIRS.Tm.tmrs{itm} = ctm;
            NIRS.Tm.n{itm} = ctm.n;
                           
            for itp=1:length(job.temp_pts)
                tic
                disp(['current : ' int2str(job.temp_pts(itp))])
                %%% Y %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Y_t0 = fnirs.d(job.temp_pts(itp),C_cs)';
                
                %ctm.Y = 'fantom';
                %Yt0 = load('Yt0.nirs','-mat');
                %Y_t0 = Yt0.Yt0;
                
                switch job.ReML_method
                    case 0
                        disp('code Huppert');
                        ctm.alg = 'HUP';
                        % on doit envoyer X et Y pas les bar !!!!!
                        warning('off','MATLAB:nearlySingularMatrix');
                        warning('off','MATLAB:singularMatrix');
                        [lambda,beta,Stats]=nirs_run_DOT_REML(Y_t0,sX,beta_prior,Qn,Qp);
                        %                         [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,Xbar*W',beta_prior,Qn,Qp);
                        %Convert to the image domain and display
                        %                         beta = W'*beta_W;
                        %beta = beta_W;
                    case 1
                        
                        disp('code spm_reml');
                        ctm.alg = 'SPM';
                        Y_bar_t0 = sparse([Y_t0; zeros(3*Nvx,1)]);%%%%%%%% SHBTB
                        %%% Y %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %                         Y_t0 = fnirs.d(job.temp_pts(itp),C_cs)';
                        
                        %Set up the extended covariance model by concatinating the measurement
                        %and parameter noise terms
                        Q=cell(length(Qn)+length(Qp),1);
                        
                        %%%% ACHTUNG : les Q sont des matrices carres du
                        %%%% grandeur le nombre de canaux...
                        %%%% On leur applique la mem svd qu a Xbar !
                        %  Y = U*S*V'*Beta
                        %  Y = U*(S*V'*Beta) --> Y=U*S*Beta2;
                        %  cov(Beta2) = S*V'*Q*V*S';
                        [U,S,V]=svd(full(Xbar),'econ');
                        
                        for idx=1:length(Qp)
                            Qp2{idx}=S*V'*Qp{idx}*V*S';
                        end
                        beta_prior=S*V'*beta_prior;
                        %%%%
                        %%%%
                        for idx=1:length(Qn)
                            %                             Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
                            Q{idx}=Qn{idx};
                        end
                        for idx2=1:length(Qp)
                            %                             Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
                            Q{idx+idx2}=Qp2{idx2};
                        end
                        
                        % sample covariance matrix Y*Y'
                        YY = (Y_t0-mean(Y_t0))*(Y_t0-mean(Y_t0))';
                        %                         [C,h,Ph,F,Fa,Fc]=spm_reml_hijacked(YY,Xbar,Q);
                        [C,h,Ph,F,Fa,Fc]=nirs_spm_reml(YY,Xbar,Q);
                        iC     = spm_inv(C);
                        iCX    = iC*Xbar;
                        Cq = spm_inv(Xbar'*iCX);
                        beta = Cq*Xbar'*iC*Y_t0;
                end
                
                betaR_HbO = reshape(beta(1:Nvx,1),VsegRR.dim);
                betaR_HbR = reshape(beta(Nvx+1:2*Nvx,1),VsegRR.dim);                
                str0 = gen_num_str(itp,4);               
                V_O = nirs_create_vol(fullfile(ctm.p,['O_' str0 '.nii']),...
                    VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbO);
                V_R = nirs_create_vol(fullfile(ctm.p,['R_' str0 '.nii']),...
                    VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbR);
                toc
            end
        end
        save(job.NIRSmat{Idx,1},'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        if exist('daate','var')
            rmdir(ctm.p);
        end
        disp(['Could not run MonteCarlo reconstruction for subject' int2str(Idx)]);
    end
end
clear out %careful, variable out was used with another meaning previously!
out.NIRSmat = job.NIRSmat;