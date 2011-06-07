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
% /Y(t0)\   /Xsens*[ext] Xsens*[ext]\  /omega\   /epsilon_channel-noise\
% |  0  | = |      1           0    |* \beta0/ + |       -omega        |
% \  0  /   \     0           1     /            \       -beta0        /
%
% Y(t0)                 = instant t0 dans le fichier .nirs
% Xsens                 = matrice sensitivite en un unique point temporel
% [ext]                 = matrice des coefficients d extinction
% epsilon_channel-noise = bruit dans les canaux
% omega (notre beta !)  = effet du paradigme

%%%% lire le BOLD du contraste correspondant et prevoir d'enregistrer sous
%%%% differents noms em fonction du contrates puisqu apres il faudra
%%%% continuer avec le GLM et le meme contraste


for Idx=1:size(job.NIRSmat,1)
    
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        
        dir_in = job.dir_in{:};
        daate = strrep(datestr(now),':','-');
        mkdir(dir_in,daate);
        
        % gets current simulation cs
        sep = strfind(dir_in,'\');
        csn = dir_in(sep(end-1)+3:sep(end)-1);
        itest=1;
        while itest<length(NIRS.Cs.n) && (isempty(strfind(csn,NIRS.Cs.n{itest})) || length(csn)~=length(NIRS.Cs.n{itest}))
            itest =itest+1;
        end
        i_cs =itest;
        cs = NIRS.Cs.mcs{i_cs};
        
        %         fid=fopen(cs.b8i,'rb');% 8bits isotropic voxel volume has voxels in the same order as sensitivity matrix
        %         m = fread(fid);% currently not used
        %         fclose(fid);
        % b8i and segR have the same orientations as there is no
        % permutation (C and Matlab conventions in nirs_run_configMC)
        VsegR = spm_vol(cs.segR);
        jobR.image_in ={cs.segR};
        jobR.out_dir = dir_in;
        jobR.out_dim = [1 1 1];
        jobR.out_dt = 'same';
        jobR.out_vxsize = job.sens_vxsize;
        jobR.out_autonaming = 0;
        cs.segRR =  nirs_resize(jobR);
        
        % segRR has the same size as the sensitivity matrix
        VsegRR = spm_vol(cs.segRR);
        YsegRR = spm_read_vols(VsegRR);
        YsegRR = reshape(YsegRR,[1 prod(VsegRR.dim)]);

%         m_c1 = zeros(size(YsegRR));
%         m_c1(YsegRR==1)=1;
%         m_c5 = zeros(size(YsegRR));
%         m_c5(YsegRR==5)=1;
        %positions for Qp
        Y_c1 = find(YsegRR==1);% mask for GM
        Y_c5 = find(YsegRR==5);% mask for skin
        clear YsegRR;
        
        % Pairs....
        C_cs = cs.C; %ancien Cmc
        NC_cs = length(C_cs); %Total number of measurements
        
        %%% X %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % p330
        % on prend l'identite pour \Omegachapeau en premiere approximation
        %%% pour l'instant pas de SVD
        load(fullfile(dir_in,'sens.mat'));

        % On veut reconstruire efficacement. Comme on a une resolution
        % proche du cm, on sous echantillonne la matrice de sensitivite...
        for Ci =1:NC_cs
            Xmci = reshape(sens(Ci,:),VsegR.dim);
            Vmc = struct('fname',fullfile(dir_in,'Xmci.nii'),...
                'dim',  VsegR.dim,...
                'dt',   VsegR.dt,...
                'pinfo',VsegR.pinfo,...
                'mat',  VsegR.mat);
            Vmc = spm_create_vol(Vmc);
            spm_write_vol(Vmc, Xmci);
            
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
        % Xmc : matrice de sensitivite retaille avec des plus gros voxels !
        Xmc = Xmc_cm;
        clear Xmc_cm Vmc Ymc sens
        Nvx = size(Xmc,2);
        NC2mi = NC_cs/2;
        
        Xwl{1} = sparse(Xmc(1:NC2mi,:));
        Xwl{2} = sparse(Xmc(NC2mi+1:end,:));
        clear Xmc
        
        beta_prior = zeros(2*Nvx,1);
        
        ext1 = GetExtinctions(NIRS.Cf.dev.wl(1,1));
        ext2 = GetExtinctions(NIRS.Cf.dev.wl(1,2));
%         %%% Matrice epsilons (coefficients d'extcintion) %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Coefficients d'extinction selon Ted Huppert's HOMer (GetExtinctions.m)
%         % à, respectivement, 690 (ligne1, HbO HbR) et 830 (ligne2, HbO HbR)
%         
%         % [ Lamda1-HbO  Lambda1-HbR  (690=lambda1 avec le CW6)
%         %   Lambda2-HbO Lambda2-HbR] (830=lambda2 avec le CW6)
%         % size : nLambda x nHb (2x2)
%         
%         % source : http://omlc.ogi.edu/spectra/hemoglobin/summary.html (au 22
%         % juillet 2008)
%         if lambda(1)==830 && lambda(1)==690
%             coeff_ext = [ 974  693.04 ;   % POUR DONNÉES PRISES AVEC LE CW5
%                 276  2051.96  ] .* 2.303 ./1e6; % en cm^-1 / (umol/L).
%         elseif lambda(1)==690 && lambda(1)==830
%             coeff_ext = [ 276  2051.96 ;   % POUR DONNÉES PRISES AVEC LE CW6
%                 974  693.04  ] .* 2.303 ./1e6; % en cm^-1 / (umol/L).
%             % Le facteur 2.303 permet d'obtenir des coefficients d'absorption
%             % mu_a lorsqu'on mutliplie par la concentration en umol/L (toujours
%             % selon le site des données compilées par Scott Prahl!)
      
        for ifnirs=1:size(NIRS.Dt.fir.pp,2)         
            fnirs = load(NIRS.Dt.fir.pp.p{1,ifnirs},'-mat');
            
            %%% X %%% une matrice par session
            Xsens = sparse([Xwl{1} zeros(size(Xwl{1}));zeros(size(Xwl{1})) Xwl{2}]);
            E11 = sparse(ext1(1,1)*eye(Nvx));
            E12 = sparse(ext1(1,2)*eye(Nvx));
            E21 = sparse(ext2(1,1)*eye(Nvx));
            E22 = sparse(ext2(1,2)*eye(Nvx));
            Egrande = sparse([E11 E12 ; E21 E22]);
            X = Xsens*Egrande;
            
            %%% covariances %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Qn : Covariance du bruit de mesure
            % huppert_direct_2008 : (section 3.6 p 12) The initial seed of R was calculated from linear regression of the data with the temporal basis.
            for idx=1:size(NIRS.Cf.dev.wl,2) % over the wavelengths
                lst=(1:NC_cs/2)+NC_cs/2*(idx-1);
                Qn{idx}=sparse(lst,lst,ones(size(lst)),NC_cs,NC_cs);
            end
            %  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
            Qp{1}=sparse(Y_c5,Y_c5,ones(length(Y_c5),1),2*Nvx,2*Nvx);         %Skin layer- HbO % Prior knowledge of location of ROI (eq 22 huppert_hierarchical_2010)
            Qp{2}=sparse(Y_c1,Y_c1,ones(length(Y_c1),1),2*Nvx,2*Nvx);         %Brain layer- HbO
            Qp{3}=sparse(Nvx+Y_c5,Nvx+Y_c5,ones(length(Y_c5),1),2*Nvx,2*Nvx); %Skin layer- HbR
            Qp{4}=sparse(Nvx+Y_c1,Nvx+Y_c1,ones(length(Y_c1),1),2*Nvx,2*Nvx); %Brain layer- HbR
            Qp{5}=sparse(1:Nvx,1:Nvx,-ones(Nvx,1),2*Nvx,2*Nvx);               % Quantifie la covariance entre HbO et HbR dans chacun des voxels
%             if ~isempty(job.subj(1,is).boldmask{:})                           % BOLD
%                 %Qp{6}=sparse(1);
%             elseif ~isfield(NIRS.Cm,'bold')
%                 %'NIRS.Cm.bold'
%                 %Qp{6}=sparse();
%             end

            for itp=1:length(job.temp_pts)
                disp(['current : ' int2str(job.temp_pts(itp))])
                %%% Y %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Y_t0 = fnirs.d(job.temp_pts(itp),C_cs)';
              
                switch job.ReML_method
                    case 0
                        disp('code Huppert');
                        meth = 'HUP';
                        % on doit envoyer X et Y pas les bar !!!!!
                        
                        
                        [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,X,beta_prior,Qn,Qp);
%                         [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,X*W',beta_prior,Qn,Qp);
                        %Convert to the image domain and display
%                         beta = W'*beta_W;
                        
                    case 1
                        disp('code spm_reml');
                        meth = 'SPM';
                        
                        %%% Y %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Y_t0 = fnirs.d(job.temp_pts(itp),C_cs)';
         
                        %Set up the extended covariance model by concatinating the measurement
                        %and parameter noise terms
                        Q=cell(length(Qn)+length(Qp),1);
                        for idx=1:length(Qn)
                            Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
                        end
                        for idx2=1:length(Qp)
                            Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
                        end
                        
                        % sample covariance matrix Y*Y'
                        YY = (Y_t0-mean(Y_t0))*(Y_t0-mean(Y_t0))';
                        [C,h,Ph,F,Fa,Fc]=spm_reml(YY,X,Q);
                        iC     = spm_inv(C);
                        iCX    = iC*X;
                        Cq = spm_inv(X'*iCX);
                        beta = Cq*X'*iC*Y_t0;
                end
                
                betaR_HbO = reshape(beta(1:Nvx,1),[VsegRR.dim]);
                betaR_HbR = reshape(beta(Nvx+1:2*Nvx,1),[VsegRR.dim]);

                %write nifti for DHbO DHbR
                V_O = struct('fname',fullfile(dir_in,daate,['D[HbO]_' meth '_t' int2str(job.temp_pts(itp)) '.nii']),...
                    'dim',  VsegRR.dim,...
                    'dt',   [16,0],...
                    'pinfo',VsegRR.pinfo,...
                    'mat',  VsegRR.mat);
                V_O = spm_create_vol(V_O);
                spm_write_vol(V_O, betaR_HbO);
                
                V_R = struct('fname',fullfile(dir_in,daate,['D[HbR]_' meth '_t' int2str(job.temp_pts(itp)) '.nii']),...
                    'dim',  VsegRR.dim,...
                    'dt',   [16,0],...
                    'pinfo',VsegRR.pinfo,...
                    'mat',  VsegRR.mat);
                V_R = spm_create_vol(V_R);
                spm_write_vol(V_R, betaR_HbR);
            end
        end
        save(job.NIRSmat{Idx,1},'NIRS');
    catch exception
        disp(exception.identifier);
        rmdir(daate);
        disp(['Could not run MonteCarlo reconstruction for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;

% % % % % % %                 case 2 % ON FAIT POUR CHAQUE LONGUEUR D'ONDE SEPAREMENT
% % % % % % %                     tic
% % % % % % %                     %%% Y %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %                     fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
% % % % % % %                     Y_t0 = fnirs.d(itp,C_cs)';
% % % % % % %
% % % % % % %                     %%% X %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %                     % p330
% % % % % % %                     % on prend l'identite pour \Omegachapeau en premiere approximation
% % % % % % %                     Kmat = load(fullfile(dir_in,'sens.mat'));
% % % % % % %                     Xmc = Kmat.sens;
% % % % % % %                     clear Kmat
% % % % % % %                     % On veut reconstruire efficacement. Comme on a une resolution
% % % % % % %                     % proche du cm, on sous echantillonne la matrice de sensitivite...
% % % % % % %                     for Ci =1:NC_cs
% % % % % % %                         Xmci = reshape(Xmc(Ci,:),VsegR.dim);
% % % % % % %                         Vmc = struct('fname',fullfile(dir_in,'Xmci.nii'),...
% % % % % % %                             'dim',  VsegR.dim,...
% % % % % % %                             'dt',   VsegR.dt,...
% % % % % % %                             'pinfo',VsegR.pinfo,...
% % % % % % %                             'mat',  VsegR.mat);
% % % % % % %                         Vmc = spm_create_vol(Vmc);
% % % % % % %                         spm_write_vol(Vmc, Xmci);
% % % % % % %
% % % % % % %                         jobR.out_autonaming = 1;
% % % % % % %                         jobR.out_prefix = 'R';
% % % % % % %                         out =  nirs_resize(jobR);
% % % % % % %                         clear Xmci Vmc
% % % % % % %                         Vmc = spm_vol(out);
% % % % % % %                         Ymc = spm_read_vols(Vmc);
% % % % % % %                         Xmc_cm(Ci,:) = reshape(Ymc,[1 prod(Vmc.dim)]);
% % % % % % %                     end
% % % % % % %                     Xmc = Xmc_cm;
% % % % % % %                     clear Xmc_cm Vmc Ymc
% % % % % % %
% % % % % % %                     NC2mi = NC_cs/2;
% % % % % % %                     %         X = sparse([[Xmc(1:NCdemi,:)*415.5 Xmc(1:NCdemi,:)*2141.8];[Xmc(NCdemi+1:end,:)*1008.0 Xmc(NCdemi+1:end,:)*778.0]]);
% % % % % % %                     Xwlt{1} = sparse([Xmc(1:NC2mi,:)*415.5 Xmc(1:NC2mi,:)*2141.8]);
% % % % % % %                     Xwlt{2} = sparse([Xmc(NC2mi+1:end,:)*1008.0 Xmc(NC2mi+1:end,:)*778.0]);
% % % % % % %
% % % % % % %                     for iwl=1:2
% % % % % % %                         Xwl = Xwlt{1,iwl};
% % % % % % %
% % % % % % %                         %% Qn : Covariance pour la longueur d'onde
% % % % % % %                         lst=(1:NC2mi);
% % % % % % %                         Qn=sparse(lst,lst,ones(size(lst)),NC2mi,NC2mi);
% % % % % % %                         %%  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
% % % % % % %                         Qp{1}=sparse(Yb8i_c5,Yb8i_c5,ones(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbO
% % % % % % %                         Qp{2}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbO
% % % % % % %                         Qp{3}=sparse(Yb8i_c5,Yb8i_c5,ones(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbR
% % % % % % %                         Qp{4}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbR
% % % % % % %
% % % % % % %                         %% On prepare les reconstructions :
% % % % % % %                         switch job.ReML_method
% % % % % % %                             case 0
% % % % % % %                                 disp('code Huppert');
% % % % % % %                                 meth = 'HUP';
% % % % % % %                                 [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,X*W',beta_prior,Qn,Qp);
% % % % % % %                                 %Convert to the image domain and display
% % % % % % %                                 beta = W'*beta_W;
% % % % % % %
% % % % % % %                             case 1
% % % % % % %                                 disp('code spm_reml');
% % % % % % %                                 meth = 'SPM';
% % % % % % %                                 %Set up the extended covariance model by concatinating the measurement
% % % % % % %                                 %and parameter noise terms
% % % % % % %                                 Q=cell(length(Qn)+length(Qp),1);
% % % % % % %                                 for idx=1:length(Qn)
% % % % % % %                                     Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
% % % % % % %                                 end
% % % % % % %                                 for idx2=1:length(Qp)
% % % % % % %                                     Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
% % % % % % %                                 end
% % % % % % %                                 % sample covariance matrix Y*Y'
% % % % % % %                                 YY = (Y_t0-mean(Y_t0))*(Y_t0-mean(Y_t0))';
% % % % % % %                                 [C,h,Ph,F,Fa,Fc]=spm_reml(YY,X,Q);
% % % % % % %                                 iC     = spm_inv(C);
% % % % % % %                                 iCX    = iC*X;
% % % % % % %                                 Cq = spm_inv(X'*iCX);
% % % % % % %                                 beta = Cq*X'*iC*Y_t0;
% % % % % % %
% % % % % % %                                 %                             case 2
% % % % % % %                                 %                                 disp('peudo inverse');
% % % % % % %                                 %                                 meth = 'PInv';
% % % % % % %                                 %                                 % Clement's version
% % % % % % %                                 %
% % % % % % %                                 %                                 % % % % %         Ybar = sparse([Y_t0;zeros(3*size(X,2),1)]);
% % % % % % %                                 %                                 % M_c1 Mask for cortex
% % % % % % %                                 %                                 M_c1_wl = [m_c1+m_c5 m_c1+m_c5];
% % % % % % %                                 %                                 M_c1_wl = [m_c1 m_c1];
% % % % % % %                                 %                                 M_c1 = sparse(diag(M_c1_wl));
% % % % % % %                                 %                                 % % % % %         Xbar = sparse([log(X) log(X)*M_c1 log(X)*M_c1; sparse(1:3*size(X,2),1:3*size(X,2),ones(3*size(X,2),1),3*size(X,2),3*size(X,2))]);
% % % % % % %                                 %                                 Ybar = Y_t0((iwl-1)*NC2mi+(1:NC2mi));
% % % % % % %                                 %
% % % % % % %                                 %                                 % begin SVD :
% % % % % % %                                 %                                 % Chapitre 26 : p330
% % % % % % %                                 %                                 % observation noise : hatOmega pour nous Qn
% % % % % % %                                 %                                 Qinv =blkdiag(Qn); % hatsigma = Q
% % % % % % %                                 %                                 Q = inv(Qinv);
% % % % % % %                                 %                                 Xbar = sparse(Xwl*M_c1);
% % % % % % %                                 %
% % % % % % %                                 %                                 %         % custom SVD
% % % % % % %                                 %                                 %                     Xbar = sparse(Q.^(1/2)*full(Xwl*M_c1));
% % % % % % %                                 %                                 %         XbarN = Xbar'*Xbar;
% % % % % % %                                 %                                 %         [v S v] = svd(XbarN,0);
% % % % % % %                                 %                                 %         S       = sparse(S);
% % % % % % %                                 %                                 %         s       = diag(S);
% % % % % % %                                 %                                 %         j       = find(s*length(s)/sum(s) >= U & s >= T);
% % % % % % %                                 %                                 %         v       = v(:,j);
% % % % % % %                                 %                                 %         u       = spm_en(Xbar*v);
% % % % % % %                                 %                                 %         S       = sqrt(S(j,j));
% % % % % % %                                 %                                 %         % replace in full matrices
% % % % % % %                                 %                                 %         %---------------------------------------------------------------------------
% % % % % % %                                 %                                 %         j      = length(j);
% % % % % % %                                 %                                 %         U      = sparse(M,j);
% % % % % % %                                 %                                 %         V      = sparse(N,j);
% % % % % % %                                 %                                 %         if j
% % % % % % %                                 %                                 %             U(p,:) = u;
% % % % % % %                                 %                                 %             V(q,:) = v;
% % % % % % %                                 %                                 %         end
% % % % % % %                                 %                                 %
% % % % % % %                                 %                                 %         %%% voir les valeurs propres
% % % % % % %                                 %                                 %         %         i=1;
% % % % % % %                                 %                                 %         %         while S(i,i)>0 && i<=40
% % % % % % %                                 %                                 %         %             disp([int2str(i) ' ieme valeur propre de S : ' num2str(S(i,i))])
% % % % % % %                                 %                                 %         %             i = i+1;
% % % % % % %                                 %                                 %         %         end
% % % % % % %                                 %                                 %         %%%
% % % % % % %                                 %                                 %         Vbar = S*V';
% % % % % % %                                 %                                 %         Xbar = Vbar;
% % % % % % %                                 %                                 %         % end SVD
% % % % % % %                                 %
% % % % % % %                                 %                                 % TIKHONOV AMELIORE
% % % % % % %                                 %                                 % Beta contient omega_space omega et beta_prior : pour Tikhonov pas
% % % % % % %                                 %                                 % besoin de le definir puisque c'est nul...
% % % % % % %                                 %                                 % Betabar = sparse(,,,beta_prior);
% % % % % % %                                 %
% % % % % % %                                 %                                 % ebar DE MEME, par contre on definit la matrice des covariances
% % % % % % %                                 %                                 % des erreurs :
% % % % % % %                                 %                                 %         coef = 0.1;%%%%% moyen de calculer ca sur les images ??????????????
% % % % % % %                                 %                                 % of course : idee : en pratique surtout au niveau des interfaces,
% % % % % % %                                 %                                 % peut etre sortir l'info de ci_ fournie par newsegment puisque
% % % % % % %                                 %                                 % c'est des cartes de probabilite....
% % % % % % %                                 %                                 %         Qs=sparse(1:2*Nvx,1:2*Nvx,coef*ones(2*Nvx,1),2*Nvx,2*Nvx); % omega_space
% % % % % % %                                 %                                 %Set up the extended covariance model by concatinating the measurement
% % % % % % %                                 %                                 %and parameter noise terms and spatial prior
% % % % % % %                                 %                                 %         Q =blkdiag(Qn{1}+Qn{2},Qs,Qp{1}+Qp{2}+Qp{3}+Qp{4},Qp{1}+Qp{2}+Qp{3}+Qp{4});
% % % % % % %                                 %
% % % % % % %                                 %
% % % % % % %                                 %                                 % on applique ensuite la formule de l'inversion :
% % % % % % %                                 %                                 %%%cas de la svd
% % % % % % %                                 %                                 %         Beta_estimate = (Xbar'*Xbar) \ (Xbar'*U'*Ybar);
% % % % % % %                                 %
% % % % % % %                                 %                                 %%% pour sauver de l'espace memoire: from Philippe
% % % % % % %                                 %                                 %%%%%% je pense qu'on n'a pas le droit d'utiliser le slash a la
% % % % % % %                                 %                                 %%%%%% place de pinv...
% % % % % % %                                 %                                 %%% si l'on voulais prendre en compte les covariances, se rapporter
% % % % % % %                                 %                                 %%% au fichier TeX :
% % % % % % %                                 %                                 %                     clear M_c1 M_c1_wl W X Xmc Yb8i_c1 Yb8i_c5 beta_prior m_c1
% % % % % % %                                 %                                 alpha =1;
% % % % % % %                                 %                                 XX =Xbar'*Xbar;
% % % % % % %                                 %                                 YY = (Xbar'*Ybar);
% % % % % % %                                 %                                 clear Xbar;
% % % % % % %                                 %                                 XXLI = sparse(XX + eye(size(XX,2)));
% % % % % % %                                 %                                 %                     METHODE 1
% % % % % % %                                 %                                 % % % % % % % % % % %                     PAS SUUUUUUUUUUUR : je pense aue
% % % % % % %                                 %                                 % c4est fqux / on va essayer avec une plus grosse taille de voxels
% % % % % % %                                 %                                 % % % % % % % % % % %                     pinvXX = XXLI \ eye(size(XXLI,1));
% % % % % % %                                 %                                 % % % % % % % % % % %                         %                 pinvXX = pinv(XX + eye(size(XX,2)));
% % % % % % %                                 %                                 % % % % % % % % % % %
% % % % % % %                                 %                                 % % % % % % % % % % %                     %save the inverse
% % % % % % %                                 %                                 % % % % % % % % % % %                     save([dir_in '\pinvXX_' int2str(iwl) '.mat'],'pinvXX','-v7.3');
% % % % % % %                                 %                                 % % % % % % % % % % %                     save([dir_in '\YY_' int2str(iwl) '.mat'],'YY','-v7.3');
% % % % % % %                                 %                                 % % % % % % % % % % %                     %free up some memory
% % % % % % %                                 %                                 % % % % % % % % % % %                     clear Xbar;
% % % % % % %                                 %                                 % % % % % % % % % % %
% % % % % % %                                 %                                 % % % % % % % % % % %                     Beta_estimate = pinvXX*YY;
% % % % % % %                                 %                                 %                     METHODE 2
% % % % % % %                                 %                                 %free up some memory
% % % % % % %                                 %                                 clear Xbar;
% % % % % % %                                 %
% % % % % % %                                 %                                 Beta_estimate = XXLI \ YY;
% % % % % % %                                 %                                 % % % % % % % % % % % % % % % % % RETOUR AU CODE APRES INVERSION % % % % % % % % % % % % %
% % % % % % %                                 %                                 % TIKHONOV AMELIORE
% % % % % % %                                 %                                 %         beta = Beta_estimate((size(Y_t0,1)+2*size(Xbar,2)+1):size(Beta_estimate,1),1);
% % % % % % %                                 %                                 beta = Beta_estimate;
% % % % % % %                         end
% % % % % % %                         beta_4d = reshape(full(beta),[VsegRR.dim 2]);
% % % % % % %
% % % % % % %                         beta_HbO = beta_4d(:,:,:,1);
% % % % % % %                         beta_HbR = beta_4d(:,:,:,2);
% % % % % % %                     end
% % % % % % %             end