function out = nirs_run_ReMLreconstruct(job)

load(job.NIRSmat{:});
% on recupere cs
dir_in = job.dir_in{:};
sep = strfind(dir_in,'\');
csn = dir_in(sep(end-1)+1:sep(end)-1);

itest=1;
while itest<length(NIRS.Cs.n) && isempty(strfind(csn,NIRS.Cs.n{itest}))
    itest =itest+1;
end
i_cs =itest;
cs = NIRS.Cs.mcs{i_cs};

% volume 8bits et voxels isotropiques
% lecture du volume dans l'ordre de la matrice de sensitivite
fid=fopen(cs.b8i,'rb');
Yb8i = fread(fid);
fclose(fid);
Nvx = length(Yb8i);

% reperage des couches :
Yb8i_c1 = find(Yb8i==1);
Yb8i_c5 = find(Yb8i==5);
m_c1 = zeros(size(Yb8i));
m_c1(find(Yb8i==1))=1;
clear Yb8i;

% Le systeme qu'on resout est le suivant :
% on est en un UNIQUE POINT TEMPOREL donc
%              ---------------------
%
%/ Y = X W betaWvx + epsilon_channel-noise
%\ betaWvx = omega
%
% la matrice de sensitivite donne les delta mua
%/ Y = X W delta mua +epsilon_channel-noise
%\ delta mua = matrice de passage de mua a Hb. omega
%
% Y = un instant donne dans le fichier .nirs
% X = matrice sensitivite avec un point temporel
% W = matrice d'ondelettes spatiales // voire matrice de passage entre
% l'espace des voxels et un espace de taille plus petite
% betaWvx = concentrations dans les voxels
% epsilon_channel-noise = bruit dans les canaux
% omega = effet du paradigme

% Pairs....
Cmc = cs.C;
NC = length(Cmc); %Total number of measurements
Cwl=[];
Cwl = [Cwl NIRS.Cf.H.C.wl(Cmc)];
wl = unique(Cwl);

WLruns =2;% are we solving for both wl or not

switch WLruns
    case 1
        %% Y
        fnirs = load(NIRS.Dt.fir.pp(1,4).p{:},'-mat');
        t0 =3000;% on choisit au pif un point temporel
        Y_t0 = fnirs.d(t0,Cmc)';
        
        %% X
        % p330
        % on prend pour \Omegachapeau I en premiere approximation
        Kmat = load(fullfile(dir_in,'sens.mat'));
        Xmc = Kmat.sens;
        % maxi = max(Kmat.sens,[],2);
        % for i=1:size(Xmc,1)
        %     Xmc(i,:) = Xmc(i,:)/maxi(i);
        % end
        % 1.2956 et 0.1940 sont les rapports entre les coefficients d'absortivite
        % de HbO et HbR pour 690 et 830.
        NCdemi = size(Xmc,1)/2;
        X = sparse([[Xmc(1:NCdemi,:)*415.5 Xmc(1:NCdemi,:)*2141.8];[Xmc(NCdemi+1:end,:)*1008.0 Xmc(NCdemi+1:end,:)*778.0]]);
        % %Omega chapeau est le bruit d'observation. Ici on le prend egal a un
        % Omegachapo = eye(NC);
        % Kbarre = Omegachapo*K;
        %% omega
        beta_prior = zeros(2*Nvx,1);
        
        %% W
        % dans le cas d'une premiere reconstruction, on prend W =Id
        W = sparse((1:2*Nvx),(1:2*Nvx),ones(1,2*Nvx)',2*Nvx,2*Nvx);
        
        %% Qn : epsilon_channel-noise
        for iwl=1:length(wl) %loop over number of wavelengths
            lst=find(Cwl==iwl); %List of all wavelength <idx>
            Qn{iwl}=sparse(lst,lst,ones(size(lst)),NC,NC);
        end
        %%  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
        % c1 : couche cortex : huppert 2
        % c5 : couche skin   : huppert 1
        %%% PREMIERE IDEE : on projette les ondelettes sur les couches en reperant
        %%% quelle valeur appartient a quelle couche dans la longue matrice colonne
        %%% (on fait comme si les voxels dans la matrice colonne etaient adjacents dans l'image)
        %%% DEUXIEME IDEE
        %- on calcule des wavelets dans un plan et on incline ce plan selon celui
        %des sources et detecteurs
        
        Qp{1}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbO
        Qp{2}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbO
        Qp{3}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbR
        Qp{4}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbR
        
        %% On prepare les reconstructions :
        
        switch job.ReML_method
            case 0
                disp('code Huppert');
                meth = 'HUP';
                [lambda,beta_W,Stats]=nirs_run_DOT_REML(Y_t0,X*W',beta_prior,Qn,Qp);
                
                %Convert to the image domain and display
                beta = W'*beta_W;
                
            case 1
                disp('code spm_reml');
                meth = 'SPM';
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
            case 2
                disp('peudo inverse');
                meth = 'PInv';
                % Clement's version
                
                % M_c1 Mask for cortex
                M_c1_wl = [m_c1;m_c1];
                M_c1 = sparse(diag(M_c1_wl));
                
                Ybar = Y_t0;
                % TIKHONOV AMELIORE
                % % % % %         Ybar = sparse([Y_t0;zeros(3*size(X,2),1)]);
                
                %%%%%Attention pour le Xbar dans le cas de Tikhonov ameliore,la SVD n'est pas implemente....
                % % % % %         Xbar = sparse([log(X) log(X)*M_c1 log(X)*M_c1; sparse(1:3*size(X,2),1:3*size(X,2),ones(3*size(X,2),1),3*size(X,2),3*size(X,2))]);
                
                % begin SVD :
                % Chapitre 26 : p330
                % observation noise : hatOmega pour nous Qn
                Qinv =blkdiag(Qn{1}+Qn{2}); % hatsigma = Q
                Q = inv(Qinv);
                Xbar = Q.^(1/2)*full(X*M_c1);
                %         % custom SVD
                %         XbarN = Xbar'*Xbar;
                %         [v S v] = svd(XbarN,0);
                %         S       = sparse(S);
                %         s       = diag(S);
                %         j       = find(s*length(s)/sum(s) >= U & s >= T);
                %         v       = v(:,j);
                %         u       = spm_en(Xbar*v);
                %         S       = sqrt(S(j,j));
                %         % replace in full matrices
                %         %---------------------------------------------------------------------------
                %         j      = length(j);
                %         U      = sparse(M,j);
                %         V      = sparse(N,j);
                %         if j
                %             U(p,:) = u;
                %             V(q,:) = v;
                %         end
                %
                %         %%% voir les valeurs propres
                %         %         i=1;
                %         %         while S(i,i)>0 && i<=40
                %         %             disp([int2str(i) ' ieme valeur propre de S : ' num2str(S(i,i))])
                %         %             i = i+1;
                %         %         end
                %         %%%
                %         Vbar = S*V';
                %         Xbar = Vbar;
                %         % end SVD
                
                % TIKHONOV AMELIORE
                % - Beta contient omega_space omega et beta_prior : non code
                % car n'intervient pas dans la regularisation
                % - ebar DE MEME, par contre on definit la matrice des covariances
                % des erreurs :
                %         coef = 0.1;%%%%% moyen de calculer ca sur les images ??????????????
                % of course : idee : en pratique surtout au niveau des interfaces,
                % peut etre sortir l'info de ci_ fournie par newsegment puisque
                % c'est des cartes de probabilite....
                %         Qs=sparse(1:2*Nvx,1:2*Nvx,coef*ones(2*Nvx,1),2*Nvx,2*Nvx); % omega_space
                %Set up the extended covariance model by concatinating the measurement
                %and parameter noise terms and spatial prior
                %         Q =blkdiag(Qn{1}+Qn{2},Qs,Qp{1}+Qp{2}+Qp{3}+Qp{4},Qp{1}+Qp{2}+Qp{3}+Qp{4});
                
                
                % on applique ensuite la formule de l'inversion :
                %%%cas de la svd
                %         Beta_estimate = (Xbar'*Xbar) \ (Xbar'*U'*Ybar);
                
                %%% pour sauver de l'espace memoire: from Philippe
                %%%%%% je pense qu'on n'a pas le droit d'utiliser le slash a la
                %%%%%% place de pinv...
                %%% si l'on voulais prendre en compte les covariances, se rapporter
                %%% au fichier TeX :
                clear M_c1 M_c1_wl W X Xmc Yb8i_c1 Yb8i_c5 beta_prior m_c1
                alpha =1;
                pinvXX = pinv(Xbar'*Xbar + alpha*eye(size(Xbar,2)));
                YY = (Xbar'*Ybar);
                %save the inverse
                save([dir_in '\pinvXX.mat'],'pinvXX','-v7.3');
                save([dir_in '\YY.mat'],'YY','-v7.3');
                %free up some memory
                clear Xbar;
                
                Beta_estimate = pinvXX*YY;
                % TIKHONOV AMELIORE
                %         beta = Beta_estimate((size(Y_t0,1)+2*size(Xbar,2)+1):size(Beta_estimate,1),1);
                beta = Beta_estimate;
        end
        %%
        % %Convert the Stats
        % %%% juste les t stats
        % Stats.tstat.Cov_beta = W'*Stats.tstat.Cov_beta*W;
        % Stats.tstat.t=beta./sqrt(diag(Stats.tstat.Cov_beta));
        % Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);
        
        V = spm_vol(cs.segR);
        %Now, display the results
        beta_4d = reshape(full(beta),[V.dim 2]);
        
        beta_HbO = beta_4d(:,:,:,1);
        beta_HbR = beta_4d(:,:,:,2);
        
        
    case 2 % ON FAIT POUR CHAQUE LONGUEUR D'ONDE SEPAREMENT
        
        %% Y
        fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
        t0 =3000;% on choisit au pif un point temporel
        Y_t0 = fnirs.d(t0,Cmc)';
        
        %% X
        % p330
        % on prend pour \Omegachapeau I en premiere approximation
        Kmat = load(fullfile(dir_in,'sens.mat'));
        Xmc = Kmat.sens;
        NCdemi = size(Xmc,1)/2;
        %         X = sparse([[Xmc(1:NCdemi,:)*415.5 Xmc(1:NCdemi,:)*2141.8];[Xmc(NCdemi+1:end,:)*1008.0 Xmc(NCdemi+1:end,:)*778.0]]);
        Xwlt{1} = sparse([Xmc(1:NCdemi,:)*415.5 Xmc(1:NCdemi,:)*2141.8]);
        Xwlt{2} = sparse([Xmc(NCdemi+1:end,:)*1008.0 Xmc(NCdemi+1:end,:)*778.0]);
        
        for iwl=1:2
            Xwl = Xwlt{1,iwl};
            
            %% Qn : Covariance pour la longueur d'onde
            lst=(1:NCdemi);
            Qn=sparse(lst,lst,ones(size(lst)),NCdemi,NCdemi);
            %%  Qp : Covariance components for the parameters (4 total- 2 per HbO/HbR {layer 1; layer II})
            Qp{1}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbO
            Qp{2}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbO
            Qp{3}=sparse(Yb8i_c5,Yb8i_c5,zeros(length(Yb8i_c5),1),2*Nvx,2*Nvx);  %Skin layer- HbR
            Qp{4}=sparse(Yb8i_c1,Yb8i_c1,ones(length(Yb8i_c1),1),2*Nvx,2*Nvx);  %Brain layer- HbR
            
            %% On prepare les reconstructions :
            switch job.ReML_method
                case 0 % HUP
                    disp('A CODER')
                case 1 % SPM
                    disp('A CODER')
                case 2
                    disp('peudo inverse');
                    meth = 'PInv';
                    % Clement's version
                    
                    % % % % %         Ybar = sparse([Y_t0;zeros(3*size(X,2),1)]);
                    % M_c1 Mask for cortex
                    M_c1_wl = [m_c1;m_c1];
                    M_c1 = sparse(diag(M_c1_wl));
                    % % % % %         Xbar = sparse([log(X) log(X)*M_c1 log(X)*M_c1; sparse(1:3*size(X,2),1:3*size(X,2),ones(3*size(X,2),1),3*size(X,2),3*size(X,2))]);
                    Ybar = Y_t0((iwl-1)*NCdemi+(1:NCdemi));
                    
                    % begin SVD :
                    % Chapitre 26 : p330
                    % observation noise : hatOmega plur nous Qn
                    Qinv =blkdiag(Qn); % hatsigma = Q
                    Q = inv(Qinv);
                    Xbar = sparse(Q.^(1/2)*full(Xwl*M_c1));
                    
                    %         % custom SVD
                    %         XbarN = Xbar'*Xbar;
                    %         [v S v] = svd(XbarN,0);
                    %         S       = sparse(S);
                    %         s       = diag(S);
                    %         j       = find(s*length(s)/sum(s) >= U & s >= T);
                    %         v       = v(:,j);
                    %         u       = spm_en(Xbar*v);
                    %         S       = sqrt(S(j,j));
                    %         % replace in full matrices
                    %         %---------------------------------------------------------------------------
                    %         j      = length(j);
                    %         U      = sparse(M,j);
                    %         V      = sparse(N,j);
                    %         if j
                    %             U(p,:) = u;
                    %             V(q,:) = v;
                    %         end
                    %
                    %         %%% voir les valeurs propres
                    %         %         i=1;
                    %         %         while S(i,i)>0 && i<=40
                    %         %             disp([int2str(i) ' ieme valeur propre de S : ' num2str(S(i,i))])
                    %         %             i = i+1;
                    %         %         end
                    %         %%%
                    %         Vbar = S*V';
                    %         Xbar = Vbar;
                    %         % end SVD
                    
                    % TIKHONOV AMELIORE
                    % Beta contient omega_space omega et beta_prior : pour Tikhonov pas
                    % besoin de le definir puisque c'est nul...
                    % Betabar = sparse(,,,beta_prior);
                    
                    % ebar DE MEME, par contre on definit la matrice des covariances
                    % des erreurs :
                    %         coef = 0.1;%%%%% moyen de calculer ca sur les images ??????????????
                    % of course : idee : en pratique surtout au niveau des interfaces,
                    % peut etre sortir l'info de ci_ fournie par newsegment puisque
                    % c'est des cartes de probabilite....
                    %         Qs=sparse(1:2*Nvx,1:2*Nvx,coef*ones(2*Nvx,1),2*Nvx,2*Nvx); % omega_space
                    %Set up the extended covariance model by concatinating the measurement
                    %and parameter noise terms and spatial prior
                    %         Q =blkdiag(Qn{1}+Qn{2},Qs,Qp{1}+Qp{2}+Qp{3}+Qp{4},Qp{1}+Qp{2}+Qp{3}+Qp{4});
                    
                    
                    % on applique ensuite la formule de l'inversion :
                    %%%cas de la svd
                    %         Beta_estimate = (Xbar'*Xbar) \ (Xbar'*U'*Ybar);
                    
                    %%% pour sauver de l'espace memoire: from Philippe
                    %%%%%% je pense qu'on n'a pas le droit d'utiliser le slash a la
                    %%%%%% place de pinv...
                    %%% si l'on voulais prendre en compte les covariances, se rapporter
                    %%% au fichier TeX :
%                     clear M_c1 M_c1_wl W X Xmc Yb8i_c1 Yb8i_c5 beta_prior m_c1
                    alpha =1;
                    XX =Xbar'*Xbar;
                    YY = (Xbar'*Ybar);
                    clear Xbar;
                    XXLI = sparse(XX + eye(size(XX,2)));
                    pinvXX = XXLI \ eye(size(XXLI,1));
                        %                 pinvXX = pinv(XX + eye(size(XX,2)));
                    
                    %save the inverse
                     save([dir_in '\pinvXX_' int2str(iwl) '.mat'],'pinvXX','-v7.3');
                    save([dir_in '\YY_' int2str(iwl) '.mat'],'YY','-v7.3');
                    %free up some memory
                    clear Xbar;
                    
                    Beta_estimate = pinvXX*YY;
                    % TIKHONOV AMELIORE
                    %         beta = Beta_estimate((size(Y_t0,1)+2*size(Xbar,2)+1):size(Beta_estimate,1),1);
                    beta = Beta_estimate;
            end
            %%
            % t-stats :
            %calcul de C(beta|Y) ::: pour l'instant on laisse tomber...
            %             pinv(Xwl'*Q*Xwl
            %             Stats.tstat.Cov_beta = Stats.tstat.Cov_beta;
            %             Stats.tstat.t=beta./sqrt(diag(Stats.tstat.Cov_beta));
            %             Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);
            
            V = spm_vol(cs.segR);
            %Now, display the results
            beta_4d = reshape(full(beta),[V.dim 2]);
            
            beta_HbO = beta_4d(:,:,:,1);
            beta_HbR = beta_4d(:,:,:,2);
        end
end

% creation de nii :
V_O = struct('fname',fullfile(dir_in,['HbO_' meth '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);

V_O = spm_create_vol(V_O);
V_O = spm_write_vol(V_O, beta_HbO);

V_R = struct('fname',fullfile(dir_in,['HbR_' meth '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);

V_R = spm_create_vol(V_R);
V_R = spm_write_vol(V_R, beta_HbR);

% superpositions : on cree des images semi transparentes rouges ou bleues
% sur les anatomiques (VOIR CHECKREG + CLIC DROIT)
spm_imcalc_ui({fullfile(dir_in,['HbO_' meth '.nii']);...
    cs.segR},...
    fullfile(dir_in,['HbO_anat_' meth '.nii']),...
    'i1+i2');

spm_imcalc_ui({fullfile(dir_in,['HbR_' meth '.nii']);...
    cs.segR},...
    fullfile(dir_in,['HbR_anat_' meth '.nii']),...
    'i1+i2');

out =1;
end