function out = nirs_run_inverse_tikhonov(job)
% Clément Bonnéry May 2011

% Le systeme qu'on resout est le suivant :
% on est en un UNIQUE POINT TEMPOREL donc
%              ---------------------
%
% / Y = Xsens*Dmua(vx) + epsilon_channel-noise
% \ Dmua(vx) = [ext]*[DHbO(vx);DHbR(vx)]
%
% on cherche a minimiser :
% argmin(||Y - Xsens*[ext]*beta||^2 + alpha*||beta - beta0||^2)
%
% Y                     = un instant donne dans le fichier .nirs
% Xsens                 = matrice sensitivite avec un point temporel
% [ext]                 = matrice des coefficients d extinctions
% Dmua                  = mu_a dans les voxels
% DHbO\R (beta)         = variation des concentrations dans les voxels
% alpha                 = hyperparametre de beta (DHbR\O)
% epsilon_channel-noise = bruit dans les canaux
if isfield(job.psel_choice,'all_points_downsampled')
    pmethod = 1;
    downfreq = job.psel_choice.all_points_downsampled.downsample_freq;
else
    pmethod = 0;
    %temporal points
    temp_pts = job.psel_choice.specific_points.temp_pts;
end
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
        
        %current tomo reconstruction
        ctm = {};
        % b8i is...
        % segR  is...
        % b8i and segR have the same orientations as there is no
        % permutation (C and Matlab conventions in nirs_run_configMC)
        
        %resize segR so that it has same size as sensitivity matrix
        VsegR = spm_vol(cs.segR);
        clear jobR
        jobR.image_in ={cs.segR};
        jobR.out_dir = cs_dir;
        jobR.out_dim = [1 1 1];
        jobR.out_dt = 'same';
        jobR.out_vxsize = job.sens_vxsize;
        jobR.out_autonaming = 0;
        temp.segRR = nirs_resize(jobR);
        %put segRR in memory, calling it YsegRR, and used only for masks
        VsegRR = spm_vol(temp.segRR);
        YsegRR = spm_read_vols(VsegRR);
        YsegRR = reshape(YsegRR,[1 prod(VsegRR.dim)]);
        %masks c1 and c5
        m_c1 = zeros(size(YsegRR));
        m_c1(YsegRR==1)=1;% mask for GM
        m_c5 = zeros(size(YsegRR));
        m_c5(YsegRR==5)=1;% mask for skin
        clear YsegRR;
        
        % Source detector pairs....
        C_cs = cs.C; %former Cmc
        NC_cs = length(C_cs); %Total number of measurements
        %load sensitivity matrix
        load(fullfile(cs_dir,'sens.mat'));
        Xmc = sens;
        clear sens
        % Need to downsample sensitivity matrix
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
            fname =  nirs_resize(jobR);
            clear Xmci Vmc
            Vmc = spm_vol(fname);
            Ymc = spm_read_vols(Vmc);
            %
            XmcR(Ci,:) = reshape(Ymc,[1 prod(Vmc.dim)]);
        end
        Xmc = XmcR;
        clear XmcR Vmc Ymc
        Nvx = size(Xmc,2);
        NC2mi = NC_cs/2; %only for 2 wavelengths; number of pairs
        
        ext1 = GetExtinctions(NIRS.Cf.dev.wl(1,1));
        ext2 = GetExtinctions(NIRS.Cf.dev.wl(1,2));
        Ext = [ext1(1,1) ext1(1,2) ; ext2(1,1) ext2(1,2)];
        %regularization parameter
        alpha =job.alpha;
        
        Xwl{1} = sparse(Xmc(1:NC2mi,:)); %sensitivity matrix for 1st wavelength
        Xwl{2} = sparse(Xmc(NC2mi+1:end,:));
        
        Msk = sparse(diag(m_c1+m_c5)); % M Mask for cortex and skin
        %need to load data files prior to OD->concentration
        pdata = NIRS.Dt.fir.pp;
        wh = length(pdata);
        while wh
            [dummy fil1] = fileparts(pdata(wh).p{1});
            if strcmp(fil1(1),'h') || strcmp(fil1(1:2),'bh')
                wh = wh-1;
            else
                pdata = pdata(wh).p;
                wh = 0;
            end
        end
        
        for ifnirs=1:length(pdata)
            Y = fopen_NIR(pdata{ifnirs},NIRS.Cf.H.C.N)';           
            %creates Tmrs : tomographical reconstruction
            if ~isfield(NIRS,'Tm')
                NIRS.Tm ={};
            end
            if pmethod
                downstep = round(NIRS.Cf.dev.fs/downfreq);
                TR = downstep/downfreq/NIRS.Cf.dev.fs;               
                temp_pts = 1:downstep:size(Y,1);
            end
            %mkdir
            switch job.tikh_method
                case 0 % Tikhonov regularization : ancienne version
                    disp('Methode a l''ancienne')
                    ctm.alg ='anci';
                case 1 % Tikhonov regularization
                    disp('pseudo inverse');
                    ctm.alg = 'PInv';
                case 2 % avec wavelets...
                case 3 % Li et al extended Tikhonov regularization (est ce que ca a de l interet alors qu il va falloir trouver deux hyperparametres)
                case 4 % Tikhonov ameliore
            end
            daate = strrep(datestr(now),':','-');
            tm_dir = ['tm_' daate '_a' ctm.alg '_' int2str(job.sens_vxsize) 'mm'];
            ctm.p = fullfile(cs_dir,tm_dir);
            if ~exist(ctm.p,'dir'),mkdir(ctm.p);end
            if pmethod, save(fullfile(ctm.p,'TR.mat'),'TR'); end
            if isfield(NIRS.Tm,'tmrs')
                itm = size(NIRS.Tm.tmrs,2)+1;
            else
                itm=1;
            end
            ctm.n = tm_dir;
            ctm.hyperparameter = alpha;
            ctm.sens_vxsize = job.sens_vxsize;
            
            [dum,namRR,extRR] = fileparts(temp.segRR);
            ctm.segRR = fullfile(ctm.p,[namRR extRR]);
            copyfile(temp.segRR,ctm.segRR);
            delete(temp.segRR);
            NIRS.Tm.tmrs{itm} = ctm;
            NIRS.Tm.n{itm} = ctm.n;
            
            for itp=1:length(temp_pts)
                %tic
                Y_t0 = Y(temp_pts(itp),C_cs)';
                Y_to{1} = Y_t0(1:NC2mi,1);
                Y_to{2} = Y_t0(NC2mi+1:end,1);
                
                %                 %%%% test fantome
                %                 Yt0 = load('Yt0.nirs','-mat');
                %                 Y_to{1} = Yt0.Yt0(1:size(Yt0.Yt0,1)/2,1);
                %                 Y_to{2} = Yt0.Yt0(size(Yt0.Yt0,1)/2+1:size(Yt0.Yt0,1),1);
                %                 %%%%
                
                % Dmua
                Dmua{1} = zeros(Nvx,1);
                Dmua{2} = zeros(Nvx,1);
                
                switch job.tikh_method
                    case 0 % Tikhonov regularization : ancienne version
                        for iwl=1:2
                            Ybar = Y_to{1,iwl};
                            X = Xwl{1,iwl};
                            Xbar = sparse(X*Msk); %??
                            
                            XX =Xbar'*Xbar;
                            YY = (Xbar'*Ybar);
                            clear Xbar;
                            XXLI = sparse(XX + alpha*eye(size(XX,2)));%eq19 huppert_2010_hierarchical
                            Dmua{iwl} = XXLI \ YY;
                        end
                    case 1 % Tikhonov regularization
                        for iwl=1:2
                            if itp==1
                                XX{iwl} =Xwl{iwl}'*Xwl{iwl};
                                %takes 27 GB of memory for 4e4 x 4e4 size
                                XXLI{iwl} = sparse(XX{iwl} + alpha*eye(size(XX{iwl},2))); % eq.19 huppert_2010_hierarchical
                            end
                            YY = (Xwl{iwl}'*Y_to{iwl});
                            Dmua{iwl} = XXLI{iwl} \ YY;
                        end
                    case 2 % avec wavelets...
                    case 3 % Li et al extended Tikhonov regularization (est ce que ca a de l interet alors qu il va falloir trouver deux hyperparametres)
                    case 4 % Tikhonov ameliore
                        % Interpretation bayesienne simple (regularization de Tikhonov avec des valeurs de covariances non nulles)
                        %                                 % Beta contient omega_space omega et beta_prior : pour Tikhonov pas
                        %                                 % besoin de le definir puisque c'est nul...
                        %                                 % Betabar = sparse(,,,beta_prior);
                        %
                        %                                 % ebar DE MEME, par contre on definit la matrice des covariances
                        %                                 % des erreurs :
                        %                                 %         coef = 0.1;%%%%% moyen de calculer ca sur les images ??????????????
                        %                                 % of course : idee : en pratique surtout au niveau des interfaces,
                        %                                 % peut etre sortir l'info de ci_ fournie par newsegment puisque
                        %                                 % c'est des cartes de probabilite....
                        %                                 %         Qs=sparse(1:2*Nvx,1:2*Nvx,coef*ones(2*Nvx,1),2*Nvx,2*Nvx); % omega_space
                        %                                 %Set up the extended covariance model by concatinating the measurement
                        %                                 %and parameter noise terms and spatial prior
                        %                                 %         Q =blkdiag(Qn{1}+Qn{2},Qs,Qp{1}+Qp{2}+Qp{3}+Qp{4},Qp{1}+Qp{2}+Qp{3}+Qp{4});
                        %
                        %
                        %                                 % on applique ensuite la formule de l'inversion :
                        %                                 %%%cas de la svd
                        %                                 %         Beta_estimate = (Xbar'*Xbar) \ (Xbar'*U'*Ybar);
                end
                
                beta = Ext \ [Dmua{1} Dmua{2}]'; %get HbO HbR from Mua
                betaR_HbO = reshape(full(beta(1,:)),VsegRR.dim);
                betaR_HbR = reshape(full(beta(2,:)),VsegRR.dim);
                %write nifti for DHbO DHbR
                str0 = gen_num_str(itp,4);
                V_O = nirs_create_vol(fullfile(ctm.p,['O_' str0 '.nii']),...
                    VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbO);
                V_R = nirs_create_vol(fullfile(ctm.p,['R_' str0 '.nii']),...
                    VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbR);
                
                %toc
                %disp(int2str(itp));
            end
        end
        save(job.NIRSmat{Idx,1},'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        try
            rmdir(ctm.p); %careful, the code might break here, so put breakpoint on line before
        end
        disp(['Could not run MonteCarlo reconstruction for subject' int2str(Idx)]);
    end
end
%clear out %careful, variable out was used with another meaning previously!
out.NIRSmat = job.NIRSmat;

