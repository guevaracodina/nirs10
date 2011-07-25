function out = nirs_run_configMC2(job)
% Manages the writing of the configuration files for either MCX or tMCimg
% FORMAT nirs_configMC_writeCFGfiles(jobW) = nirs_run_configMC(MC_nam,NIRSmat,NewDirCopyNIRS,mcim_cfg,MC_CUDAchoice,MC_configdir,MC_parameters)
% MC_nam         - Name of the MC simulation
% NIRSmat        - NIRS matrix
% NewDirCopyNIRS - New directory for NIRSmat if specified
% mcim_cfg       - Segmented image that will be binarised
% MC_CUDAchoice  - MC algorithm
% MC_configdir   - directory for configuration files
% MC_parameters  - parameters for the simulation
%_______________________________________________________________________
%
% Either the last image from which a ROI as been taken or the selected
% image can be processed. Whereas in the second case all the channels will
% be kept, in the first, the user can choose.
% The chosen image is then transformed from anisotropic voxels space to
% isotropic voxels space and transformed to binary image. All the positions
% are also corrected : the positions are first changed from mm to voxels.
% The probe is then fitted in the rigid space with isotropic voxels.
% Last, directions are calculated.
%
% BEWARE : directions must be in voxel space (one must apply rotations on
% directions calculated in mm). They must point towards the brain.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clement Bonnery

%%%%%%% WARNING: possible erreur lorsqu on choisit latestROI que l on se
%%%%%%% trouve dans le bon repertoire... (buildROI2 laisse la possibilite de creer un sous directory peut etre que c est pas une si bonne idee !!!!! a verifier !!!!!)

for Idx=1:size(job.NIRSmat,1)
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        
        if ~isfield(NIRS,'Cs'), NIRS.Cs ={}; end
        if isfield(NIRS.Cs,'mcs')
            i_cs = size(NIRS.Cs.mcs,2)+1;
        else
            i_cs=1;
        end
        
        % cs current simulation
        cs ={};
        G={};%%% juste pour le code
        P={};%%% juste pour le code
        if isfield(job.mcim_cfg,'mcim_in')% image segmentee de l'anatomique de base
            G.roi =0;% image choisie
            G.seg = job.mcim_cfg.mcim_in{:};
        else
            G.roi =1; %ROI from temp
            try
                G.seg = NIRS.Cs.temp.segR{:};%% ici segR pour ROI
            catch
                G.seg = NIRS.Cs.temp.segR;%% ici segR pour ROI
            end
        end
        G.seg_p = G.seg(1:max(strfind(G.seg,filesep))-1);
        
        cs.alg = job.MC_CUDAchoice;%1=MCX ; 2=tMCimg
        if cs.alg==1, alg_nam='MCX'; else alg_nam ='tMC';end
        
        cs.par = job.MC_parameters;
        cs.NSinit = NIRS.Cf.H.S.N;
%         cs.Ckpt = NIRS.Cs.temp.Ckpt;
        cs.Pkpt = NIRS.Cs.temp.Pkpt;
        cs.NSkpt = NIRS.Cs.temp.NSkpt;
        cs.NDkpt = NIRS.Cs.temp.NDkpt;
        
        P.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
        P.Pfp_rmm = NIRS.Cs.temp.Pfp_roi_rmm;
        P.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
        P.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
        
        if ~G.roi && isfield(NIRS.Cf.H.P,'void')
            Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
        else
            Pvoid = {};
        end
        cs.Pvoid = Pvoid;
        
        Sr = cs.par.radiis * ones(cs.NSkpt,1);
        cs.Dr = cs.par.radiid * ones(cs.NDkpt,1);
        cs.r = [Sr' cs.Dr' zeros(1,size(cs.Pkpt,1) -(cs.NSkpt+cs.NDkpt))];
        
        daate = strrep(datestr(now),':','-');
        
        %%% on definit un nom
        if job.pve_cfg==1 %%%% Thresholded BOLD image is concidered as layer 6
            %%%% it might be also a perturbation included by user... any use ??
            % on cherche toutes les sessions qui pourraient mener a un
            % calcul de PVE
            last = size(NIRS.Dt.fir.pp,2);
            NSess = size(NIRS.Dt.fir.pp(1,last).p,1);
            
            for iSess =1:NSess
                arun=1;
                
                if arun==0
                    clear xSPM_boldmask xSPM_boldmask_sorted
                    % on genere une simulation MonteCarlo par contraste BOLD ///
                    [dummy,xSPM,SPM] = spm_results_ui('Setup');
                    [dir,dummy] = fileparts(job.NIRSmat{:});
                    save(fullfile(dir,'xSPM.mat'),'xSPM');
                    NIRS.Dt.fmri.xSPM = fullfile(dir,'xSPM.mat');
                    save(fullfile(dir,'SPM.mat'),'SPM');
                    NIRS.SPM = fullfile(dir,'SPM.mat');
                else%%%% cas ou on a deja la matrice xSPM
                    load('D:\Users\Clément\Projet_ReML\donnees\test_MCX\S531\xSPM.mat');
                end
                
                cs_title = xSPM.title;
                csn = [alg_nam '_' cs_title '_' daate];
                cs.dir =fullfile(G.seg_p,csn);
                if ~exist(cs.dir,'dir')%Directory for configuration files
                    mkdir(cs.dir);
                end
                
                try
                    jobP.boldmask = xSPM.Z;
                    jobP.XYZmm = xSPM.XYZmm;
                    jobP.M = xSPM.M;
                    jobP.T1seg = G.seg;%NIRS.Dt.ana.T1seg;
                    outP = nirs_MCsegment_PVE(jobP);
                    G.seg = outP;
                    cs.pve_cfg = 1;
                    
                    % fin de config et ecriture
                    jobF.G = G;
                    jobF.G.dir = cs.dir;
                    jobF.G.voxelSize =cs.par.voxelSize;
                    jobF.G.alg = cs.alg;
                    jobF.P = P;
                    outF = nirs_configMC_toMCframe(jobF);%images,positions and directions
                    cs.segR = outF.G.segR;
                    cs.b8i = outF.G.b8i;
                    cs.ROIlimits = outF.G.ROIlimits;
                    cs.P = outF.P;
                    cs.Pfp_rmiv = outF.P.Pfp_rmiv;
                    cs.Pwd_rmm = outF.P.Pwd_rmm;
                    cs.nummed = 11;
                    
                    NIRS.Cs.mcs{i_cs} = cs;
                    NIRS.Cs.n{i_cs} = csn;
                    
                    jobW = cs;
                    jobW.wl_dev = NIRS.Cf.dev.wl;
                    nirs_configMC_writeCFGfiles2(jobW);
                catch
                    disp(['PVE failed for session ' iSess]);
                end
            end
        else
            csn = [alg_nam '_' daate];
            cs.dir =fullfile(G.seg_p,csn);
            %Directory for configuration files
            if ~exist(cs.dir,'dir')
                mkdir(cs.dir);
            end
            cs.pve_cfg = 0;
            
            % fin de config et ecriture
            jobF.G = G;
            jobF.P = P;
            outF = nirs_configMC_toMCframe(jobF);%images,positions and directions
            cs.segR = outF.G.segR;
            cs.b8i = outF.G.b8i;
            cs.P = outF.P;
            cs.Pfp_rmiv = outF.P.Pfp_rmiv;
            cs.Pwd_rmm = outF.P.Pwd_rmm;
            cs.nummed =6;
            
            NIRS.Cs.mcs{i_cs} = cs;
            NIRS.Cs.n{i_cs} = csn;
            
            jobW = cs;
            jobW.wl_dev = NIRS.Cf.dev.wl;
            nirs_configMC_writeCFGfiles2(jobW);
        end
        
        newNIRSlocation = fullfile(cs.dir,'NIRS.mat');
        save(newNIRSlocation,'NIRS');
        job.NIRSmat{Idx,1} = newNIRSlocation;
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not run MonteCarlo configuration for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;

% %%%% old version
% function out = nirs_run_configMC(job)
% % Manages the writing of the configurqtion files for either MCX or tMCimg
% % FORMAT nirs_configMC_writeCFGfiles(jobW) = nirs_run_configMC(MC_nam,NIRSmat,NewDirCopyNIRS,mcim_cfg,MC_CUDAchoice,MC_configdir,MC_parameters)
% % MC_nam         - Name of the MC simulation
% % NIRSmat        - NIRS matrix
% % NewDirCopyNIRS - New directory for NIRSmat if specified
% % mcim_cfg       - Segmented image that will be binarised
% % MC_CUDAchoice  - MC algorithm
% % MC_configdir   - directory for configuration files
% % MC_parameters  - parameters for the simulation
% %_______________________________________________________________________
% %
% % Either the last image from which a ROI as been taken or the selected
% % image can be processed. Whereas in the second case all the channels will
% % be kept, in the first, the user can choose.
% % The chosen image is then transformed from anisotropic voxels space to
% % isotropic voxels space and transformed to binary image. All the positions
% % are also corrected : the positions are first changed from mm to voxels.
% % The probe is then fitted in the rigid space with isotropic voxels.
% % Last, directions are calculated.
% %
% % BEWARE : directions must be in voxel space (one must apply rotations on
% % directions calculated in mm). They must point towards the brain.
% %_______________________________________________________________________
% % Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% % Clement Bonnery
%
% %%%%%%% WARNING: possible erreur lorsqu on choisit latestROI que l on se
% %%%%%%% trouve dans le bon repertoire... (buildROI2 laisse la possibilite de creer un sous directory peut etre que c est pas une si bonne idee !!!!! a verifier !!!!!)
%
% prmts.gmPpties_l1 =    [0.0186   11.1   0.9   1.4];
% prmts.wmPpties_l1 =    [0.0186   11.1   0.9   1.4];
% prmts.csfPpties_l1 =   [0.0026   0.10   0.9   1.4];
% prmts.skullPpties_l1 = [0.0136   8.60   0.9   1.4];
% prmts.scalpPpties_l1 = [0.0191   6.60   0.9   1.4];
%
% prmts.gmPpties_l2 =    [0.0178   12.5   0.9   1.4];
% prmts.wmPpties_l2 =    [0.0178   12.5   0.9   1.4];
% prmts.csfPpties_l2 =   [0.0004   0.10   0.9   1.4];
% prmts.skullPpties_l2 = [0.0101   10.0   0.9   1.4];
% prmts.scalpPpties_l2 = [0.0159   8.00   0.9   1.4];
%
% %Overwrite to fix path
% job.MC_configdir = 'MC';
%
% for Idx=1:size(job.NIRSmat,1)
%     %Load NIRS.mat information
%     try
%         NIRS = [];
%         load(job.NIRSmat{Idx,1});
%
%         if ~isfield(NIRS,'Cs')
%             NIRS.Cs ={};
%         end
%
%         %%%% Thresholded BOLD image is concidered as layer 6
%         %%%% it might be also a perturbation included by user... any use ??
%          if job.pve_cfg==1
%              % on cherche toutes les sessions qui pourraient mener a un
%              % calcul de PVE
%              Nsess = size(NIRS.Dt.fir.pp.p,1);
%
%              % on genere une simulation MonteCarlo par contraste BOLD ///
%             [hReg,xSPM,SPM] = spm_results_ui('Setup');
%             [dir,dummy] = fileparts(job.NIRSmat{:});
%             save(fullfile(dir,'SPM.mat'),'SPM');
%             NIRS.Dt.fmri.xSPM = fullfile(dir,'SPM.mat');
%             save(fullfile(dir,'xSPM.mat'),'xSPM');
%             NIRS.SPM = fullfile(dir,'xSPM.mat');
%
%             xSPM_boldmask(1,:) = xSPM.Z;
% %             xSPM_boldmask(2:4,:) = xSPM.XYZmm;
%             xSPM_boldmask(2,:) =(1:length(xSPM.Z));
%
%             xSPM_boldmask_sorted = sortrows(xSPM_boldmask')';
%             i=1;
%             while sum(xSPM_boldmask_sorted(1,end+1-i:end)) < 0.8*sum(xSPM_boldmask(1,:)), i = i+1;end
%             level = xSPM_boldmask_sorted(1,i);
%
%             boldmask = zeros(1,length(xSPM.Z));
%             boldmask(xSPM_boldmask(1,:)>level)=1;
%
%             parameters.bold.boldmask = boldmask;
%             parameters.bold.XYZmm = xSPM.XYZmm;
%          end
%
%         if isfield(NIRS.Cs,'mcs')
%             i_cs = size(NIRS.Cs.mcs,2)+1;
%             if ~sum(strcmp(NIRS.Cs.n,job.MC_nam))
%                 csn = job.MC_nam;
%             else
%                 csn = [job.MC_nam strrep(datestr(now),':','-') '_' xSPM.title];
%             end
%         else
%             i_cs=1;
%             csn = job.MC_nam;
%         end
%
%         % cs current simulation
%         cs.alg = job.MC_CUDAchoice;%1=MCX ; 2=tMCimg
%         cs.dir = [job.MC_configdir csn];
%         cs.par = job.MC_parameters;
%         cs.NSinit = NIRS.Cf.H.S.N;
%
%         if isfield(job.mcim_cfg,'mcim_in')% image segmentee de l'anatomique de base
%             roi =0;% image choisie
%             cs.seg = job.mcim_cfg.mcim_in{:};
%             cs.p = cs.seg(1:max(strfind(cs.seg,filesep))-1);
%             cs.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
%             cs.Pfp_rmm = NIRS.Cs.temp.Pfp_roi_rmm;
%             cs.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
%             cs.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
%             cs.Pkpt = NIRS.Cs.temp.Pkpt;
%             cs.NSkpt = NIRS.Cs.temp.NSkpt;
%             cs.NDkpt = NIRS.Cs.temp.NDkpt;
%         else
%             roi =1;
%             try
%                 cs.seg = NIRS.Cs.temp.segR{:}; %ROI from temp
%             catch
%                 cs.seg = NIRS.Cs.temp.segR;
%             end
%             cs.p = cs.seg(1:max(strfind(cs.seg,filesep))-1);
%             cs.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
%             cs.Pfp_rmm = NIRS.Cs.temp.Pfp_roi_rmm;
%             cs.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
%             cs.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
%             cs.Pkpt = NIRS.Cs.temp.Pkpt;
%             cs.NSkpt = NIRS.Cs.temp.NSkpt;
%             cs.NDkpt = NIRS.Cs.temp.NDkpt;
%         end
%
%         if ~exist(fullfile(cs.p,[job.MC_configdir csn]),'dir')%Directory for configuration files
%             mkdir(fullfile(cs.p,[job.MC_configdir csn]));
%         end
%
%         NIRS.Cs.mcs{i_cs} = cs;
%         NIRS.Cs.n{i_cs} = csn;
%         save(job.NIRSmat{1,1},'NIRS');
%
%         parameters.voxelSize = job.MC_parameters.voxelSize;
%         % transform image from anisotropic voxels space to isotropic voxels space
%         jobRS.image_in = {cs.seg};
%         jobRS.out_dir = fullfile(cs.p,[job.MC_configdir csn]);
%         jobRS.out_dim = [1 1 1];
%         jobRS.out_vxsize = parameters.voxelSize;
%         jobRS.out_dt = 'same';
%         jobRS.out_autonaming = 0;
%         jobRS.out_prefix = 'prefix';
%         outRS =nirs_resize(jobRS);
%
%         clear NIRS
%
%         V_rmiv = spm_vol(outRS);
%         Y_rmiv = spm_read_vols(V_rmiv);
%         Y8_rmiv = uint8(Y_rmiv);
%
%         load(job.NIRSmat{1,1});
%         cs = NIRS.Cs.mcs{i_cs};
%
%         V = spm_vol(cs.seg);
%         % dim = V.dim;% on est ici avec l'ancienne image qui permet d'obtenir les bonnes positions en voxels
%         inv_mat = spm_imatrix(V.mat);
%
%         % Transform also P positions and directions
%         NP = size(cs.Pfp_rmv,2);
%         % Positions : Transform MNI mm -> MNI isotropic voxels
%         for i=1:size(cs.Pfp_rmm,2)
%         Pfp_ancienne_rmv(:,i) = V.mat\[cs.Pfp_rmm(:,i);1];
%         Pfp_ancienne_rmiv(:,i) = abs(inv_mat(7:9)').*(Pfp_ancienne_rmv(1:3,i)/parameters.voxelSize);
%         end
%
%         Pfp_ancienne_rmiv = round(Pfp_ancienne_rmiv);
%
%         % Directions
%         Pd_rmm = cs.Pp_rmm - cs.Pp_c1_rmm;
%         Pwd_rmm = zeros(3,NP);
%         for iP=1:NP
%             temp_dir = Pd_rmm(:,iP);
%             lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
%             Pwd_rmm(:,iP) = temp_dir/lgth;
%         end
%
%           if job.MC_CUDAchoice==1
%             jobF.Pp_rmm = cs.Pp_rmm;
%             jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
%             jobF.NP = NP;
%             jobF.image_in = {outRS};
%             jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
%             jobF.lby = 'configMC_MCX';
%             outF = nirs_fit_probe(jobF);
%             Pfp_ancienne_rmiv = outF{1};
%           elseif job.MC_CUDAchoice==2 % pour tMC, les points doivent etre dans le volume (peut etre est-ce juste une question de direction de lq propagation...)
%               jobF.Pp_rmm = cs.Pp_rmm;
%             jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
%             jobF.NP = NP;
%             jobF.image_in = {outRS};
%             jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
%             jobF.lby = 'configMC_tMC';
%             outF = nirs_fit_probe(jobF);
%             Pfp_ancienne_rmiv = outF{1};
%           end
%
%         % on prepare la sauvegarde de cs
%         cs.Pfp_rmiv = Pfp_ancienne_rmiv;
%         cs.Pwd_rmm = Pwd_rmm;
%
%         % 8bits .bin image
%         [dummy,id,dummy2] = fileparts(V_rmiv.fname);
%         dim_rmiv = V_rmiv.dim;
%         n_b8i = ['vol8bit_' id '.bin'];
%
%         cs.segR = outRS;
%         cs.b8i = fullfile(cs.p,[job.MC_configdir csn],n_b8i);
%
%         cs.numTimeGates = job.MC_parameters.numTimeGates;
%         cs.deltaT = job.MC_parameters.deltaT;
%         %Path for configuration files
%         NIRS.Cs.mcs{i_cs} = cs;
%         save(job.NIRSmat{1,1},'NIRS');
%
%         fid = fopen(fullfile(cs.p,[job.MC_configdir csn],n_b8i),'wb');
%         fwrite(fid, Y8_rmiv, 'uint8');
%         fclose(fid);
%
%         if ~roi && isfield(NIRS.Cf.H.P,'void')
%             Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
%         else
%             Pvoid = {};
%         end
%
%         Sr = job.MC_parameters.radiis * ones(cs.NSkpt,1);
%         Dr = job.MC_parameters.radiid * ones(cs.NDkpt,1);
%
%         parameters.nphotons = job.MC_parameters.nphotons;
%         parameters.seed = job.MC_parameters.seed;
%         parameters.modulationFreq = job.MC_parameters.modulationFreq;
%
%         parameters.numTimeGates = job.MC_parameters.numTimeGates;
%         parameters.deltaT = job.MC_parameters.deltaT;
%
%         % Create a .cfg or .inp file for each optode and each wavelength
%         for iwl = 1:size(NIRS.Cf.dev.wl,2)
%             if NIRS.Cf.dev.wl(iwl) == 830 %830
%                 parameters.gmPpties           = prmts.gmPpties_l1;
%                 parameters.wmPpties           = prmts.wmPpties_l1;
%                 parameters.csfPpties          = prmts.csfPpties_l1;
%                 parameters.skullPpties        = prmts.skullPpties_l1;
%                 parameters.scalpPpties        = prmts.scalpPpties_l1;
%                 parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l1+ prmts.gmPpties_l1;
%             elseif NIRS.Cf.dev.wl(iwl) == 690 %690
%                 parameters.gmPpties           = prmts.gmPpties_l2;
%                 parameters.wmPpties           = prmts.wmPpties_l2;
%                 parameters.csfPpties          = prmts.csfPpties_l2;
%                 parameters.skullPpties        = prmts.skullPpties_l2;
%                 parameters.scalpPpties        = prmts.scalpPpties_l2;
%                 parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l2+ prmts.gmPpties_l2;
%             end
%
%
%             %%%% attention la on est avec l'image resizee...........
%             jobW.algo = job.MC_CUDAchoice;
%             jobW.n_b8i = n_b8i;
%             jobW.mc_dir = fullfile(cs.p,[job.MC_configdir csn]);
%             jobW.ROIlimits = [1 1 1; dim_rmiv];
%
%             jobW.parameters = parameters;
%             jobW.NSinit = cs.NSinit;
%             jobW.NS = cs.NSkpt;
%             jobW.ND = cs.NDkpt;
%             jobW.Pvoid = Pvoid;
%
%             for i=1:size(cs.Pwd_rmm,2)
%                 Pwd_rmiv(:,i) = V_rmiv.mat(1:3,1:3)\cs.Pwd_rmm(:,i);
% %                 Pfp_tMC(:,i) = V_rmiv.mat*[Pfp_ancienne_rmiv(:,i);1];
%             end
%
%             if job.MC_CUDAchoice==1
%                 %%% MCX en voxel /////!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 % comme on travaille sur les donnees sur des voxels isotropiques de
%                 % 1mm, il suffit de diviser par la taille des voxels (dans le cas
%                 % ou on serait avec des voxels de taille custom, il faudrait
%                 % diviser la taille des voxels voulu par la taille courante)
%                 P.p =Pfp_ancienne_rmiv;%/(parameters.voxelSize); DEFINITIFFFFFF
%                 P.wd = -Pwd_rmiv(1:3,:); %to point toward the inside of brain
%
%             elseif job.MC_CUDAchoice==2
%                 % MonteCarlo in a particular frame. Positions must be in mm but the
%             % origin is the same as the origin of the voxel frame (these positions
%             % don't respect SPM conventions) %%definitif ////
%              P.p = parameters.voxelSize*Pfp_ancienne_rmiv;
%              P.wd = -V_rmiv.mat(1:3,1:3)*Pwd_rmm;
%             end
%
%             P.Pkpt = cs.Pkpt;
%
%             P.r = [Sr' Dr' zeros(1,NP -(cs.NSkpt+cs.NDkpt))];
%             P.Sr = job.MC_parameters.radiis;
%             P.Dr = job.MC_parameters.radiid;
%             jobW.P =P;
%             jobW.wl = NIRS.Cf.dev.wl(iwl);
%
%             out_void = nirs_configMC_writeCFGfiles(jobW);
%         end
%         newNIRSlocation = fullfile(cs.p,[job.MC_configdir csn],'NIRS.mat');
%         save(newNIRSlocation,'NIRS');
%         job.NIRSmat{Idx,1} = newNIRSlocation;
%     catch exception
%         disp(exception.identifier);
%         disp(exception.stack(1));
%         disp(['Could not run MonteCarlo configuration for subject' int2str(Idx)]);
%     end
% end
% out.NIRSmat = job.NIRSmat;