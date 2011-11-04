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
        
        cs.mu_subj = job.mu_subj;
        cs.par = job.MC_parameters;
        cs.NSinit = NIRS.Cf.H.S.N;
        cs.Ckpt = NIRS.Cs.temp.Ckpt;
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
        
        daate = strrep(datestr(now),':','-');
        
        % CONFIGURE SIMULATION %     

        if sum(isfield(job.pve_cfg,{'pve_bold','pve_asl','pve_anat'}))
            %%%% Thresholded BOLD image is considered as layer 6
            
            % on cherche toutes les sessions qui pourraient mener a un        
            % calcul de PVE
            last = size(NIRS.Dt.fir.pp,2);
            NSess = size(NIRS.Dt.fir.pp(1,last).p,1);
            do_results_ui = 1;
            cs.pve_cfg = 1;
            cs.nummed = 12;
            
            if isfield(job.pve_cfg,'pve_anat')
                %%%% it might be also a perturbation included by user, not
                %%%% necessarily taken from BOLD results
                do_results_ui = 0;
                NSess = 1; %size(job.pve_cfg.pve_anat)
            end
            
        elseif isfield(job.pve_cfg,'no_pve')           
            cs.pve_cfg = 0;          
            NSess = 1; % by default only 1 MC simulation        
            cs.nummed = 6; % nombre de couches dans l'image
            do_results_ui = 0;
        end           
            
        for iSess = 1:NSess
            arun=1; % already run               
            if do_results_ui
                if arun==0 %%% prevoir la sauvegarde des differentes matrices xSPM
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
                jobP.xSPM = xSPM;               

                % Name for current simulation
                cs_title = xSPM.title;
                cs.dir = fullfile(G.seg_p,csn);
                csn = [alg_nam '_' cs_title '_' daate];

            else                   
                % Name for current simulation
                csn = [alg_nam '_' daate];   
                cs.dir = fullfile(G.seg_p,csn);

            end
            % FOR EACH SIMULATION : WRITE CONFIGURATION FILES (.inp)
            % Directory for configuration files
            if ~exist(cs.dir,'dir')
                mkdir(cs.dir);
            end

            if cs.pve_cfg
                try % Resample the BOLD or user-specified mask to the space
                    % and resolution of the medium for the simulation

                    % User-specified mask (if given)
                    if isfield(job.pve_cfg,'pve_anat') && ~isempty(job.pve_cfg.pve_anat)
                        jobP.user_mask = job.pve_cfg.pve_anat{1};
                    end

                    jobP.T1seg = G.seg;%NIRS.Dt.ana.T1seg;
                    jobP.cs_dir = cs.dir;

                    outP = nirs_MCsegment_PVE(jobP);
                    G.seg = outP;


                catch
                    disp(['PVE failed for session ' iSess]);
                end
            end
            
            % Fin de config
            jobF.G = G;
            jobF.G.dir = cs.dir;
            jobF.G.voxelSize =cs.par.voxelSize;
            jobF.G.alg = cs.alg;
            jobF.P = P;
            outF = nirs_configMC_toMCframe(jobF);%images,positions and directions
            cs.segR = outF.G.segR;
            cs.b8i = outF.G.b8i;
            cs.n_b8i = outF.G.n_b8i;
            cs.ROIlimits = outF.G.ROIlimits;
            cs.P = outF.P;
            cs.Pvoid = Pvoid;
            cs.Pfp_rmiv = outF.P.Pfp_rmiv;
            cs.Pwd_rmm = outF.P.Pwd_rmm;

            Sr = cs.par.radiis * ones(cs.NSkpt,1);
            Dr = cs.par.radiid * ones(cs.NDkpt,1);
            cs.P.r = [Sr' Dr' zeros(1,size(cs.Pkpt,1) -(cs.NSkpt+cs.NDkpt))];

            NIRS.Cs.mcs{i_cs} = cs;
            NIRS.Cs.n{i_cs} = csn;

            jobW = cs;
            jobW.wl_dev = NIRS.Cf.dev.wl;

            % FOR EACH SIMULATION : WRITE CONFIGURATION FILES (.inp)
            % Directory for configuration files
            if ~exist(cs.dir,'dir')
                mkdir(cs.dir);
            end
            % Write the files
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
