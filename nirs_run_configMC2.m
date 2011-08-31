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
        
        %%% on definit un nom
        if sum(isfield(job.pve_cfg,{'pve_bold','pve_asl','pve_anat'})) %%%% Thresholded BOLD image is concidered as layer 6
            %%%% it might be also a perturbation included by user... any use ??
            % on cherche toutes les sessions qui pourraient mener a un
            
            % calcul de PVE
            last = size(NIRS.Dt.fir.pp,2);
            NSess = size(NIRS.Dt.fir.pp(1,last).p,1);
            
            for iSess =1:NSess
                arun=1;
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
                
                cs_title = xSPM.title;
                csn = [alg_nam '_' cs_title '_' daate];
                cs.dir =fullfile(G.seg_p,csn);
                if ~exist(cs.dir,'dir')%Directory for configuration files
                    mkdir(cs.dir);
                end
                
                try
                    jobP.xSPM = xSPM;
                    jobP.T1seg = G.seg;%NIRS.Dt.ana.T1seg;
                    jobP.cs_dir = cs.dir;
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
                    cs.n_b8i = outF.G.n_b8i;
                    cs.ROIlimits = outF.G.ROIlimits;
                    cs.P = outF.P;
                    cs.Pvoid = Pvoid;
                    cs.Pfp_rmiv = outF.P.Pfp_rmiv;
                    cs.Pwd_rmm = outF.P.Pwd_rmm;
                    cs.nummed = 12;
                    
                    Sr = cs.par.radiis * ones(cs.NSkpt,1);
                    Dr = cs.par.radiid * ones(cs.NDkpt,1);
                    cs.P.r = [Sr' Dr' zeros(1,size(cs.Pkpt,1) -(cs.NSkpt+cs.NDkpt))];
            
                    NIRS.Cs.mcs{i_cs} = cs;
                    NIRS.Cs.n{i_cs} = csn;
                    
                    jobW = cs;
                    jobW.wl_dev = NIRS.Cf.dev.wl;
                    nirs_configMC_writeCFGfiles2(jobW);
                catch
                    disp(['PVE failed for session ' iSess]);
                end
            end
        elseif isfield(job.pve_cfg,'no_pve')
            csn = [alg_nam '_' daate];
            cs.dir =fullfile(G.seg_p,csn);
            %Directory for configuration files
            if ~exist(cs.dir,'dir')
                mkdir(cs.dir);
            end
            cs.pve_cfg = 0;
            
            % fin de config et ecriture
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
            cs.nummed = 6;
            
            Sr = cs.par.radiis * ones(cs.NSkpt,1);
            Dr = cs.par.radiid * ones(cs.NDkpt,1);
            cs.P.r = [Sr' Dr' zeros(1,size(cs.Pkpt,1) -(cs.NSkpt+cs.NDkpt))];
            
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