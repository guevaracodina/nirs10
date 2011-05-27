function out = nirs_run_configMC(job)
% Manages the writing of the configurqtion files for either MCX or tMCimg
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
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clement Bonnery
%Overwrite to fix path
job.MC_configdir = 'MC';

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});

        if isfield(NIRS.Cs,'mcs')
            i_cs = size(NIRS.Cs.mcs,2)+1;
            if ~sum(strcmp(NIRS.Cs.n,job.MC_nam))
                csn = job.MC_nam;
            else
                csn = [job.MC_nam strrep(datestr(now),':','-')];
            end
        else
            i_cs=1;
            csn = job.MC_nam;
        end

        if ~exist(fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]),'dir')%Directory for configuration files
            mkdir(fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]));
        end

        % cs current simulation
        cs.alg = job.MC_CUDAchoice;%1=MCX ; 2=tMCimg
        cs.dir = [job.MC_configdir csn];
        cs.par = job.MC_parameters;

        if isfield(job.mcim_cfg,'mcim_in')% image segmentee de l'anatomique de base
            roi =0;% image choisie
            cs.seg = job.mcim_cfg.mcim_in{:};
            cs.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
            cs.Pfp_rmm = NIRS.Cs.temp.Pfp_roi_rmm;
            cs.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
            cs.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
            cs.Pkpt = NIRS.Cs.temp.Pkpt;
            cs.NSkpt = NIRS.Cs.temp.NSkpt;
            cs.NDkpt = NIRS.Cs.temp.NDkpt;
        else
            roi =1;
            try cs.seg = NIRS.Cs.temp.segR{:}; %ROI from temp
            catch
                cs.seg = NIRS.Cs.temp.segR;
            end
            cs.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
            cs.Pfp_rmm = NIRS.Cs.temp.Pfp_roi_rmm;
            cs.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
            cs.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
            cs.Pkpt = NIRS.Cs.temp.Pkpt;
            cs.NSkpt = NIRS.Cs.temp.NSkpt;
            cs.NDkpt = NIRS.Cs.temp.NDkpt;
        end

        NIRS.Cs.mcs{i_cs} = cs;
        NIRS.Cs.n{i_cs} = csn;
        save(job.NIRSmat{1,1},'NIRS');

        parameters.voxelSize = job.MC_parameters.voxelSize;
        % transform image from anisotropic voxels space to isotropic voxels space
        jobRS.image_in = {cs.seg};
        jobRS.out_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]);
        jobRS.out_dim = [1 1 1];
        jobRS.out_vxsize = parameters.voxelSize;
        jobRS.out_dt = 'same';
        jobRS.out_autonaming = 0;
        jobRS.out_prefix = 'prefix';
        outRS =nirs_resize(jobRS);

        clear NIRS

        V_rmiv = spm_vol(outRS);
        Y_rmiv = spm_read_vols(V_rmiv);

        % % % % % % % % % % % % % % % % %%%%%% attention on N'inverse JAMAIS !!!!!!!!!!!!!!!!!!!!!!!!!
        % % % % % % % % % % % % % % % % if job.MC_CUDAchoice==1
        % % % % % % % % % % % % % % % % %     Y_rmiv=permute(Y_rmiv,[2,1,3]);
        % % % % % % % % % % % % % % % % elseif job.MC_CUDAchoice==3
        % % % % % % % % % % % % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % % % % % % % % % % %     % IL FAUT FAIRE DEUX CHEMINS CAR POUR UN IL FAUT INVERSER MAIS PAS POUR L AUTRE
        % % % % % % % % % % % % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % % % % % % % % % % % end
        Y8_rmiv = uint8(Y_rmiv);

        load(job.NIRSmat{1,1});
        cs = NIRS.Cs.mcs{i_cs};

        V = spm_vol(cs.seg);
        % dim = V.dim;% on est ici avec l'ancienne image qui permet d'obtenir les bonnes positions en voxels
        inv_mat = spm_imatrix(V.mat);
        % scalings = diag(inv_mat(7:9));

        % Transform also P positions and directions
        NP = size(cs.Pfp_rmv,2);
        % Positions : Transform MNI mm -> MNI isotropic voxels
        for i=1:size(cs.Pfp_rmm,2)
        Pfp_ancienne_rmv(:,i) = V.mat\[cs.Pfp_rmm(:,i);1];
        Pfp_ancienne_rmiv(:,i) = abs(inv_mat(7:9)').*(Pfp_ancienne_rmv(1:3,i)/parameters.voxelSize);
        end

        Pfp_ancienne_rmiv = round(Pfp_ancienne_rmiv);

        %%%%%%% appliquer ca sur le 8bit !!!!!!!!!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if job.MC_CUDAchoice==1 % MCX : sources et detecteurs doivent etre DANS
        % le volume
        %%%%%%%%%%%%%%%%%%%%
        % pour MCX, ca n'a pas d'effet... peut etre on projette pas sur la bonne
        % image!!!!!
        %%%%%%%%%%%%%%%%%%%%
        %%%%%%%pour tMCimg : je sais pas trop...
            jobF.Pp_rmm = cs.Pp_rmm;
            jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
            jobF.NP = NP;
            jobF.image_in = {outRS};
            jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
            jobF.lby = 'configMC';
            outF = nirs_fit_probe(jobF);
            Pfp_ancienne_rmiv = outF{1};
        % end

        % Directions
        Pd_rmm = cs.Pp_rmm - cs.Pp_c1_rmm;
        Pwd_rmm = zeros(3,NP);
        for iP=1:NP
            temp_dir = Pd_rmm(:,iP);
            lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
            Pwd_rmm(:,iP) = temp_dir/lgth;
        end

        % on prepare la sauvegarde de cs
        cs.Pfp_rmiv = Pfp_ancienne_rmiv;
        cs.Pwd_rmm = Pwd_rmm;

        % 8bits .bin image
        [~,id,~] = fileparts(V_rmiv.fname);
        dim_rmiv = V_rmiv.dim;
        n_b8i = ['vol8bit_' id '.bin'];

        cs.segR = outRS;
        cs.b8i = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn],n_b8i);

        cs.numTimeGates = job.MC_parameters.numTimeGates;
        cs.deltaT = job.MC_parameters.deltaT;
        %Path for configuration files
        NIRS.Cs.mcs{i_cs} = cs;
        save(job.NIRSmat{1,1},'NIRS');

        fid = fopen(fullfile(NIRS.Dt.s.p,[job.MC_configdir csn],n_b8i),'wb');
        fwrite(fid, Y8_rmiv, 'uint8');
        fclose(fid);


        %Generate the .cfg (tMCimg) or .inp (CUDA MCX) configuration files
        %Required data: optodes.positionsInROI.s, optodes.directions.s, radii.s and
        %for sources, and same for detectors, with .d
        if ~roi && isfield(NIRS.Cf.H.P,'void')
            Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
        else
            Pvoid = {};
        end

        Sr = job.MC_parameters.radiis * ones(cs.NSkpt,1);
        Dr = job.MC_parameters.radiid * ones(cs.NDkpt,1);

        parameters.nphotons = job.MC_parameters.nphotons;
        parameters.seed = job.MC_parameters.seed;
        parameters.modulationFreq = job.MC_parameters.modulationFreq;

        parameters.numTimeGates = job.MC_parameters.numTimeGates;
        parameters.deltaT = job.MC_parameters.deltaT;

        % Create a .cfg or .inp file for each optode and each wavelength
        for iwl = 1:size(NIRS.Cf.dev.wl,2)
            %%%%%%%%%%%% il faut vérifier les valeurs des longueurs d'ondes
            %%%%%%%%%%%% utilisées dans le code...
            if iwl == 2 % 830 pour le code attend que 830 soit la premiere longueur d'onde JE CROIS
                parameters.gmPpties = job.MC_parameters.gmPpties_l1;
                parameters.wmPpties = job.MC_parameters.wmPpties_l1;
                parameters.csfPpties = job.MC_parameters.csfPpties_l1;
                parameters.skullPpties = job.MC_parameters.skullPpties_l1;
                parameters.scalpPpties = job.MC_parameters.scalpPpties_l1;
                parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l1+job.MC_parameters.gmPpties_l1;
            else %set them all to _l2 even if more than two wavelengths
                parameters.gmPpties = job.MC_parameters.gmPpties_l2;
                parameters.wmPpties = job.MC_parameters.wmPpties_l2;
                parameters.csfPpties = job.MC_parameters.csfPpties_l2;
                parameters.skullPpties = job.MC_parameters.skullPpties_l2;
                parameters.scalpPpties = job.MC_parameters.scalpPpties_l2;
                parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l2+job.MC_parameters.gmPpties_l2;
            end

            %%%% attention la on est avec l'image resizee...........
            jobW.algo = job.MC_CUDAchoice;
            jobW.n_b8i = n_b8i;
            jobW.mc_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]);    
            jobW.ROIlimits = [1 1 1; dim_rmiv];

            jobW.parameters = parameters;
            jobW.NS = cs.NSkpt;
            jobW.ND = cs.NDkpt;
            jobW.Pvoid = Pvoid;

            if job.MC_CUDAchoice==1
                %%% MCX en voxel /////!!!!!!!!!!!!!!!!!!!!!!!!!!
                % comme on travaille sur les donnees sur des voxels isotropiques de
                % 1mm, il suffit de diviser par la taille des voxels (dans le cas 
                % ou on serait avec des voxels de taille custom, il faudrait 
                % diviser la taille des voxels voulu par la taille courante)
                P.p = Pfp_ancienne_rmiv;%/(parameters.voxelSize);
            elseif job.MC_CUDAchoice==2
                % MonteCarlo in a particular frame. Positions must be in mm but the
            % origin is the same as the origin of the voxel frame (these positions 
            % don't respect SPM conventions)
            P.p = parameters.voxelSize*Pfp_ancienne_rmiv;
            end
            
            for i=1:size(cs.Pwd_rmm,2)
                Pwd_rmiv(:,i) = V_rmiv.mat(1:3,1:3)\cs.Pwd_rmm(:,i);           
            end
            P.wd = -Pwd_rmiv(1:3,:); %to point toward the inside of brain
            P.r = [Sr' Dr' zeros(1,NP -(cs.NSkpt+cs.NDkpt))];
            jobW.P =P;
            jobW.wl = NIRS.Cf.dev.wl(iwl);

            out_void = nirs_configMC_writeCFGfiles(jobW);
        end
        newNIRSlocation = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn],'NIRS.mat');
        save(newNIRSlocation,'NIRS');
        job.NIRSmat{Idx,1} = newNIRSlocation;
    catch exception
        disp(exception.identifier);
        disp(['Could not run MonteCarlo configuration for subject' int2str(Idx)]);
    end  
end
out.NIRSmat = job.NIRSmat;