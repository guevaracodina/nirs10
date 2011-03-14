function out = nirs_run_configMC(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%______________________________________________________________________

%Usage: Generate .bin segmentation volume (Reduce the .nii segmented volume
%with int16 precisions on intensities to a .nii volume with uint8 precision)
%Then generate .cfg or .inp configuration files

%Note to clarify: using index nothing for anisotropic, not for interpolated, and iso
%for isotropic

% Ce que ca fait :
% 1 : on travaille dans le domaine des voxels
% En fait on reste dans l'espace rigide originel de l'IRM. La seule chose
% est que l'on passe des voxels parall�pip�diques aux voxels
% isotropiques mais en gardant les m�mes orientations !!!
% 2 : changer de uint16 a uint8

load(job.NIRSmat{1,1});

if isfield(NIRS.Cs,'mcs')
    i_cs = size(NIRS.Cs.mcs,2)+1;
    if ~sum(strcmp(NIRS.Cs.n,job.MC_nam))
        csn = job.MC_nam;
    else
        csn = [job.MC_nam '1'];
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

if isfield(job.mcim_cfg,'mcim_in')
    cs.seg = job.image_in{1,1};
    cs.Pn = NIRS.Cf.H.P.n;
else
    roi =1;
    cs.seg = NIRS.Cs.temp.segR; %ROI from temp
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

% transform image from anisotropic voxels space to isotropic voxels space
jobRS.image_in = {cs.seg};
jobRS.out_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]);
jobRS.out_dim = [1 1 1];
jobRS.out_dt = 'same';
jobRS.out_autonaming = 0;
jobRS.out_prefix = 'prefix';
outRS =nirs_resize(jobRS);

clear NIRS

V_rmiv = spm_vol(outRS);
Y_rmiv = spm_read_vols(V_rmiv);

if job.MC_CUDAchoice==1
    Y_rmiv=permute(Y_rmiv,[2,1,3]);
elseif job.MC_CUDAchoice==3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IL FAUT FAIRE DEUX CHEMINS CAR POUR UN IL FAUT INVERSER MAIS PAS POUR L AUTRE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Y8_rmiv = uint8(Y_rmiv);% en fait on est dans des voxels de 1 mm d'ou le racourci

load(job.NIRSmat{1,1});
cs = NIRS.Cs.mcs{i_cs};

V = spm_vol(cs.seg);
dim = V.dim;
inv_mat = spm_imatrix(V.mat);
scalings = diag(inv_mat(7:9));

% Transform also P positions and directions
NP = size(cs.Pfp_rmv,2);
% %Transform MNI voxels -> MNI isotropic voxels
% Pfp_rmiv = scalings * cs.Pfp_rmv;

% Positions
for i=1:size(cs.Pfp_rmm,2)
Pfp_ancienne_rmv(:,i) = V.mat\[cs.Pfp_rmm(:,i);1];
Pfp_ancienne_rmiv(:,i) = abs(inv_mat(7:9)').*Pfp_ancienne_rmv(1:3,i);
end

Pfp_ancienne_rmiv = round(Pfp_ancienne_rmiv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% on verifie qu'on n'a pas de pb apres le resizing
jobF.Pp_rmm = cs.Pp_rmm;
jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
jobF.NP = NP;
jobF.image_in = {outRS};
jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
jobF.lby = 'configMC';
outF = nirs_fit_probe(jobF);
Pfp_ancienne_rmiv = outF{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directions
Pd_rmm = cs.Pp_rmm - cs.Pp_c1_rmm;
%%%%test
Pwd_rmm = zeros(3,NP);
for iP=1:NP
    temp_dir = Pd_rmm(:,iP);
    lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
    Pwd_rmm(:,iP) = temp_dir/lgth;
end

% on prepare la suvegarde de cs
cs.Pfp_rmiv = Pfp_ancienne_rmiv;
cs.Pwd_rmm = Pwd_rmm;
% %%%%
% R = V.mat(1:3,1:3); %no translations
% %Transform MNI mm -> MNI voxels
% Pd_rmv = R\Pd_rmm;
% %Transform Voxels -> Isotropic Voxels Space
% Pd_rmiv = scalings * Pd_rmv;
% %normalize directions
% Pwd_rmiv = zeros(3,NP);
% for iP=1:NP
%     temp_dir = Pd_rmiv(:,iP);
%     lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
%     Pwd_rmiv(:,iP) = temp_dir/lgth;
% end

% 8bits .bin image
dim_rmiv = ceil((dim-1) * abs(scalings));

[~,id,~] = fileparts(V.fname);
n = ['vol8bit' id '-' num2str(dim_rmiv(1)) 'x' num2str(dim_rmiv(2)) 'x' num2str(dim_rmiv(3)) '.bin'];

cs.segR = outRS;
cs.b8i = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn],n);
NIRS.Cs.mcs{i_cs} = cs;
save(job.NIRSmat{1,1},'NIRS');

fid = fopen(fullfile(NIRS.Dt.s.p,[job.MC_configdir csn],n),'wb');
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
parameters.voxelSize = job.MC_parameters.voxelSize;
parameters.numTimeGates = job.MC_parameters.numTimeGates;
parameters.deltaT = job.MC_parameters.deltaT;

% Create a .cfg or .inp file for each optode and each wavelength
for iwl = 1:size(NIRS.Cf.dev.wl,2)
    %%%%%%%%%%%% il faut v�rifier les valeurs des longueurs d'ondes
    %%%%%%%%%%%% utilis�es dans le code...
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
    
    jobW.algo = job.MC_CUDAchoice;
    jobW.n_id = n;
    jobW.mc_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir csn]);
    
    jobW.n = cs.b8i;
    
    jobW.dim_rmiv = dim_rmiv;
    jobW.ROIlimits = [1 1 1; V.dim];
    
    jobW.parameters = parameters;
    jobW.NS = cs.NSkpt;
    jobW.ND = cs.NDkpt;
    jobW.Pvoid = Pvoid;
    
    P.p = Pfp_ancienne_rmiv;%cs.Pfp_rmm;%abs(round(Pfp_rmiv));
    P.wd = Pwd_rmm;%Pwd_rmiv;
    P.r = [Sr' Dr' zeros(1,NP -(cs.NSkpt+cs.NDkpt))];
    jobW.P =P;
    jobW.wl = NIRS.Cf.dev.wl(iwl);
    
    out = nirs_configMC_writeCFGfiles(jobW);
end

%
% function out = nirs_run_configMC(job)
% %_______________________________________________________________________
% % Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
% %                    �cole Polytechnique de Montr�al
% %______________________________________________________________________
%
% %Usage: Generate .bin segmentation volume (Reduce the .nii segmented volume
% %with int16 precisions on intensities to a .nii volume with uint8 precision)
% %Then generate .cfg or .inp configuration files
%
% %Note to clarify: using index nothing for anisotropic, not for interpolated, and iso
% %for isotropic
%
% % Ce que ca fait :
% % 1 : on travaille dans le domaine des voxels
% % En fait on reste dans l'espace rigide originel de l'IRM. La seule chose
% % est que l'on passe des voxels parall�pip�diques aux voxels
% % isotropiques mais en gardant les m�mes orientations !!!
% % 2 : changer de uint16 a uint8
%
% % try
% %     NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
% %     NewDirCopyNIRS = 1;
% % catch
% %     NewDirCopyNIRS = 0;
% % end
%
% load(job.NIRSmat{1,1});
%
% if ~exist(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]),'dir')%Directory for configuration files
%     mkdir(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]));
% end
%
% % on test pour voir combien il existe de simulation deja roulee :
% if isfield(NIRS.Cs,'mcs'), cs_n = size(NIRS.Cs.mcs,2)+1; else cs_n=1; end
% cs.alg = job.MC_CUDAchoice;%1=MCX ; 2=tMCimg
% cs.dir = [job.MC_configdir job.MC_nam];
% cs.par = job.MC_parameters;
%
% %%%%%%%%%%%%%%%%%%% bOITEUX
% if isfield(job.mcim_cfg,'mcim_in')
%     cs.seg = job.image_in{1,1};
%     cs.Pn = NIRS.Cf.H.P.n;
% else
%     cs.seg = NIRS.Cs.temp.segR;
%     cs.Pfp_rmv = NIRS.Cs.temp.Pfp_roi_rmv;
%     cs.Pp_rmm = NIRS.Cs.temp.Pp_roi_rmm;
%     cs.Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm;
%     cs.Pkpt = NIRS.Cs.temp.Pkpt;
%     cs.NSkpt = NIRS.Cs.temp.NDkpt;
%     cs.NDkpt = NIRS.Cs.temp.NSkpt;
% end
%
% NIRS.Cs.mcs{cs_n} = cs;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ATTENTION ON DOIT VERROUILLER POUR NE PAS CREER DEUX SIMS AVEC LE MEME NOM
% NIRS.Cs.n{cs_n,1} = job.MC_nam;
% save(job.NIRSmat{1,1},'NIRS');
%
% jobR.image_in = {cs.seg};
% jobR.out_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]);
% jobR.out_dim = [1 1 1];
% jobR.out_dt = 'same';
% jobR.out_autonaming = 0;
% jobR.out_prefix = 'prefix';
% outR =nirs_resize(jobR);
%
% clear NIRS
%
% V_rmiv = spm_vol(outR);
% Y_rmiv = spm_read_vols(V_rmiv);
%
% if job.MC_CUDAchoice==1
%     Y_rmiv=permute(Y_rmiv,[2,1,3]);
% elseif job.MC_CUDAchoice==3
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % IL FAUT FAIRE DEUX CHEMINS CAR POUR UN IL FAUT INVERSER MAIS PAS POUR L AUTRE
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% Y8_rmiv = uint8(Y_rmiv);% en fait on est dans des voxels de 1 mm d'ou le racourci
%
% load(job.NIRSmat{1,1});
% cs = NIRS.Cs.mcs{cs_n};
%
% V = spm_vol(cs.seg);
% dim = V.dim;
% inv_mat = spm_imatrix(V.mat);
% scalings = diag(inv_mat(7:9));
%
% dim_rmiv = ceil((dim-1) * abs(scalings));
%
% [~,id,~] = fileparts(V.fname);
% if strcmp(id(1:3),'roi')
%     n_id = ['roi_' id(21:end)];
% else
%     n_id = id(17:end);
% end
% n = [n_id '-vol8bit-' num2str(dim_rmiv(1)) 'x' num2str(dim_rmiv(2)) 'x' num2str(dim_rmiv(3)) '.bin'];
%
% cs.segR = outR;
% cs.b8i = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam],n);
% NIRS.Cs.mcs{cs_n} = cs;
% save(job.NIRSmat{1,1},'NIRS');
%
% fid = fopen(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam],n),'wb');
% fwrite(fid, Y8_rmiv, 'uint8');
% fclose(fid);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Generate the .cfg (tMCimg) or .inp (CUDA MCX) configuration files
%
% %Required data: optodes.positionsInROI.s, optodes.directions.s, radii.s and
% %for sources, and same for detectors, with .d
%
% % Source and detector positions
% % Pn = cs.Pn;
% % NS = NIRS.Cf.H.S.N;
% % ND = NIRS.Cf.H.D.N;
% % Cid = NIRS.Cf.H.C.id;
% % NQ = NIRS.Cf.H.Q.N;
% % Cre =[];
% % %on extrait les sources
% % Sn_roi = Pn(Pn<=NS);
% % %on reconstruit les paires
% % CidS = Cid(2,:);
% % for i=1:length(Sn_roi)
% %     Cpossible{i} = Cid(1,CidS==Sn_roi(i));
% %     Cposs = Cpossible{i};
% %     for j=1:length(Cposs)
% %         if sum(Pn==Cid(3,Cposs(j)))
% %             Cre = [Cre Cposs(j)];
% %         end
% %     end
% % end
% %
% % if Cre==[]
% %     disp('aucune paire');
% % end
%
% %Transform MNI mm -> MNI voxels
% % % Pfp_rmv = V.mat\[Pfp_rmm; ones(1,NP)];
% Pfp_rmv = cs.Pfp_rmv;
% Pp_rmm = cs.Pp_rmm;
% Pp_c1_rmm = cs.Pp_c1_rmm;
% %%% attentionm on red2finit tous les pqrq;etres des P
% NP = size(Pfp_rmv,2);
%
% if isfield(NIRS.Cf.H.P,'void')
%     Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
% else
%     Pvoid = {};
% end
%
% for i=1:NIRS.Cf.H.P.N
%     if ~sum(cs.Pkpt==i)
%         Pvoid(i) = 1;
%     end
% end
%
% %Transform MNI voxels -> MNI isotropic voxels
% scalings4 = zeros(4);
% scalings4(4,4)=1;
% scalings4(1:3,1:3) = scalings;
% Pfp_rmiv = scalings4 * [Pfp_rmv;zeros(1,NP)];
%
% %%% Directions
% Pd_rmm = Pp_rmm - Pp_c1_rmm;
%
% R = V.mat(1:3,1:3); %no translations
% %Transform MNI mm -> MNI voxels
% Pd_rmv = R\Pd_rmm;
% %Transform Voxels -> Isotropic Voxels Space
% Pd_rmiv = scalings * Pd_rmv;
% %normalize directions
% Pwd_rmiv = zeros(3,NP);
%
% for iP=1:NP
%     temp_dir = Pd_rmiv(:,iP);
%     lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
%     Pwd_rmiv(:,iP) = temp_dir/lgth;
% end
%
% Sr = job.MC_parameters.radiis * ones(cs.NSkpt,1);
% Dr = job.MC_parameters.radiid * ones(NIRS.Cf.H.D.N,1);
%
% parameters.nphotons = job.MC_parameters.nphotons;
% parameters.seed = job.MC_parameters.seed;
% parameters.modulationFreq = job.MC_parameters.modulationFreq;
% parameters.voxelSize = job.MC_parameters.voxelSize;
% parameters.numTimeGates = job.MC_parameters.numTimeGates;
% parameters.deltaT = job.MC_parameters.deltaT;
%
% % Create a .cfg or .inp file for each optode and each wavelength
% for iwl = 1:size(NIRS.Cf.dev.wl,2)
%     %%%%%%%%%%%% il faut v�rifier les valeurs des longueurs d'ondes
%     %%%%%%%%%%%% utilis�es dans le code...
%     if iwl == 2 % 830 pour le code attend que 830 soit la premiere longueur d'onde JE CROIS
%         parameters.gmPpties = job.MC_parameters.gmPpties_l1;
%         parameters.wmPpties = job.MC_parameters.wmPpties_l1;
%         parameters.csfPpties = job.MC_parameters.csfPpties_l1;
%         parameters.skullPpties = job.MC_parameters.skullPpties_l1;
%         parameters.scalpPpties = job.MC_parameters.scalpPpties_l1;
%         parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l1+job.MC_parameters.gmPpties_l1;
%     else %set them all to _l2 even if more than two wavelengths
%         parameters.gmPpties = job.MC_parameters.gmPpties_l2;
%         parameters.wmPpties = job.MC_parameters.wmPpties_l2;
%         parameters.csfPpties = job.MC_parameters.csfPpties_l2;
%         parameters.skullPpties = job.MC_parameters.skullPpties_l2;
%         parameters.scalpPpties = job.MC_parameters.scalpPpties_l2;
%         parameters.perturbationPpties = job.MC_parameters.perturbationPpties_l2+job.MC_parameters.gmPpties_l2;
%     end
%
%     % tMCimg
%     jobW.n_id = n_id;
%     jobW.mc_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]);
%
%     jobW.n = cs.b8i;
%
%     jobW.dim_rmiv = dim_rmiv;
%     jobW.ROIlimits = [1 1 1; V.dim];
%
%     jobW.parameters = parameters;
%     jobW.NS = cs.NSkpt;
%     jobW.ND = cs.NDkpt;
%     jobW.Pvoid = Pvoid;
%
%     P.p = abs(round(Pfp_rmiv));
%     P.wd = Pwd_rmiv;
%     P.r = [Sr' Dr' zeros(1,NIRS.Cf.H.P.N-(cs.NSkpt+cs.NDkpt))];
%     jobW.P =P;
%
%     jobW.wl = NIRS.Cf.dev.wl(iwl);
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %&% en fin de fichier
%     %         if NewDirCopyNIRS
%     %             newNIRSlocation = fullfile(dir2,'NIRS.mat');
%     %             save(newNIRSlocation,'NIRS');
%     %             job.NIRSmat{Idx,1} = newNIRSlocation;
%     %         else
%     %             save(job.NIRSmat{Idx,1},'NIRS');
%     %         end
%
%
%     out = nirs_configMC_writeCFGfiles(jobW);
% end