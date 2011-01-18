function out = nirs_run_configMC(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: Generate .bin segmentation volume (Reduce the .nii segmented volume
%with int16 precisions on intensities to a .nii volume with uint8 precision)
%Then generate .cfg or .inp configuration files

%Note to clarify: using index nothing for anisotropic, not for interpolated, and iso
%for isotropic

% Ce que ca fait :
% 1 : on travaille dans le domaine des voxels
% En fait on reste dans l'espace rigide originel de l'IRM. La seule chose
% est que l'on passe des voxels parallépipédiques aux voxels
% isotropiques mais en gardant les mêmes orientations !!!
% 2 : changer de uint16 a uint8

load(job.NIRSmat{1,1});
NIRS.Cs = {};

if ~exist(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]),'dir')%Directory for configuration files
    mkdir(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]));
end

% on test pour voir co;bien il esxiste de simulation deja roulee :
if isfield(NIRS.Cs,'mcs'), cs_n = size(NIRS.Cs.mcs,1)+1; else cs_n=1; end
cs.alg = job.MC_CUDAchoice;%1=MCX ; 2=tMCimg
cs.dir = [job.MC_configdir job.MC_nam];
cs.par = job.MC_parameters;
cs.seg = job.image_in{1,1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATTENTION ON DOIT VERROUILLER POUR NE PAS CREER DEUX SIMS AVEC LE MEME NOM
NIRS.Cs.mcs{cs_n,1} = cs;
NIRS.Cs.n{cs_n,1} = job.MC_nam;
save(job.NIRSmat{1,1},'NIRS');

jobR.image_in = job.image_in;
jobR.out_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]);
jobR.out_dim = [1 1 1];
jobR.out_dt = 'same';
jobR.out_autonaming = 0;
jobR.out_prefix = 'prefix';
outR =nirs_resize(jobR);

clear NIRS

V_rmiv = spm_vol(outR);
Y_rmiv = spm_read_vols(V_rmiv);

if job.MC_CUDAchoice==1
    Y_rmiv=permute(Y_rmiv,[2,1,3]);
end
Y8_rmiv = uint8(Y_rmiv);% en fait on est dans des voxels de 1 mm d'ou le racourci

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENTION SI CAS 3 ALORS IL FAUT FAIRE DEUX CHEMINS CAR POUR UN IL FAUT
% INVERSER MAIS PAS POUR L AUTRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(job.NIRSmat{1,1});
cs = NIRS.Cs.mcs{cs_n,:};

V = spm_vol(cs.seg);
dim = V.dim;
inv_mat = spm_imatrix(V.mat);
scalings = diag(inv_mat(7:9));

dim_rmiv = ceil((dim-1) * abs(scalings));

NS = NIRS.Cf.H.S.N;
ND = NIRS.Cf.H.D.N;
NQ = NIRS.Cf.H.Q.N;

[~,id,~] = fileparts(V.fname);
if strcmp(id(1:3),'roi')
    n_id = ['roi_' id(21:end)];
else
    n_id = id(17:end);
end
n = [n_id '-vol8bit-' num2str(dim_rmiv(1)) 'x' num2str(dim_rmiv(2)) 'x' num2str(dim_rmiv(3)) '.bin'];

cs.segR = outR;
cs.b8i = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam],n);
NIRS.Cs.mcs{cs_n,:} = cs;
save(job.NIRSmat{1,1},'NIRS');

fid = fopen(fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam],n),'wb');
fwrite(fid, Y8_rmiv, 'uint8');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the .cfg (tMCimg) or .inp (CUDA MCX) configuration files

%Required data: optodes.positionsInROI.s, optodes.directions.s, radii.s and
%for sources, and same for detectors, with .d

% Source and detector positions
if isfield(NIRS.Cf.H.P,'void')
    Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
else
    Pvoid = {};
end
NP = NIRS.Cf.H.P.N;
Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;

%Transform MNI mm -> MNI voxels
Pfp_rmv = V.mat\[Pfp_rmm; ones(1,NP)];
%Transform MNI voxels -> MNI isotropic voxels
scalings4 = zeros(4);
scalings4(4,4)=1;
scalings4(1:3,1:3) = scalings;
Pfp_rmiv = scalings4 * Pfp_rmv;

%%% Directions
Pd_rmm = Pp_rmm - Pp_c1_rmm;

R = V.mat(1:3,1:3); %no translations
%Transform MNI mm -> MNI voxels
Pd_rmv = R\Pd_rmm;
%Transform Voxels -> Isotropic Voxels Space
Pd_rmiv = scalings * Pd_rmv;
%normalize directions
Pwd_rmiv = zeros(3,NP);

for iP=1:NP
    temp_dir = Pd_rmiv(:,iP);
    lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
    Pwd_rmiv(:,iP) = temp_dir/lgth;
end


Sr = job.MC_parameters.radiis * ones(NIRS.Cf.H.S.N,1);
Dr = job.MC_parameters.radiid * ones(NIRS.Cf.H.D.N,1);

parameters.nphotons = job.MC_parameters.nphotons;
parameters.seed = job.MC_parameters.seed;
parameters.modulationFreq = job.MC_parameters.modulationFreq;
parameters.voxelSize = job.MC_parameters.voxelSize;
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
    
    % tMCimg
    jobW.n_id = n_id;
    jobW.mc_dir = fullfile(NIRS.Dt.s.p,[job.MC_configdir job.MC_nam]);

    jobW.n = cs.b8i;
    
    jobW.dim_rmiv = dim_rmiv;
    jobW.ROIlimits = [1 1 1; V.dim];
    
    jobW.parameters = parameters;
    jobW.NS = NS;
    jobW.ND = ND;
    jobW.Pvoid = Pvoid;
    
    P.p = abs(round(Pfp_rmiv));
    P.wd = Pwd_rmiv;
    P.r = [Sr' Dr' zeros(1,NP-(NS+ND))];
    jobW.P =P;
    
    jobW.wl = NIRS.Cf.dev.wl(iwl);
    
    out = nirs_configMC_writeCFGfiles(jobW);
end