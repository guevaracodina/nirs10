function out = nirs_configMC_createCS(job)
% Prepare images, positions and directions !
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Clément Bonnéry 07/2011

cs = job.cs;

% transform image from anisotropic voxels space to isotropic voxels space
jobRS.image_in = {cs.seg};
jobRS.out_dir = fullfile(cs.p,cs.dir);
jobRS.out_dim = [1 1 1];
jobRS.out_vxsize = cs.par.voxelSize;
jobRS.out_dt = 'same';
jobRS.out_autonaming = 0;
jobRS.out_prefix = 'prefix';
outRS =nirs_resize(jobRS);

cs.segR = outRS;

V_rmiv = spm_vol(cs.segR);
Y_rmiv = spm_read_vols(V_rmiv);
Y8_rmiv = uint8(Y_rmiv);

V = spm_vol(cs.seg);
% dim = V.dim;% on est ici avec l'ancienne image qui permet d'obtenir les bonnes positions en voxels
inv_mat = spm_imatrix(V.mat);

% Transform also P positions and directions
NP = size(cs.Pfp_rmv,2);
% Positions : Transform MNI mm -> MNI isotropic voxels
for i=1:size(cs.Pfp_rmm,2)
    Pfp_ancienne_rmv(:,i) = V.mat\[cs.Pfp_rmm(:,i);1];
    Pfp_ancienne_rmiv(:,i) = abs(inv_mat(7:9)').*(Pfp_ancienne_rmv(1:3,i)/cs.par.voxelSize);
end

Pfp_ancienne_rmiv = round(Pfp_ancienne_rmiv);

% Directions
Pd_rmm = cs.Pp_rmm - cs.Pp_c1_rmm;
Pwd_rmm = zeros(3,NP);
for iP=1:NP
    temp_dir = Pd_rmm(:,iP);
    lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
    Pwd_rmm(:,iP) = temp_dir/lgth;
end

if cs.alg==1
    jobF.Pp_rmm = cs.Pp_rmm;
    jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
    jobF.NP = NP;
    jobF.image_in = {cs.segR};
    jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
    jobF.lby = 'configMC_MCX';
    outF = nirs_fit_probe(jobF);
    Pfp_ancienne_rmiv = outF{1};
elseif cs.alg==2 % pour tMC, les points doivent etre dans le volume (peut etre est-ce juste une question de direction de lq propagation...)
    jobF.Pp_rmm = cs.Pp_rmm;
    jobF.Pp_c1_rmm = cs.Pp_c1_rmm;
    jobF.NP = NP;
    jobF.image_in = {cs.segR};
    jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
    jobF.lby = 'configMC_tMC';
    outF = nirs_fit_probe(jobF);
    Pfp_ancienne_rmiv = outF{1};
end
cs.Pfp_rmiv = Pfp_ancienne_rmiv;
cs.Pwd_rmm = Pwd_rmm;

% 8bits .bin image
[dummy,id,dummy2] = fileparts(V_rmiv.fname);
dim_rmiv = V_rmiv.dim;

cs.n_b8i = ['vol8bit_' id '.bin'];
cs.b8i = fullfile(cs.p,cs.dir,cs.n_b8i);

fid = fopen(fullfile(cs.p,cs.dir,cs.n_b8i),'wb');
fwrite(fid, Y8_rmiv, 'uint8');
fclose(fid);

cs.numTimeGates = cs.par.numTimeGates;
cs.deltaT = cs.par.deltaT;

if ~cs.roi && isfield(NIRS.Cf.H.P,'void')
    Pvoid = NIRS.Cf.H.P.void;% Keep track of non-existent sources/detectors, to exclude them explicitly later
else
    Pvoid = {};
end

cs.Pvoid = Pvoid;

Sr = cs.par.radiis * ones(cs.NSkpt,1);
Dr = cs.par.radiid * ones(cs.NDkpt,1);


cs.ROIlimits = [1 1 1; dim_rmiv];

for i=1:size(cs.Pwd_rmm,2)
    Pwd_rmiv(:,i) = V_rmiv.mat(1:3,1:3)\cs.Pwd_rmm(:,i);
    %                 Pfp_tMC(:,i) = V_rmiv.mat*[Pfp_ancienne_rmiv(:,i);1];
end

if cs.alg==1
    %%% MCX en voxel /////!!!!!!!!!!!!!!!!!!!!!!!!!!
    % comme on travaille sur les donnees sur des voxels isotropiques de
    % 1mm, il suffit de diviser par la taille des voxels (dans le cas
    % ou on serait avec des voxels de taille custom, il faudrait
    % diviser la taille des voxels voulu par la taille courante)
    P.p =Pfp_ancienne_rmiv;%/(parameters.voxelSize); DEFINITIFFFFFF
    P.wd = -Pwd_rmiv(1:3,:); %to point toward the inside of brain
    
elseif cs.alg==2
    % MonteCarlo in a particular frame. Positions must be in mm but the
    % origin is the same as the origin of the voxel frame (these positions
    % don't respect SPM conventions) %%definitif ////
    P.p = parameters.voxelSize*Pfp_ancienne_rmiv;
    P.wd = -V_rmiv.mat(1:3,1:3)*Pwd_rmm;
end


G.alg = cs.alg;
G.cs_dir = cs.dir;
G.n_b8i = cs.n_b8i;
G.ROIlimits = cs.ROIlimits;
G.par =cs.par;
G.pve_cfg =cs.pve_cfg;

P.NSinit = cs.NSinit;
P.NS = cs.NSkpt;
P.ND = cs.NDkpt;
P.Pvoid = cs.Pvoid;
P.Pkpt = cs.Pkpt;

P.r = [Sr' Dr' zeros(1,NP -(cs.NSkpt+cs.NDkpt))];
P.Sr = cs.par.radiis;
P.Dr = cs.par.radiid;

out.G = G;
out.P = P;
end