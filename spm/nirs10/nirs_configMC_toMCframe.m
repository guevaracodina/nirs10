function out = nirs_configMC_toMCframe(job)
% Prepare images, positions and directions !
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Clément Bonnéry 07/2011

G = job.G;
P = job.P;

% transform image from anisotropic voxels space to isotropic voxels space

jobRS.image_in = {G.seg};
jobRS.out_dir = G.dir;
jobRS.out_dim = [1 1 1];
jobRS.out_vxsize = G.voxelSize;
jobRS.out_dt = 'same';
jobRS.out_autonaming = 0;
jobRS.out_prefix = 'prefix';
outRS =nirs_resize(jobRS);


V_rmiv = spm_vol(outRS);
Y_rmiv = spm_read_vols(V_rmiv);
Y8_rmiv = uint8(Y_rmiv);

G.segR = outRS;

%%%%% calcul des positions
V = spm_vol(G.seg);
% dim = V.dim;% on est ici avec l'ancienne image qui permet d'obtenir les bonnes positions en voxels
inv_mat = spm_imatrix(V.mat);

% Transform also P positions and directions
NP = size(P.Pfp_rmv,2);
% Positions : Transform MNI mm -> MNI isotropic voxels
for i=1:size(P.Pfp_rmm,2)
    Pfp_ancienne_rmv(:,i) = V.mat\[P.Pfp_rmm(:,i);1];
    Pfp_ancienne_rmiv(:,i) = abs(inv_mat(7:9)').*(Pfp_ancienne_rmv(1:3,i)/G.voxelSize);
end

%%%% MODIFICATION 29 07 2011 %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%
% Pfp_ancienne_rmiv = round(Pfp_ancienne_rmiv);
Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
%%%% MODIFICATION 29 07 2011 %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%

% Directions
Pd_rmm = P.Pp_rmm - P.Pp_c1_rmm;
Pwd_rmm = zeros(3,NP);
for iP=1:NP
    temp_dir = Pd_rmm(:,iP);
    lgth = (temp_dir(1)^2 + temp_dir(2)^2 + temp_dir(3)^2)^(1/2);
    Pwd_rmm(:,iP) = temp_dir/lgth;
end

if G.alg==1
    jobF.Pwd_rmm = Pwd_rmm;
    jobF.Pp_rmm = P.Pp_rmm;
    jobF.Pp_c1_rmm = P.Pp_c1_rmm;
    jobF.NP = NP;
    jobF.image_in = {outRS};
    jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
    jobF.lby = 'configMC_MCX';
    outF = nirs_fit_probe(jobF);
    Pfp_ancienne_rmiv = outF{1};
elseif G.alg==2 % pour tMC, les points doivent etre dans le volume (peut etre est-ce juste une question de direction de la propagation...)
    jobF.Pwd_rmm = Pwd_rmm;
    jobF.Pp_rmm = P.Pp_rmm;
    jobF.Pp_c1_rmm = P.Pp_c1_rmm;
    jobF.NP = NP;
    jobF.image_in = {outRS};
    jobF.Pfp_ancienne_rmiv = Pfp_ancienne_rmiv;
    jobF.lby = 'configMC_tMC';
    outF = nirs_fit_probe(jobF);
    Pfp_ancienne_rmiv = outF{1};
end

% 8bits .bin image
[dummy,id,dummy2] = fileparts(V_rmiv.fname);
dim_rmiv = V_rmiv.dim;
G.ROIlimits = [1 1 1; dim_rmiv];

G.n_b8i = ['vol8bit_' id '.bin'];
G.b8i = fullfile(G.dir,G.n_b8i);

fid = fopen(G.b8i,'wb');
fwrite(fid, Y8_rmiv, 'uint8');
fclose(fid);


for i=1:size(Pwd_rmm,2)
    Pwd_rmiv(:,i) = V_rmiv.mat(1:3,1:3)\Pwd_rmm(:,i);
    %                 Pfp_tMC(:,i) = V_rmiv.mat*[Pfp_ancienne_rmiv(:,i);1];
end

if G.alg==1
    %%% MCX en voxel /////!!!!!!!!!!!!!!!!!!!!!!!!!!
    % comme on travaille sur les donnees sur des voxels isotropiques de
    % 1mm, il suffit de diviser par la taille des voxels (dans le cas
    % ou on serait avec des voxels de taille custom, il faudrait
    % diviser la taille des voxels voulu par la taille courante)
    P.p =Pfp_ancienne_rmiv;
    P.wd = -Pwd_rmiv(1:3,:); %to point toward the inside of brain
    
elseif G.alg==2
    % MonteCarlo in a particular frame. Positions must be in mm but the
    % origin is the same as the origin of the voxel frame (these positions
    % don't respect SPM conventions) %%definitif ////
    P.p = G.voxelSize*Pfp_ancienne_rmiv;
    P.wd = -V_rmiv.mat(1:3,1:3)*Pwd_rmm;% towards inside
end

P.Pfp_rmiv = Pfp_ancienne_rmiv;
P.Pwd_rmm = Pwd_rmm;

out.G = G;
out.P = P;
end