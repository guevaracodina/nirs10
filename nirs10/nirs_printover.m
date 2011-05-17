function out = nirs_printover(job)
% - Tu as une image fonctionnelle dans un référentiel R_fonc et une image 
% anatomique dans un référentiel R_ana, je pense que les deux sont dans MNI
% et en mm.
% - Tu travailles en mm (hyper important à mon avis) !!! Tu utilises la 
% normalisation pour être dans les mêmes espaces pour les deux images, le 
% code te produit alors la transformation qui permet de faire coller au 
% mieux les deux images (là c'est le code de coregistration que tu utilises).
% - Tu devrais à ce point savoir comment passer d'une image à l'autre (en mm).
% Il te suffit alors d'utiliser le changement de tailles des voxels pour 
% faire correspondre les tailles des voxels dans chacune des images 
% (ré échantillonnage et interpolation de l'image sur le nouvel espace).

%% cas du résultat d'une simulation MC qu'on veut superposer à l'anatomique
image_in = job.image_in; % pour l'instant on prend l'image segmentee
outMC = job.outMCfile;
% anatomique sous forme de NIFTI qui a servi de base pour les fichiers cfg
% de la simulation !!
% résultat simulation (matrice 3D) dans l'espace des voxels isotropiques


%%%% fonction de Frederic : nirs_read_2pt(file,nx,ny,nz,nt)
% Read and sum 2pt files
% Any other initializations would go here

fid = fopen(outMC, 'rb');
Yf=zeros(nx*ny*nz,1);

% nt = points temporels
for index=1:nt
    Yf=Yf+fread(fid,nx*ny*nz,'double');
end
fclose(fid);

Yf=reshape(Yf,[ny nx nz]);
% Some values are negative... should know why
Yf=(Yf>0).*Yf;

% on revient dans l'espace non-isotropique :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXTRAIT DE run_configMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inv_mat = spm_imatrix(V.mat); %parameters of inverse transform
scalings = diag(inv_mat(7:9)); 

% Passage de isovoxel a voxel

%Define dimension of new volume
%New volume has size 1mm X 1mm X 1mm
dim_rmiv = ceil((dim-1) * abs(scalings));
%Create linear spaces for initial volume (voxel) data and final (in mm) [en fait voxels d'un mm]
x_mm2v = linspace(0,dim_rmiv(1)-1,dim(1));
y_mm2v = linspace(0,dim_rmiv(2)-1,dim(2));
z_mm2v = linspace(0,dim_rmiv(3)-1,dim(3));
x_rmiv = 0:1:dim_rmiv(1)-1;
y_rmiv = 0:1:dim_rmiv(2)-1;
z_rmiv = 0:1:dim_rmiv(3)-1;
JE pense au4il fqut pqsser pqr les mm pour qvoir des points d4interpolations


%Make grids, to use with 3-dimensional interpolation
[xg_rmiv,yg_rmiv,zg_rmiv] = meshgrid(x_rmiv,y_rmiv,z_rmiv);
[xg_rmv,yg_rmv,zg_rmv] = meshgrid(x_rmv,y_rmv,z_rmv);

Y=permute(Y,[2,1,3]);
Y_rmiv = interp3(xg_rmiv,yg_rmiv,zg_rmiv,Y,xg_rmv,yg_rmv,zg_rmv,'nearest');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mise sous forme de nifti :
Vf.dim = [ny nx nz];

V = read_vol(image_in);

V = struct('fname',fullfile('D:\Users\Clément\test_tMCimg\MCconfig',[file(1:end-4) '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);

V = spm_create_Yf(V);
V = spm_write_vol(V, Yf);


end