function out = nirs_resize(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Clément Bonnéry 11/2010

P = job.image_in{:};
if iscell(P), P = P{:};end
PO_dim = job.out_dim;
PO_dt = job.out_dt;
PO_autonaming = job.out_autonaming;
try PO_prefix = job.out_prefix; catch, end
PO_dir = job.out_dir;

V = spm_vol(P);
dimi = V.dim;

if PO_dim == [1 1 1];%isotropic voxels
    vxsize = job.out_vxsize;
    
    inv_mat = spm_imatrix(V.mat);
    scalings = diag(inv_mat(7:9)); %see spm_matrix for definition
    
    dimf = ceil((dimi-1)/vxsize * abs(scalings));
    for i =1:3, RZS(:,i) = vxsize*(V.mat(:,i)./norm(V.mat(:,i))); end
    rV.mat = [RZS V.mat(:,4)];
else
    dimf = PO_dim;
    rV.mat =  V.mat*diag([V.dim./dimf 1]);
end

% Create linear spaces for initial and final volume
xi = linspace(0,dimf(1)-1,dimi(1));
yi = linspace(0,dimf(2)-1,dimi(2));
zi = linspace(0,dimf(3)-1,dimi(3));
xf = 0:1:dimf(1)-1;
yf = 0:1:dimf(2)-1;
zf = 0:1:dimf(3)-1;

% Make grids, to use with 3-dimensional interpolation
[xi1,yi1,zi1] = meshgrid(xi,yi,zi);
[xf1,yf1,zf1] = meshgrid(xf,yf,zf);

Y = spm_read_vols(V);
Y = permute(Y,[2,1,3]);

rY = interp3(xi1,yi1,zi1,Y,xf1,yf1,zf1,'nearest');

% preparation et ecriture
if PO_autonaming==0
    PO_prefix = [int2str(dimf(1,1)) 'x' int2str(dimf(1,2)) 'x' int2str(dimf(1,3)) '_'];
end
[dummy,nam,ext] = fileparts(V.fname);
rV.fname = fullfile(PO_dir,[PO_prefix nam ext]);

switch PO_dt
    case 'same'
        rV.dt = V.dt;
    case 'uint16'
        rV.dt   = [16 0]; %float32 not V.dt which is int16
end
rV.dim  = dimf;
rV.pinfo = V.pinfo;

rY=permute(rY,[2,1,3]);
spm_write_vol(rV,rY);

out = rV.fname;
end