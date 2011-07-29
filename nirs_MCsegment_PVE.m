function out = nirs_MCsegment_PVE(job)

% boldmask = job.boldmask;
% XYZmm = job.XYZmm;
% mat = job.M;
boldmask = job.xSPM.Z;
XYZmm = job.xSPM.XYZmm;
% mat = job.xSPM.M;
xSPM = job.xSPM;

[dir,name,ext] = fileparts(job.T1seg);
V = spm_vol(job.T1seg);
Y = spm_read_vols(V);

% %%%%%% ici je fais le choix de reconstruire l image BOLD brute et de la
% %%%%%% resizer avec le module, Sinon on aurait pu utiliser V.mat de la
% %%%%%% T1seg et interpoler mais il fallait ecrire le code de l
% %%%%%% interpolation
% Ybold = zeros(xSPM.DIM');
% for i=1:size(xSPM.XYZ,2)
%     Ybold(xSPM.XYZ(1,i),xSPM.XYZ(2,i),xSPM.XYZ(3,i))=xSPM.Z(1,i);
% end
% 
% % Vbold = struct('fname',fullfile(job.cs_dir,'boldmask.nii'),...
% %     'dim',  xSPM.DIM',...
% %     'dt',   V.dt,...
% %     'pinfo',V.pinfo,...
% %     'mat',  xSPM.M);
% % Vbold = spm_create_vol(Vbold);
% % spm_write_vol(Vbold,Ybold);
% 
% Ybold_vxT1 = V.mat\(xSPM.M*Ybold);
% 
% V = spm_vol(outRS);
% Yb = spm_read_vols(V);
% threshold=0;
% Yf = Y + 6*(Yb>threshold);
% %%%%%% 
   
% % 
%%%%%% autre methode< On utilise V.mat et comme ca on est direct dans le
%%%%%% bon espace mais il faut interoler point a point (on peut pas utiliser resize)
% % count=0;
for i=1:size(XYZmm,2)
    XYZ_vxT1seg(:,i) = V.mat\[XYZmm(:,i);1];
end
Yb = zeros(V.dim);
%%%%%%%%%% interpolation
xi = XYZ_vxT1seg(1,:);
yi = XYZ_vxT1seg(2,:);
zi = XYZ_vxT1seg(3,:);

xf = 0:1:V.dim(1)-1;
yf = 0:1:V.dim(2)-1;
zf = 0:1:V.dim(3)-1;
% Make grids, to use with 3-dimensional interpolation
% [xi1,yi1,zi1] = meshgrid(xi,yi,zi);
[xf1,yf1,zf1] = meshgrid(xf,yf,zf);

rY = interp3(xi,yi,zi,xSPM.Z,xf1,yf1,zf1,'nearest');
Yf = Y + 6*(rY>0);
% disp([int2str(count) ' voxels (among ' int2str(size(xSPM.XYZ,2)) ')do not belong to the ROI...']);

Vpve = struct('fname',fullfile(dir,['PVE_' name ext]),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);
Vpve = spm_create_vol(Vpve);
spm_write_vol(Vpve,Yf);

out = fullfile(dir,['PVE_' name ext]);
end

