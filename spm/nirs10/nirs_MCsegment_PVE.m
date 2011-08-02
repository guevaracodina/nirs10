function out = nirs_MCsegment_PVE(job)
% Cl2;ent Bonnéry 01/08/2011

xSPM = job.xSPM;

[dir,name,ext] = fileparts(job.T1seg);
V = spm_vol(job.T1seg);
Y = spm_read_vols(V);

Ybold = zeros(xSPM.DIM');
for i=1:size(xSPM.XYZ,2)
    Ybold(xSPM.XYZ(1,i),xSPM.XYZ(2,i),xSPM.XYZ(3,i))=xSPM.Z(1,i);
end

% following lines equivalent to rY_vxT1 = V.mat\(xSPM.M*rY);
affine = V.mat\xSPM.M;
T = maketform('affine',affine');
R = makeresampler('nearest','bound');
rY_vxT1 = tformarray(Ybold,T',R,[1 2 3],[1 2 3],V.dim,[],[]);

% segments the BOLD image
Yf = Y + 6*(rY_vxT1>0);

pve_n = fullfile(job.cs_dir,['PVE' name ext]);
Vpve = struct('fname',pve_n,...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);
Vpve = spm_create_vol(Vpve);
spm_write_vol(Vpve,Yf);

out = pve_n;
end

% % % 
% % Vbold = struct('fname',fullfile(job.cs_dir,'boldmask.nii'),...
% %     'dim',  xSPM.DIM',...
% %     'dt',   V.dt,...
% %     'pinfo',V.pinfo,...
% %     'mat',  xSPM.M);
% % Vbold = spm_create_vol(Vbold);
% % spm_write_vol(Vbold,Ybold);
% % 
% xi = 0:1:xSPM.DIM(1)-1;
% yi = 0:1:xSPM.DIM(2)-1;
% zi = 0:1:xSPM.DIM(3)-1;
% % 
% % % 
% % % V = spm_vol(outRS);
% % % Yb = spm_read_vols(V);
% % % threshold=0;
% % % Yf = Y + 6*(Yb>threshold);
% % % %%%%%% 
% %    
% % % % 
% % %%%%%% autre methode< On utilise V.mat et comme ca on est direct dans le
% % %%%%%% bon espace mais il faut interoler point a point (on peut pas utiliser resize)
% % % % count=0;
% % for i=1:size(XYZmm,2)
% %     XYZ_vxT1seg(:,i) = V.mat\[XYZmm(:,i);1];
% %     Ybold_vxT1seg(XYZ_vxT1seg(1,i),XYZ_vxT1seg(2,i),XYZ_vxT1seg(3,i))=xSPM.Z(1,i);
% % end
% % Yb = zeros(V.dim);
% % %%%%%%%%%% interpolation
% % xi = XYZ_vxT1seg(1,:);
% % yi = XYZ_vxT1seg(2,:);
% % zi = XYZ_vxT1seg(3,:);
% 
% xf = 0:1:V.dim(1)-1;
% yf = 0:1:V.dim(2)-1;
% zf = 0:1:V.dim(3)-1;
% % Make grids, to use with 3-dimensional interpolation
% [xi1,yi1,zi1] = meshgrid(xi,yi,zi);
% [xf1,yf1,zf1] = meshgrid(xf,yf,zf);
% 
% Ybold_p = permute(Ybold,[2,1,3]);
% rY_p = interp3(xi1,yi1,zi1,Ybold_p,xf1,yf1,zf1,'nearest');
% rY = permute(rY_p,[2,1,3]);