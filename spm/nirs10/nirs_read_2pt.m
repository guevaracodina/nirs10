function vol = nirs_read_2pt(file,nx,ny,nz,nt)
% Read and sum 2pt files
% Any other initializations would go here
fid = fopen(file, 'rb');
[pth,nme,~] = fileparts(file);
vol=zeros(nx*ny*nz,1);
for index=1:nt
    vol=vol+fread(fid,nx*ny*nz,'double');
end
fclose(fid);

vol=reshape(vol,[nx ny nz]);%vol=reshape(vol,[ny nx nz]);
% Some values are negative... should know why
eps = 0.00001;
vol=log(abs((vol>0).*vol)+eps);

Vbase = spm_vol('D:\Users\Clément\Projet_ReML\donnees\test_roi\MCconfigHopingpongPINGPOUNG\99x95x80_roi_00044_segmented_s201007051500-0002-00001-000160-01.nii');

% mise sous forme de nifti :
V.dim = [nx ny nz];%V.dim = [ny nx nz];
V.mat = Vbase.mat;
V.dt = [4,0];
V.pinfo = [1;0;352];

V = struct('fname',fullfile(pth,[nme '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);

V = spm_create_vol(V);
V = spm_write_vol(V, vol);