function vol = nirs_read_2pt(file,nx,ny,nz,nt)
% Read and sum 2pt files
% Any other initializations would go here
fid = fopen(file, 'rb');
vol=zeros(nx*ny*nz,1);
for index=1:nt
    vol=vol+fread(fid,nx*ny*nz,'double');
end
fclose(fid);

vol=reshape(vol,[ny nx nz]);

% Some values are negative... should know why
vol=(vol>0).*vol;

%%% mise sous forme de nifti :

V.dim = [ny nx nz];

V.mat = eye(4);
V.mat(1:3,4) = [round(ny/2);round(nx/2);round(nz/2)];
V.dt = [4,0];
V.pinfo = [1;0;352];


V = struct('fname',fullfile('D:\Users\Clément\test_tMCimg\MCconfig',[file(1:end-4) '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);

V = spm_create_vol(V);
V = spm_write_vol(V, vol);