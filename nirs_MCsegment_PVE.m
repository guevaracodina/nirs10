function out = nirs_MCsegment_PVE(job)

boldmask = job.boldmask;
XYZ = job.XYZ;

[dir,name,ext] = fileparts(job.T1seg);
V = spm_vol(job.T1seg);
Y = spm_read_vols(V);

for i=1:size(boldmask,2)
    Y(XYZ(1,i),XYZ(2,i),XYZ(3,i)) = Y(XYZ(1,i),XYZ(2,i),XYZ(3,i))+6;
end

Vpve = struct('fname',fullfile(dir,['PVE_' name ext]),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);
Vpve = spm_create_vol(Vpve);
spm_write_vol(Vpve,Y);

out = fullfile(dir,['PVE_' name ext]);
end

