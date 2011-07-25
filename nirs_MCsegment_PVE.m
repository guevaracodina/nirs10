function out = nirs_MCsegment_PVE(job)

boldmask = job.boldmask;
XYZmm = job.XYZmm;
mat = job.M;

[dir,name,ext] = fileparts(job.T1seg);
V = spm_vol(job.T1seg);
Y = spm_read_vols(V);

% load(fullfile(dir,[name '-roi2raw.mat']));
count=0;
for i=1:size(XYZmm,2)
    XYZ(:,i) = round(V.mat\[XYZmm(:,i);1]);%- mat_roi2raw(:,4);
    try
        Y(XYZ(1,i),XYZ(2,i),XYZ(3,i)) = Y(XYZ(1,i),XYZ(2,i),XYZ(3,i))+6;
    catch
        count = count+1;
    end
end
disp([int2str(count) ' voxels (among ' int2str(size(XYZmm,2)) ')do not belong to the ROI...'])

Vpve = struct('fname',fullfile(dir,['PVE_' name ext]),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);
Vpve = spm_create_vol(Vpve);
spm_write_vol(Vpve,Y);

out = fullfile(dir,['PVE_' name ext]);
end

