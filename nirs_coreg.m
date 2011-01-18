function out = nirs_coreg(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: using previously created NIRS.mat containing information of
%optode and fiducial positions in one system of coordinates, this function
%will 1) normalize an anatomical image to an atlas (the T1.nii in SPM
%templates) that includes the fiducials (nasion, auricular left and right),
%2) apply the inverse normalization transform on the atlas fiducials,
%3) find a rotation that matches the subject fiducials to the atlas
%4) output an updated NIRS.mat that contains the transformations and
%the transformed coordinates

%Load NIRS.mat information
NIRS = job.NIRS;
%Store T1 file location
NIRS.anatT1 = job.anatT1{1,1};
NIRS.anatT1_template = job.anatT1_template{1,1};

%Various options that we don't make available to the user in the GUI
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {NIRS.anatT1};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {NIRS.anatT1};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {NIRS.anatT1_template};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50 78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';

spm_jobman('run_nogui',matlabbatch);

%Recreate name of _sn.mat file just created by spm_normalise, and load it
[pth,nam] = spm_fileparts(deblank(NIRS.anatT1));
sn_filename  = fullfile(pth,[nam,'_sn.mat']);
NIRS.wT1 = load(sn_filename);

%There are two physical objects: the MRI image, and the NIRS construct
%for the MRI image, it can be normalized to a standard atlas (Talairach Tournoux) or not
%for the MRI image, it can be in mm or in voxels
%for the NIRS construct, it comes in a DGT system
%In addition, we may want to project on the cortex or not
%Thus, upon alignment, quantities can be expressed
%   a) normalized or not
%   b) in mm or voxels
%   c) in MNI or DGT coordinates
%   d) surface positions can be on cortex or on skin / helmet / head
%For the MonteCarlo simulation, we want MNI in mm, not normalized, on skin

%Matrix Q maps unnormalized world coordinates to normalized world coordinates
%Affine maps unnormalized voxel coordinates to normalized voxel coordinates
%Note that transformations always act on the right, on column vectors of
%coordinates, returning column vectors
%Hence the need to transpose our matrices of coordinates in NIRS.mat
%Below: label with temp_ all the transposed coordinates to avoid confusion
Q = (NIRS.wT1.VG.mat/NIRS.wT1.Affine)/NIRS.wT1.VF.mat;

NIRS.FidPos_wMNI = [ job.nasion_wMNI; job.AL_wMNI; job.AR_wMNI ];
%call temp_... all the transposed positions
temp_FidPos_wMNI = [NIRS.FidPos_wMNI'; [1 1 1]];
temp_FidPos_MNI = Q\temp_FidPos_wMNI;
NIRS.FidPos_MNI = transpose(temp_FidPos_MNI(1:3,:));
y = NIRS.FidPos_MNI';

if iscell(NIRS.FidPos)
    x = cell2mat(NIRS.FidPos);
else
    x = NIRS.FidPos';
end

% best orientation is looked for
[s R t] = abs_orientation(x,y);     % y = s*R(x) + t
estY = zeros(size(y));
for iterData = 1:NIRS.n_Fid
    estY(:,iterData) = s*R*x(:,iterData) + t;
end;

err = y - estY;
errVal = sum(err(:).^2);
NIRS.errValofCoreg_mm2 = errVal;

if iscell(NIRS.SrcPos)
    temp_OptPos = [cell2mat(NIRS.SrcPos)' cell2mat(NIRS.DetPos)'];
else
    temp_OptPos = [NIRS.SrcPos' NIRS.DetPos'];
end
temp_OptPos_MNI = zeros(size(temp_OptPos));

optVoid = zeros(NIRS.n_Opt,1);
for OptIdx = 1:NIRS.n_Opt
    %check for Void sources (no data)
    if temp_OptPos(1,OptIdx) == 0 && temp_OptPos(2,OptIdx) == 0 && temp_OptPos(3,OptIdx) == 0
        optVoid(OptIdx,1) = 1;
    else
        temp_OptPos_MNI(:,OptIdx) = s*R*temp_OptPos(:,OptIdx) + t;
    end
end;

%MNI coordinates of optodes
NIRS.SrcPos_MNI = temp_OptPos_MNI(1:3,1:NIRS.n_Src)';
NIRS.DetPos_MNI = temp_OptPos_MNI(1:3,NIRS.n_Src+1:end)';

%unnormalized -> normalized, for optodes
temp_OptPos_wMNI = Q * [temp_OptPos_MNI;ones(1,NIRS.n_Opt)];          %% unit : mm
%projection on cortex - used to calculate normal directions of optodes to
%head surface, for Monte Carlo simulation
temp_OptPos_wMNI_cortex = projection_CS(temp_OptPos_wMNI);
temp_OptPos_MNI_cortex = zeros(4,NIRS.n_Opt);
for iterOpt = 1:NIRS.n_Opt
    if ~optVoid(iterOpt,1)
        %inversion: normalized -> unnormalized
        temp_OptPos_MNI_cortex(:,iterOpt) = Q\temp_OptPos_wMNI_cortex(:,iterOpt);
        %temp_OptPos_MNI_cortex(:,iterOpt) = temp_OptPos_MNI_cortex(1:3)';
        %inversion: MNI -> DGT
        %OptPos_cortex(:,iterOpt) = (1/s)*R'*(temp_OptPos_MNI_cortex(1:3) -t);
    end
end;
%temp_OptPos_wMNI_cortex = Q*temp_OptPos_MNI_cortex; %wT1_info.mat\temp_OptPos_MNI_cortex;
temp_OptPos_MNI_cortex = temp_OptPos_MNI_cortex(1:3,:);
%Normalize directions
temp_OptDir = temp_OptPos_MNI_cortex - temp_OptPos_MNI;
for iterOpt=1:NIRS.n_Opt
    temp_OptDir(:,iterOpt) = temp_OptDir(:,iterOpt);%/(sum(temp_OptDir(:,iterOpt).^2))^(1/2);
end
NIRS.SrcDir_MNI = temp_OptDir(:,1:NIRS.n_Src)';
NIRS.DetDir_MNI = temp_OptDir(:,NIRS.n_Src+1:NIRS.n_Opt)';
NIRS.SrcPos_MNI_cortex = temp_OptPos_MNI_cortex(:,1:NIRS.n_Src)';
NIRS.DetPos_MNI_cortex = temp_OptPos_MNI_cortex(:,NIRS.n_Src+1:NIRS.n_Opt)';
%To keep track of inexistent sources or detectors
NIRS.SrcVoid = optVoid(1:NIRS.n_Src);
NIRS.DetVoid = optVoid(NIRS.n_Src+1:end);

save(fullfile(NIRS.subj_path,'NIRS'),'NIRS');
out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');


function [s R t] = abs_orientation(x,y)
% y = sR(x) + t
% find s,  R,  t

N = size(x,2);

% centroid
cntrX = (sum(x')/N)';
cntrY = (sum(y')/N)';

% new coordinate
xNew = x - cntrX*ones(1,N);
yNew = y - cntrY*ones(1,N);

% make 'NN' matrix
NN = zeros(4,4);
for iterData = 1:N
    NN = NN + quaternion2(xNew(:,iterData))'*quaternion1(yNew(:,iterData));
end;

% find maximun eigenvector
[eigv eigd] = eig(NN);
qtn = eigv(:,4);

% estimate R
R = qtn2rtm(qtn);

% find scale factor
s = 0;
for iterData = 1:N
    ry_dot_R_rx = sum(yNew(:,iterData) .* (R * xNew(:,iterData)));
    s = s + ry_dot_R_rx;
end;
s = s / sum(xNew(:).^2);

% find translation factor
t = cntrY - s * R*cntrX;

function [R] = quaternion1(r)
if(length(r) == 3)
    R = [0 -r(1) -r(2) -r(3);...
        r(1) 0 -r(3) r(2);...
        r(2) r(3) 0 -r(1);...
        r(3) -r(2) r(1) 0;];
elseif(length(r) == 4)
    R = [r(1) -r(2) -r(3) -r(4);...
        r(2) r(1) -r(4) r(3);...
        r(3) r(4) r(1) -r(2);...
        r(4) -r(3) r(2) r(1)];
else
    disp('Error @ quaternion1');
end;

function [R] = quaternion2(r)
if(length(r) == 3)
    R = [0 -r(1) -r(2) -r(3);...
        r(1) 0 r(3) -r(2);...
        r(2) -r(3) 0 r(1);...
        r(3) r(2) -r(1) 0];
elseif(length(r) == 4)
    R = [r(1) -r(2) -r(3) -r(4);...
        r(2) r(1) r(4) -r(3);...
        r(3) -r(4) r(1) r(2);...
        r(4) r(3) -r(2) r(1)];
else
    disp('Error @ quaternion2');
end;

function [R] = qtn2rtm(q)
% quaternion to rotation matrix
R = zeros(3,3);
s = q(1);
x = q(2);
y = q(3);
z = q(4);
R(1,1) = s^2 + x^2 - y^2 - z^2;
R(1,2) = 2 * (x*y - s*z);
R(1,3) = 2 * (x*z + s*y);
R(2,1) = 2 * (y*x + s*z);
R(2,2) = s^2 - x^2 + y^2 - z^2;
R(2,3) = 2 * (y*z - s*x);
R(3,1) = 2 * (z*x - s*y);
R(3,2) = 2 * (z*y + s*x);
R(3,3) = s^2 - x^2 - y^2 + z^2;