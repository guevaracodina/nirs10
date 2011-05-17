function out = nirs_run_buildroi2(job)
% Builds the smallest ROI containing all the selected source-detector pairs
% FORMAT NIRSmat = nirs_run_buildroi2(NIRS, image_in, keepAllChannels, outputprefix)
% NIRS            - NIRS matrix
% image_in        - segmented image from which the ROI will be extracted
% keepAllChannels - Channels to keep in the ROI
% outputprefix    - Prefix of the ROI image name
%_______________________________________________________________________
%
% The image will be croped with respect to the selected source-detector 
% pairs (channels) to are to be kept for the MonteCarlo simulation. You
% only have to enter the channels numbers for the first wavelength.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2011-03

load(job.NIRSmat{:});

V = spm_vol(job.image_in{:});
Y = spm_read_vols(V);
[dir,name] = fileparts(V.fname);
NS = NIRS.Cf.H.S.N;
Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
Cid = NIRS.Cf.H.C.id;
NC = NIRS.Cf.H.C.N;
wl = NIRS.Cf.dev.wl;
nc = NC/length(wl);
try
    Ckpt = job.keepAllChannels.keepChannels;
    for iwl=2:length(wl)
        Ckpt = [Ckpt Ckpt+(iwl-1)*nc];
    end
catch
    Ckpt = Cid(1,:);
end
Skpt = unique(Cid(2,Ckpt));
Dkpt = unique(Cid(3,Ckpt));
Pkpt=[Skpt Dkpt+NS];

%%%  vérifier que ça a de l'intérêt de garder tout ça
for i=1:size(Pkpt,2)
    Pfp_roi_rmm(:,i) = Pfp_rmm(:,Pkpt(i));
    Pp_roi_rmm(:,i) = Pp_rmm(:,Pkpt(i));
    Pp_roi_c1_rmm(:,i) = Pp_c1_rmm(:,Pkpt(i));
end

Pfp_roi_rmvtemp = zeros(4,size(Pfp_roi_rmm,2));
for i=1:size(Pfp_roi_rmm,2)
    Pfp_roi_rmvtemp(:,i) = V.mat\[Pfp_roi_rmm(:,i);1];
end
Pfp_roi_rmv = Pfp_roi_rmvtemp(1:3,:);

bbv(1,1) = min(Pfp_roi_rmv(1,:));
bbv(1,2) = max(Pfp_roi_rmv(1,:));

bbv(2,1) = min(Pfp_roi_rmv(2,:));
bbv(2,2) = max(Pfp_roi_rmv(2,:));

bbv(3,1) = min(Pfp_roi_rmv(3,:));
bbv(3,2) = max(Pfp_roi_rmv(3,:));

bbv = round(bbv);

% the size of the plotted image can be bigger than the size read in
% the header, in such case the value kept is the one of header
marge = 20;
for i =1:3
    bbv(i,1) = max(1,bbv(i,1)-marge);
    bbv(i,2) = min(V.dim(i),bbv(i,2)+marge);
end

% With resizing
%  mat_roi2raw is the transformation matrix from the voxel space
%  of the ROI to the voxel space of the raw image
% Translation : de l'image brute a la roi
mat_roi2raw = eye(4) + [zeros(4,3) [bbv(:,1);0]];

V_roi.dim = [bbv(1,2)-bbv(1,1)+1 bbv(2,2)-bbv(2,1)+1 bbv(3,2)-bbv(3,1)+1];
% V_roi.mat is the mat from voxel ROI to the raw image space in mm
V_roi.mat = V.mat*mat_roi2raw;

V_roi = struct('fname',fullfile(dir,[job.output_prefix,name,'.nii']),...
    'dim',  V_roi.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V_roi.mat);

Y_roi = Y(bbv(1,1):bbv(1,2),bbv(2,1):bbv(2,2),bbv(3,1):bbv(3,2));

V_roi = spm_create_vol(V_roi);
V_roi = spm_write_vol(V_roi, Y_roi);

% Sources and detectors positions must be updated also
for i=1:size(Pkpt,2)
    PfpR_roi_rmv(:,i) = [Pfp_roi_rmv(:,i);1] - mat_roi2raw(:,4);
end

%then saved
if ~isfield(NIRS,'Cs'), NIRS.Cs={}; end
if isfield(NIRS.Cs,'temp'), clear NIRS.Cs.temp; end
NIRS.Cs.temp.Pkpt = Pkpt;
NIRS.Cs.temp.NSkpt = size(Skpt,2);
NIRS.Cs.temp.NDkpt = size(Dkpt,2);
NIRS.Cs.temp.Pfp_roi_rmv = PfpR_roi_rmv(1:3,:);
NIRS.Cs.temp.Pfp_roi_rmm = Pfp_roi_rmm;
NIRS.Cs.temp.Pp_roi_rmm = Pp_roi_rmm;
NIRS.Cs.temp.Pp_roi_c1_rmm = Pp_roi_c1_rmm;
NIRS.Cs.temp.segR = fullfile(dir,[job.output_prefix,name,'.nii']);
save(job.NIRSmat{:},'NIRS');

out.NIRSmat = job.NIRSmat;
end