function out = nirs_run_coreg(job)
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

for Idx=1:size(job.NIRSmat,1)

    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        %Allow user-specified image of subject to overwrite previous
        %anatomical image in NIRS.mat; unlikely to ever happen
        if isempty(job.anatT1{1,1})        
            try
                NIRS.Dt.ana.T1;
            catch
                disp('Could not find an anatomical image');
            end
        else
            %Store T1 file location
            NIRS.Dt.ana.T1 = job.anatT1{1,1};
        end
        try 
            tmpf = job.anatT1_template{1,1};
            if spm_existfile(tmpf)
                NIRS.Dt.ana.tT1 = tmpf;
            else
                [DirSPM,~,~] = fileparts(which('spm')); 
                NIRS.Dt.ana.tT1 = fullfile(DirSPM,'templates','T1.nii');  
            end
        end
        %Various options that we don't make available to the user in the GUI
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {NIRS.Dt.ana.T1};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {NIRS.Dt.ana.T1};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {NIRS.Dt.ana.tT1};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50
                                                                     78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';

        spm_jobman('run_nogui',matlabbatch);

        %Recreate name of _sn.mat file just created by spm_normalise, and load it
        [pth,nam] = spm_fileparts(deblank(NIRS.Dt.ana.T1));
        sn_filename  = fullfile(pth,[nam,'_sn.mat']);
        NIRS.Dt.ana.wT1 = load(sn_filename);

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
        %Affine maps unnormalized voxelcoordinates to normalized voxel coordinates
        %Note that transformations always act on the right, on column vectors of
        %coordinates, returning column vectors
        %Hence the need to transpose our matrices of coordinates in NIRS.mat
        %Below: label with temp_ all the transposed coordinates to avoid confusion
        Q = (NIRS.Dt.ana.wT1.VG.mat/NIRS.Dt.ana.wT1.Affine)/NIRS.Dt.ana.wT1.VF.mat;

        NIRS.Cf.H.F.w.m.mm.p  = [ job.nasion_wMNI;
            job.AL_wMNI;
            job.AR_wMNI ];
        Fp_wmm = NIRS.Cf.H.F.w.m.mm.p;
        %call temp_... all the transposed positions
        temp_Fp_wmm = [Fp_wmm'; [1 1 1]]; 
        temp_Fp_rmm = Q\temp_Fp_wmm;
        NIRS.Cf.H.F.r.m.mm.p = transpose(temp_Fp_rmm(1:3,:));

        y = NIRS.Cf.H.F.r.m.mm.p'; %bizarre - why is y transposed but not x???
        x = NIRS.Cf.H.F.r.o.mm.p'; %PP I transposed x too

        [s R t] = abs_orientation(x,y); % y = s*R(x) + t
        estY = zeros(size(y));
        for Fi = 1:3
            estY(:,Fi) = s*R*x(:,Fi) + t;
        end;
        err = y - estY;
        errVal = sum(err(:).^2);
        %double2str does not exist in my version of Matlab
        disp(['Error Value for subject ' int2str(Idx) ': ' num2str(errVal)]);
        NIRS.Dt.pro.errValofCoreg_mm2 = errVal;

        try Sp_rom = NIRS.Cf.H.S.r.o.mm.p'; end %PP transposed
        try Dp_rom = NIRS.Cf.H.D.r.o.mm.p'; end %PP transposed
        try Qp_rom = NIRS.Cf.H.Q.r.o.mm.p'; end%PP transposed
        try 
            Pp_rom = [Sp_rom Dp_rom Qp_rom];
        catch 
            Pp_rom = [Sp_rom Dp_rom];
        end
        NP = size(Pp_rom,2);
        NIRS.Cf.H.P.N = NP;

        Pp_rmm = zeros(size(Pp_rom));
        Pvoid = zeros(1,size(Pp_rom,2));

        for Pi = 1:NP
            %check for Void sources (no data)
            if Pp_rom(1,Pi) == 0 && Pp_rom(2,Pi) == 0 && Pp_rom(3,Pi) == 0
                Pvoid(Pi,1) = 1;
            else
                Pp_rmm(:,Pi) = s*R*Pp_rom(:,Pi) + t;
            end
        end;

        %Save MNI coordinates of optodes
        NIRS.Cf.H.P.r.m.mm.p = Pp_rmm'; %PP transposed

        %unnormalized -> normalized, for optodes
        Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm
        %projection on cortex - used to calculate normal directions of optodes to
        %head surface, for Monte Carlo simulation
        Pp_c1_wmm = projection_CS(Pp_wmm);
        Pp_c1_rmm = zeros(4,NP);
        for Pi = 1:NP   
            if ~Pvoid(1,Pi)
                %inversion: normalized -> unnormalized
                Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi);
            end
        end; 
        Pp_c1_rmm = Pp_c1_rmm(1:3,:);

        %Save MNI coordinates of optodes on cortex (c1)
        NIRS.Cf.H.P.r.m.mm.c1.p = Pp_c1_rmm'; %PP transposed
        NIRS.Cf.H.P.void = Pvoid'; %PP transposed
        save(job.NIRSmat{Idx},'NIRS');
    catch
        disp(['Coregistration failed for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;


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
for(iterData = 1:N)
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