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
viewer_ON = job.View6Projections;
Save6Projections = job.Save6Projections;
ForceReprocess = job.ForceReprocess;
template_mode = job.template_mode;
% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{iSubj,1} = newNIRSlocation;
        % SPATIAL NORMALIZATION OF ANATOMICAL IMAGE %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'coregOK') || job.force_redo)
            % Allow user-specified image of subject to overwrite previous
            % anatomical image in NIRS.mat; unlikely to ever happen
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
                    [DirSPM,dummy,dummy2] = fileparts(which('spm'));
                    NIRS.Dt.ana.tT1 = fullfile(DirSPM,'templates','T1.nii');
                end
            end
            if ~template_mode || (template_mode && iSubj == 1)
                [dirT1, fil, ext] = fileparts(NIRS.Dt.ana.T1);
                fwT1 = [dirT1 filesep 'w' fil ext(1:4)];
                if spm_existfile(fwT1) && ~ForceReprocess
                    disp(['Spatial normalisation already run in file ' fwT1 ' - skipping']);
                else
                    tmp_file = fullfile(dirT1,['m' fil ext]);
                    if exist(tmp_file,'file')
                        src_file = tmp_file;
                    else
                        src_file = NIRS.Dt.ana.T1;
                    end
                    %Various options that we don't make available to the user in the GUI
                    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {src_file};
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
                    spm_jobman('run',matlabbatch);
                    
                end
                
                %Recreate name of _sn.mat file just created by spm_normalise, and load it
                [pth,nam] = spm_fileparts(deblank(NIRS.Dt.ana.T1));
                sn_filename  = fullfile(pth,[nam '_sn.mat']);
                if ~exist(sn_filename,'file')
                    sn_filename  = fullfile(pth,['m' nam '_sn.mat']);
                end
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
                
                % FIT ROM TO RMM POSITIONS, USING FIDUCIALS %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (isfield(NIRS.Dt.fir,'stax') && strcmp(NIRS.Dt.fir.stax.n,'Brainsight(c)')) ||...
                        ~isfield(NIRS.Dt.fir,'stax')
                    
                    % Positions of fiducial points
                    if job.fid_in_subject_MNI
                        %quick fix to use previously obtained coordinates in subject MNI
                        %coordinates
                        NIRS.Cf.H.F.r.m.mm.p  = [job.nasion_wMNI' job.AL_wMNI' job.AR_wMNI'];
                    else
                        NIRS.Cf.H.F.w.m.mm.p  = [job.nasion_wMNI' job.AL_wMNI' job.AR_wMNI'];
                        Fp_wmm = NIRS.Cf.H.F.w.m.mm.p;
                        %call temp_... all the transposed positions
                        temp_Fp_wmm = [Fp_wmm; [1 1 1]];
                        temp_Fp_rmm = Q\temp_Fp_wmm;
                        NIRS.Cf.H.F.r.m.mm.p = temp_Fp_rmm(1:3,:);
                    end
                    y = NIRS.Cf.H.F.r.m.mm.p;
                    x = NIRS.Cf.H.F.r.o.mm.p;
                    
                    % Compute optimal rigid transformation to match fiducial
                    % positions to those of the normalized atlas
                    [s R t] = abs_orientation(x,y); % y = s*R(x) + t
                    estY = zeros(size(y));
                    for Fi = 1:3
                        estY(:,Fi) = s*R*x(:,Fi) + t;
                    end;
                    % Store coregistration error
                    err = y - estY;
                    errVal = sum(err(:).^2);
                    errValMax = (max(sum(err.^2,1))^0.5)/10;
                    err %show error on each fiducial
                    disp(['Error Value for subject ' int2str(iSubj) ': ' num2str(errVal)]);
                    disp(['Worst coregistration error: ' num2str(errValMax) ' cm']);
                    NIRS.Dt.pro.errValofCoreg_mm2 = errVal;
                    NIRS.Dt.pro.errValofCoreg_mm2_all = err;
                    NIRS.Dt.pro.errValofCoreg_cm_worst = errValMax;
                    % Apply the same transformation to all points in order to
                    % achieve coregistration
                    try Sp_rom = NIRS.Cf.H.S.r.o.mm.p; end
                    try Dp_rom = NIRS.Cf.H.D.r.o.mm.p; end
                    try Qp_rom = NIRS.Cf.H.Q.r.o.mm.p; end
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
                        % check for Void sources (no data)
                        if Pp_rom(1,Pi) == 0 && Pp_rom(2,Pi) == 0 && Pp_rom(3,Pi) == 0
                            Pvoid(Pi,1) = 1;
                        else
                            Pp_rmm(:,Pi) = s*R*Pp_rom(:,Pi) + t;
                        end
                    end;
                    
                    % Save MNI coordinates of optodes
                    NIRS.Cf.H.P.r.m.mm.p = Pp_rmm;
                    
                else % For vitamin markers based coregistration, the positions were
                    % obtained in detectVitamins module; those are already defined
                    % in the MNI space (i voxels relative to anatomical image)
                    Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
                    % Initialize other variables used later in the code
                    NP = size(Pp_rmm,2);
                    Pvoid = zeros(1,NP);
                    
                end
                
                % unnormalized -> normalized, for optodes
                Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm
                
                % PROJECT OPTODE POSITIONS ON CORTEX SURFACE %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Projection on cortex - used to calculate normal directions of
                % optodes to head surface, for Monte Carlo simulation
                Pp_c1_wmm = projection_CS(Pp_wmm);
                Pp_c1_rmm = zeros(4,NP);
                for Pi = 1:NP
                    if ~Pvoid(1,Pi)
                        %inversion: normalized -> unnormalized
                        Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi);
                    end
                end;
                Pp_c1_rmm = Pp_c1_rmm(1:3,:);
                
                % Save Normalized coordinates of optodes
                NIRS.Cf.H.P.w.m.mm.p = Pp_wmm;
                NIRS.Cf.H.P.w.m.mm.c1.p = Pp_c1_wmm;
                % Save MNI coordinates of optodes on cortex (c1)
                NIRS.Cf.H.P.r.m.mm.c1.p = Pp_c1_rmm;
                NIRS.Cf.H.P.void = Pvoid;
                
                % PROJECT OPTODE POSITIONS ON SKIN SURFACE %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Calculation of fitted positions on the skin surface
                % Fitted positions (positions are fitted with respect to the scalp)
                
                %%%%%%%%%%%% CB: A RETRAVAILLER : SI ON VEUT VRAIMENT QUE CA SOIT
                %%%%%%%%%%%% UTILE, IL FAUT POUVOIR RECUPERER CE QUI EST AFFICHE
                %%%%%%%%%%%% SUR LE PDF DE LA NORMALISATION : en fait, lorsqu'on
                %%%%%%%%%%%% nor;qlise on obtient une ;atrice de passage entre
                %%%%%%%%%%%% l'image anato;ique du sujet nor;alisee et la te;plate.
                %%%%%%%%%%%% Cette matrice peut etre utilisee ici pour faire le
                %%%%%%%%%%%% lien entre les deux images. Mais si on a une autre
                %%%%%%%%%%%% image que la template, ca devient impossible...
                if isempty(job.segT1_4fit{1,1})
                    try
                        fsegT1_4fit = NIRS.Dt.ana.T1seg;
                    catch
                        fsegT1_4fit = NIRS.Dt.ana.T1;
                        %                 disp('Could not find a segmented image to fit positions on scalp.');
                    end
                else
                    % Store segmented image file location
                    NIRS.Dt.ana.T1seg = job.segT1_4fit{1,1};
                    fsegT1_4fit = NIRS.Dt.ana.T1seg;
                end
                
                jobe.Pp_rmm = Pp_rmm;
                jobe.Pp_c1_rmm = Pp_c1_rmm;
                jobe.NP = NP;
                jobe.image_in = {fsegT1_4fit};
                jobe.lby = 'coreg';
                out2 = nirs_fit_probe(jobe); %careful, this is not the same out as the
                %out of the function with NIRS.mat!
                Pfp_rmm = out2{1};
                % from MNI real space (mm) to MNI voxel space
                V_4fit = spm_vol(fsegT1_4fit);
                Pfp_rmv = [];
                for i=1:NP
                    Pfp_rmv(:,i) = V_4fit.mat\[Pfp_rmm(:,i);1];
                end
                Pfp_rmv = Pfp_rmv(1:3,:);
                % Save MNI and voxel coordinates of optodes, fitted on skin surface
                NIRS.Cf.H.P.r.m.mm.fp = Pfp_rmm;
                NIRS.Cf.H.P.r.m.vx.fp = Pfp_rmv;
                
                % Save all changes to NIRS matrix
                save(job.NIRSmat{iSubj},'NIRS');
                NIRS =[];
                
                % GENERATE TOPO DATA %
                %%%%%%%%%%%%%%%%%%%%%%
                
                load(job.NIRSmat{iSubj});
                if job.GenDataTopo
                    [pth2,dummy,dummy2] = fileparts(job.NIRSmat{iSubj,1});
                    [dir1,fil1,ext1] = fileparts(NIRS.Dt.ana.T1);
                    fwT1 = fullfile(dir1,['w' fil1 ext1]);
                    NIRS.Dt.ana.fwT1 = fwT1;
                    try
                        %V0 = spm_vol(NIRS.Dt.ana.T1);
                        V = spm_vol(fwT1);
                        wT1_info.mat = V.mat;%[-1.5 0 0 79; 0 1.5 0 -113; 0 0 1.5 -51; 0 0 0 1]; %V.mat;
                        wT1_info.dim = V.dim;% [105,126,91]; %V.dim;
                        %let's use positions of optodes on cortex
                        %loop over channel ids
                        Nch0 = size(NIRS.Cf.H.C.id,2)/2;
                        ch_MNIw_vx = zeros(4,Nch0); %%%% la notation serait Cp_rmv
                        ch_MNI_vx = zeros(4,Nch0);
                        %ch_MNI_vx_skin = zeros(4,Nch0);
                        %number of sources
                        Ns = NIRS.Cf.H.S.N;
                        for i=1:Nch0
                            %indices of source and detector
                            Si = NIRS.Cf.H.C.id(2,i);
                            Di = NIRS.Cf.H.C.id(3,i)+Ns;
                            posmm = [(Pp_c1_rmm(:,Si)+Pp_c1_rmm(:,Di))/2;1];
                            pos = NIRS.Dt.ana.wT1.VF.mat\posmm;
                            ch_MNI_vx(:,i) = pos;
                            %posmm_skin = [(Pp_rmm(:,Si)+Pp_rmm(:,Di))/2;1];
                            %pos_skin = NIRS.Dt.ana.wT1.VF.mat\posmm_skin;
                            %ch_MNI_vx_skin(:,i) = pos_skin;
                            posw = V.mat\(Q*posmm);
                            ch_MNIw_vx(:,i) = posw; %%%% la notation serait Cp_rmv
                        end
                        [rend, rendered_MNI] = render_MNI_coordinates(ch_MNIw_vx,ch_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,job.render_choice);
                        %[rend, rendered_MNI] = render_MNI_coordinates(ch_MNIw_vx,ch_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,job.render_choice);
                        
                        %%%% la notation serait Cp_rmv
                        for kk = 1:6
                            rendered_MNI{kk}.ren = rend{kk}.ren;
                        end
                        rend_file = fullfile(pth2,'TopoData.mat');
                        save(rend_file, 'rendered_MNI');
                        NIRS.Dt.ana.rend = rend_file;
                        
                        %Viewer from NIRS_SPM
                        %viewer_ON = 0;
                        if viewer_ON || Save6Projections
                            %rendered_MNI = varargin{1};
                            Nch = size(rendered_MNI{1}.rchn,1);
                            load Split
                            for kk=1:6
                                fh0(kk) = figure;
                                %kk = 4; % dorsal view
                                %brain = rend{kk}.ren;
                                brain = rendered_MNI{kk}.ren;
                                brain = brain.* 0.5;
                                sbar = linspace(0, 1, 128);
                                sbrain = ((-sbar(1) + sbar(64))/(0.5)).* brain + sbar(1);
                                sbrain(1,1) = 1;
                                rchn = rendered_MNI{kk}.rchn;
                                cchn = rendered_MNI{kk}.cchn;
                                for jj = 1:Nch
                                    if rchn(jj) ~= -1 && cchn(jj) ~= -1 %% updated 2009-02-25
                                        if rchn(jj) < 6 || cchn(jj) < 6
                                            sbrain(rchn(jj), cchn(jj)) = 0.9; % 0.67
                                        else
                                            sbrain(rchn(jj)-5:rchn(jj)+5, cchn(jj)-5:cchn(jj)+5) = 0.9;
                                        end
                                    end
                                end
                                imagesc(sbrain);
                                colormap(split);
                                axis image;
                                axis off;
                                
                                for jj = 1:Nch
                                    if rchn(jj) ~= -1 && cchn(jj) ~= -1 %% updated 2009-02-25
                                        text(cchn(jj)-5, rchn(jj), num2str(jj), 'color', 'r');
                                    end
                                end
                                try
                                    if Save6Projections
                                        [side_hemi spec_hemi] = nirs_get_brain_view(kk);
                                        filen1 = fullfile(pth2,[spec_hemi '.fig']);
                                        saveas(fh0(kk),filen1,'fig');
                                        filen2 = fullfile(pth2,[spec_hemi '.tif']);
                                        print(fh0(kk), '-dtiffn', filen2);
                                        close(fh0(kk));
                                    end
                                end
                            end
                        end
                    catch
                        disp('Could not create TopoData.mat file');
                    end
                end
                NIRS.flags.coregOK = 1;
                if template_mode
                    NIRScoreg = NIRS;
                end
            else
                %copy coregistration results for each additional subject
                if ~exist('NIRScoreg','var')
                    tmpNIRS = NIRS;
                    load(job.NIRSmat{1,1})
                    NIRScoreg = NIRS;
                    NIRS = tmpNIRS;
                end
                NIRS.Dt.ana.rend = NIRScoreg.Dt.ana.rend; %enough to have just TopoData
                NIRS.flags.coregOK = 1;
            end
        end
        save(job.NIRSmat{iSubj,1},'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Coregistration failed for the ' int2str(iSubj) 'th subject. for ' job.NIRSmat{iSubj,1}]);
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