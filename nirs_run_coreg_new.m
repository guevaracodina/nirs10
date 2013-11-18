function out = nirs_run_coreg_new(job)
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
%And more...
ForceReprocessNormalization = job.ForceReprocess;
if isfield(job.render_choice,'render_template')
    render_template = 1;
else
    render_template = 0; %render to subject
end
if isfield(job,'cortex_projection_method')
    CPM = job.cortex_projection_method;
    if isfield(CPM,'project_liom')
        cpm = 1;
    else
        if isfield(CPM,'project_liom_Ke')
            cpm = 2;
        else
            if isfield(CPM,'project_Korean')
                cpm = 3;
            end
        end
    end
else
    cpm = 3;
end
if isfield(job,'OutputSkinFigs')
    OutputSkinFigs = job.OutputSkinFigs;
else
    OutputSkinFigs = 1;
end
ERad = job.radius_channel;
ForceReprocessSegmentation = 0;
coreg_projected_channels_on_cortex_rather_than_midpoint = 1;
use_fSeg = 1; %option to use Clément's segmented image to extract c1; otherwise default to c1 image
% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
        if ~isfield(NIRS,'flags'), NIRS.flags = []; end
        job.NIRSmat{iSubj,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'coregOK') || job.force_redo)
            [dir_coreg,dummy] = fileparts(newNIRSlocation);
            %Normalization of anatomical image
            anatT1 = NIRS.Dt.ana.T1;
            [dirT1, fil, ext] = fileparts(anatT1);
            fwT1 = fullfile(dirT1,['w' fil ext(1:4)]);
            fc1 =  fullfile(dirT1,['c1' fil ext(1:4)]);
            %check if segmentation is required
            if ~(spm_existfile(fc1) && ~ForceReprocessSegmentation)
                nirs_batch_segment(newNIRSlocation);
                NIRS.Dt.ana.T1seg = fullfile(dirT1,['00021_segmented_' fil ext]); %This name needs to be
                %generalized or simplified -- as such, it will lead to a
                %bug when people select options in Clément's segmentation
            end
            try
                fSeg = NIRS.Dt.ana.T1seg;
            catch
                NIRS.Dt.ana.T1seg = fullfile(dirT1,['00021_segmented_' fil ext]);
            end
            
            if ~(spm_existfile(fwT1) && ~ForceReprocessNormalization)
                nirs_batch_normalize(anatT1,dir_coreg);
            end
            %Recreate name of _sn.mat file just created by spm_normalise, and load it
            [pth,nam] = spm_fileparts(deblank(anatT1));
            sn_filename  = fullfile(pth,[nam '_sn.mat']);
            if ~exist(sn_filename,'file')
                sn_filename  = fullfile(pth,['m' nam '_sn.mat']);
            end
            NIRS.Dt.ana.wT1 = load(sn_filename);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 1: coregistration: finding the affine transformation
            % relating the 2 coordinate systems
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
            wT1 = NIRS.Dt.ana.wT1;
            Q = (wT1.VG.mat/wT1.Affine)/wT1.VF.mat;
            % FIT ROM TO RMM POSITIONS, USING FIDUCIALS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(job.nasion_wMNI(:)) > 3
                job_MNI_fiducials = [job.nasion_wMNI(:,iSubj) job.AL_wMNI(:,iSubj) job.AR_wMNI(:,iSubj)];
            else
                job_MNI_fiducials = [job.nasion_wMNI job.AL_wMNI job.AR_wMNI];
            end
            % Positions of fiducial points
            if job.fiducial_MNI_choice
                %use specified coordinates in subject MNI coordinates
                NIRS.Cf.H.F.r.m.mm.p  = job_MNI_fiducials;
                Fp_rmm = NIRS.Cf.H.F.r.m.mm.p;
                temp_Fp_wmm = Q*[Fp_rmm; [1 1 1]];
                NIRS.Cf.H.F.w.m.mm.p = temp_Fp_wmm(1:3,:);
            else
                NIRS.Cf.H.F.w.m.mm.p  = job_MNI_fiducials;
                Fp_wmm = NIRS.Cf.H.F.w.m.mm.p;
                temp_Fp_rmm = Q\[Fp_wmm; [1 1 1]];
                NIRS.Cf.H.F.r.m.mm.p = temp_Fp_rmm(1:3,:);
            end
            y = NIRS.Cf.H.F.r.m.mm.p;
            x = NIRS.Cf.H.F.r.o.mm.p;
            
            % Compute optimal rigid transformation to match fiducial
            % positions from one coordinate system to the other
            [s R t] = abs_orientation(x,y); % y = s*R(x) + t
            estY = zeros(size(y));
            for Fi = 1:3 %assume 3 fiducials
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
            Sp_rom = NIRS.Cf.H.S.r.o.mm.p;
            Dp_rom = NIRS.Cf.H.D.r.o.mm.p;
            try
                Qp_rom = NIRS.Cf.H.Q.r.o.mm.p;
                Pp_rom = [Sp_rom Dp_rom Qp_rom];
            catch %No Q-type points
                Pp_rom = [Sp_rom Dp_rom];
            end
            NP = size(Pp_rom,2);
            NIRS.Cf.H.P.N = NP;
            
            Pp_rmm = zeros(size(Pp_rom));
            Pvoid = zeros(1,size(Pp_rom,2));
            
            for Pi = 1:NP
                % check for Void sources (no data)
                if Pp_rom(1,Pi) == 0 && Pp_rom(2,Pi) == 0 && Pp_rom(3,Pi) == 0
                    Pvoid(Pi) = 1;
                else
                    Pp_rmm(:,Pi) = s*R*Pp_rom(:,Pi) + t;
                end
            end;
            
            % Save MNI coordinates of optodes -- coregistration is now done
            NIRS.Cf.H.P.r.m.mm.p = Pp_rmm;
            NIRS.Cf.H.P.r.m.mm.fp = Pp_rmm; %fitted on skin
            
            % unnormalized -> normalized, for optodes
            Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm
            Ns = NIRS.Cf.H.S.N;
            % PROJECT OPTODE POSITIONS ON CORTEX SURFACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %render to subject
            if isfield(NIRS.Dt.ana,'T1seg') && use_fSeg
                fSeg = NIRS.Dt.ana.T1seg;
            else
                fSeg = [];
            end
            Nch0 = size(NIRS.Cf.H.C.id,2)/2;
            
            Pch_rmm = zeros(3,Nch0);
            Pch_wmm = zeros(4,Nch0); %different convention of keeping the 4th row
            for i=1:Nch0
                %indices of source and detector
                Si = NIRS.Cf.H.C.id(2,i);
                Di = NIRS.Cf.H.C.id(3,i)+Ns;
                Pch_rmm(:,i) = (Pp_rmm(:,Si)+Pp_rmm(:,Di))/2;
                Pch_wmm(:,i) = (Pp_wmm(:,Si)+Pp_wmm(:,Di))/2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Step 2: Cortex projection method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch cpm
                case 1 %PP
                    p_cutoff = 0.8;
                    %project sources and detectors onto cortex
                    Pp_c1_rmm = nirs_coreg_optodes_unnormalized(Pp_rmm,...
                        Pvoid,p_cutoff,fc1,fSeg,0,dir_coreg);
                    Pp_c1_wmm = zeros(4,NP);
                    
                    %Improve the projection by checking if positions have moved
                    %sufficiently -- if not, use the line to the center of the
                    %brain as the best direction to the cortex
                    %typical distance:
                    Pdist = pdist2(Pp_c1_rmm(1:3,:)',Pp_rmm');
                    Pdist = diag(Pdist);
                    Pdist_m = Pdist(logical(1-Pvoid));
                    Mdist = max(13,mean(Pdist_m)); %ensure to move by at least 13 mm toward the cortex
                    corr_factor = 1.4; %1.5; %to ensure that all optodes will be on cortex
                    %dcutoff = 10; %in millimeters
                    Zcenter_brain = 10;
                    for Pi=1:NP
                        if ~Pvoid(Pi)
                            cPs = Pp_rmm(:,Pi);
                            %cPc = Pp_c1_rmm(:,Pi);
                            %if Pdist(Pi) < dcutoff
                            dcent = pdist2(cPs',[0 0 Zcenter_brain]);
                            Pp_c1_rmm(:,Pi) = [cPs*(1-corr_factor*Mdist/dcent);1];
                            %end
                        end
                        %if pdist2(cPs,cP
                    end
                    
                    %project channels onto cortex
                    Pch_c1_rmm = nirs_coreg_optodes_unnormalized(Pch_rmm,...
                        Pvoid,p_cutoff,fc1,fSeg,1,dir_coreg);
                    %Pdist = pdist2(Pch_c1_rmm(1:3,:)',Pch_rmm');
                    %Pdist = diag(Pdist);
                    for Pi=1:Nch0
                        cPs = Pch_rmm(:,Pi);
                        %cPc = Pp_c1_rmm(:,Pi);
                        %if Pdist(Pi) < dcutoff
                        dcent = pdist2(cPs',[0 0 0]);
                        Pch_c1_rmm(:,Pi) = [cPs*(1-corr_factor*Mdist/dcent);1];
                        %end
                        %if pdist2(cPs,cP
                    end
                    Pch_c1_wmm = Q*Pch_c1_rmm;
                    for Pi = 1:NP
                        if ~Pvoid(Pi)
                            %inversion: unnormalized -> normalized
                            Pp_c1_wmm(:,Pi) = Q*Pp_c1_rmm(:,Pi);
                        end
                    end;
                case 2
                    Pp_c1_wmm = nirs_coreg_optodes(Pp_wmm,...
                        Pvoid,job.coreg_choice.coreg_c1,NIRS.Dt.ana.T1);
                    if isempty(Pp_c1_wmm)
                        disp('Coregistration of optodes onto wc1 failed. Now try to coregister onto template using projection_CS');
                        Pp_c1_wmm = projection_CS(Pp_wmm);
                    end
                    
                    Pp_c1_rmm = zeros(4,NP);
                    for Pi = 1:NP
                        if ~Pvoid(Pi)
                            %inversion: normalized -> unnormalized
                            Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi);
                        end
                    end;
                case 3
                    Pp_c1_wmm = projection_CS(Pp_wmm);%Using Korean template
                    Pch_c1_wmm = projection_CS(Pch_wmm);
                    Pp_c1_rmm = zeros(4,NP);
                    for Pi = 1:NP
                        if ~Pvoid(Pi)
                            %inversion: normalized -> unnormalized
                            Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi);
                        end
                    end
                    for Pi=1:Nch0
                        Pch_c1_rmm(:,Pi) = Q\Pch_c1_wmm(:,Pi);
                    end
            end
            
            Pp_c1_rmm = Pp_c1_rmm(1:3,:);
            
            % Save Normalized coordinates of optodes
            NIRS.Cf.H.P.w.m.mm.p = Pp_wmm;
            NIRS.Cf.H.P.w.m.mm.c1.p = Pp_c1_wmm;
            % Save MNI coordinates of optodes on cortex (c1)
            NIRS.Cf.H.P.r.m.mm.c1.p = Pp_c1_rmm;
            NIRS.Cf.H.P.void = Pvoid;
            
            % OPTODE POSITIONS ON SKIN SURFACE %
            %PP: very long, just skip -- looks like an infinite loop!
            oldClement_fitSkin = 0;
            if oldClement_fitSkin
                if ~isempty(fSeg)
                    jobe.Pp_rmm = Pp_rmm;
                    jobe.Pp_c1_rmm = Pp_c1_rmm;
                    jobe.NP = NP;
                    jobe.image_in = {fSeg};
                    jobe.lby = 'coreg';
                    %PP: I have not checked this function
                    out2 = nirs_fit_probe(jobe); %careful, this is not the same out as the
                    %out of the function with NIRS.mat!
                    Pfp_rmm = out2{1};
                    % from MNI real space (mm) to MNI voxel space
                    V_4fit = spm_vol(fSeg);
                    Pfp_rmv = [];
                    for i=1:NP
                        Pfp_rmv(:,i) = V_4fit.mat\[Pfp_rmm(:,i);1];
                    end
                    Pfp_rmv = Pfp_rmv(1:3,:);
                    % Save MNI and voxel coordinates of optodes, fitted on skin surface
                    NIRS.Cf.H.P.r.m.mm.fp = Pfp_rmm;
                    NIRS.Cf.H.P.r.m.vx.fp = Pfp_rmv;
                end
            end
            % Save all changes to NIRS structure
            save(newNIRSlocation,'NIRS');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 3: GENERATE TOPO DATA %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            V = spm_vol(fwT1);
            wT1_info.mat = V.mat;
            wT1_info.dim = V.dim;
            %let's use positions of optodes on cortex
            %loop over channel ids
            Nch0 = size(NIRS.Cf.H.C.id,2)/2;
            ch_MNIw_vx = zeros(4,Nch0);
            ch_MNI_vx = zeros(4,Nch0);
            src_MNIw_vx = zeros(4,Ns);
            src_MNI_vx = zeros(4,Ns);
            Nd = NIRS.Cf.H.D.N;
            det_MNIw_vx = zeros(4,Nd);
            det_MNI_vx = zeros(4,Nd);
            %same for skin
            ch_MNIw_vx_skin = zeros(4,Nch0);
            ch_MNI_vx_skin = zeros(4,Nch0);
            src_MNIw_vx_skin = zeros(4,Ns);
            src_MNI_vx_skin = zeros(4,Ns);
            det_MNIw_vx_skin = zeros(4,Nd);
            det_MNI_vx_skin = zeros(4,Nd);
            
            Vfc1 = spm_vol(fc1);
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
                posw = V.mat\(Q*posmm);
                ch_MNIw_vx(:,i) = posw;
                %For the skin:
                posmm_skin = [(Pp_rmm(:,Si)+Pp_rmm(:,Di))/2;1];
                pos_skin = NIRS.Dt.ana.wT1.VF.mat\posmm_skin;
                ch_MNI_vx_skin(:,i) = pos_skin;
                posw_skin = V.mat\(Q*posmm_skin);
                ch_MNIw_vx_skin(:,i) = posw_skin;
            end
            
            if coreg_projected_channels_on_cortex_rather_than_midpoint
                ch_MNI_vx = Vfc1.mat\Pch_c1_rmm;
                ch_MNIw_vx = V.mat\Pch_c1_wmm;
            else
                for i=1:Nch0
                    %indices of source and detector
                    Si = NIRS.Cf.H.C.id(2,i);
                    Di = NIRS.Cf.H.C.id(3,i)+Ns;
                    poswmm_t = (Pp_c1_wmm(:,Si)+Pp_c1_wmm(:,Di))/2;%Normalised mm coordinates
                    pos_wvx_t = V.mat\poswmm_t;%Normalised voxel coordinates
                    ch_MNI_vx(:,i) = Vfc1.mat\[(Pp_c1_rmm(:,Si)+Pp_c1_rmm(:,Di))/2;1]; %NIRS.Dt.ana.wT1.Affine\pos_wvx_t;%Un-normalised vx coordinates
                    ch_MNIw_vx(:,i) = pos_wvx_t;
                end
            end
            
            for i=1:Ns
                src_MNI_vx(:,i) = Vfc1.mat\[Pp_c1_rmm(:,i);1];
                src_MNIw_vx(:,i) = V.mat\Pp_c1_wmm(:,i);
                %skin
                src_MNI_vx_skin(:,i) = Vfc1.mat\[Pp_rmm(:,i);1];
                src_MNIw_vx_skin(:,i) = V.mat\Pp_wmm(:,i);
            end
            for i=1:Nd
                det_MNI_vx(:,i) = Vfc1.mat\[Pp_c1_rmm(:,i+Ns);1];
                det_MNIw_vx(:,i) = V.mat\Pp_c1_wmm(:,i+Ns);
                %skin
                det_MNI_vx_skin(:,i) = Vfc1.mat\[Pp_rmm(:,i+Ns);1];
                det_MNIw_vx_skin(:,i) = V.mat\Pp_wmm(:,i+Ns);
            end
            
            %PP: the function render_MNI_coordinates has not been checked
            rendered_MNI = render_MNI_coordinates_new(ch_MNIw_vx,...
                ch_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0,ERad);
            
            %Additional projections
            %channels on skin
            if OutputSkinFigs
                rendered_MNI_skin = render_MNI_coordinates_new(ch_MNIw_vx_skin,...
                    ch_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,0,fSeg,dir_coreg,1,ERad);
            end
            %sources on cortex
            rendered_MNI_src = render_MNI_coordinates_new(src_MNIw_vx,...
                src_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0,ERad);
            %detectors on cortex
            rendered_MNI_det = render_MNI_coordinates_new(det_MNIw_vx,...
                det_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0,ERad);
            if OutputSkinFigs
                %sources on skin
                rendered_MNI_src_skin = render_MNI_coordinates_new(src_MNIw_vx_skin,...
                    src_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,0,fSeg,dir_coreg,1,ERad);
                %detectors on skin
                rendered_MNI_det_skin = render_MNI_coordinates_new(det_MNIw_vx_skin,...
                    det_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,0,fSeg,dir_coreg,1,ERad);
            end
            rend_file = fullfile(dir_coreg,'TopoData.mat');
            save(rend_file, 'rendered_MNI');
            NIRS.Dt.ana.rend = rend_file;
            %save coordinates
            try
                if render_template
                    NIRS.Cf.H.C.w.m.vx.c1.p = ch_MNIw_vx;
                    NIRS.Cf.H.C.w.m.vx.fp = ch_MNIw_vx_skin;
                    NIRS.Cf.H.S.w.m.vx.c1.p = src_MNIw_vx;
                    NIRS.Cf.H.S.w.m.vx.fp = src_MNIw_vx_skin;
                    NIRS.Cf.H.D.w.m.vx.c1.p = det_MNIw_vx;
                    NIRS.Cf.H.D.w.m.vx.fp = det_MNIw_vx_skin;
                    NIRS.Cf.H.C.w.m.mm.c1.p = V.mat*ch_MNIw_vx;
                    NIRS.Cf.H.C.w.m.mm.fp = V.mat*ch_MNIw_vx_skin;
                    NIRS.Cf.H.S.w.m.mm.c1.p = V.mat*src_MNIw_vx;
                    NIRS.Cf.H.S.w.m.mm.fp = V.mat*src_MNIw_vx_skin;
                    NIRS.Cf.H.D.w.m.mm.c1.p = V.mat*det_MNIw_vx;
                    NIRS.Cf.H.D.w.m.mm.fp = V.mat*det_MNIw_vx_skin;
                end
                NIRS.Cf.H.C.r.m.vx.c1.p = ch_MNI_vx;
                NIRS.Cf.H.C.r.m.vx.fp = ch_MNI_vx_skin;
                NIRS.Cf.H.S.r.m.vx.c1.p = src_MNI_vx;
                NIRS.Cf.H.S.r.m.vx.fp = src_MNI_vx_skin;
                NIRS.Cf.H.D.r.m.vx.c1.p = det_MNI_vx;
                NIRS.Cf.H.D.r.m.vx.fp = det_MNI_vx_skin;
                NIRS.Cf.H.C.r.m.mm.c1.p = Vfc1.mat*ch_MNI_vx;
                NIRS.Cf.H.C.r.m.mm.fp = Vfc1.mat*ch_MNI_vx_skin;
                NIRS.Cf.H.S.r.m.mm.c1.p = Vfc1.mat*src_MNI_vx;
                NIRS.Cf.H.S.r.m.mm.fp = Vfc1.mat*src_MNI_vx_skin;
                NIRS.Cf.H.D.r.m.mm.c1.p = Vfc1.mat*det_MNI_vx;
                NIRS.Cf.H.D.r.m.mm.fp = Vfc1.mat*det_MNI_vx_skin;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 4: output various figures
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dir_extra_coreg = fullfile(dir_coreg,'extra_coreg');
            if ~exist(dir_extra_coreg,'dir'), mkdir(dir_extra_coreg); end
            nirs_brain_project_2d(NIRS,dir_coreg,rendered_MNI,[],'r','','',[],0);
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI,[],'r','','',[],1); %just copy the files
            if OutputSkinFigs
                nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_skin,[],'r','','skin',[],0);
            end
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_src,rendered_MNI_det,'b','g','SD',Pvoid,0);
            if OutputSkinFigs
                nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_src_skin,rendered_MNI_det_skin,'b','g','SD_skin',Pvoid,0);
            end
            NIRS.jobCoreg = job; %Required when coregistration needs to be redone in the GLM
            NIRS.flags.coregOK = 1;
            save(newNIRSlocation,'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Coregistration failed for the ' int2str(iSubj) 'th subject. for ' job.NIRSmat{iSubj,1}]);
    end
end
out.NIRSmat = job.NIRSmat;