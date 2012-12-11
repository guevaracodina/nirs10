function out = nirs_run_coreg_new(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
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
useKoreanCortexProjection = 1;
ForceReprocessSegmentation = 0;
coreg_projected_channels_on_cortex_rather_than_midpoint = 1;
use_fSeg = 1; %option to use Cl�ment's segmented image to extract c1; otherwise default to c1 image
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
                %bug when people select options in Cl�ment's segmentation 
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
            job_MNI_fiducials = [job.nasion_wMNI' job.AL_wMNI' job.AR_wMNI'];
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
            % positions to those of the normalized atlas
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
            
            % unnormalized -> normalized, for optodes
            Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm
            
            % PROJECT OPTODE POSITIONS ON CORTEX SURFACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~render_template
                %render to subject
                if isfield(NIRS.Dt.ana,'T1seg') && use_fSeg
                    fSeg = NIRS.Dt.ana.T1seg;
                else
                    fSeg = [];
                end
                p_cutoff = 0.8;
                %project sources and detectors onto cortex
                Pp_c1_rmm = nirs_coreg_optodes_unnormalized(Pp_rmm,...
                    Pvoid,p_cutoff,fc1,fSeg,0,dir_coreg);
                Pp_c1_wmm = zeros(4,NP);
                
                Nch0 = size(NIRS.Cf.H.C.id,2)/2;
                Ns = NIRS.Cf.H.S.N;
                Pch_rmm = zeros(3,Nch0);
                for i=1:Nch0
                    %indices of source and detector
                    Si = NIRS.Cf.H.C.id(2,i);
                    Di = NIRS.Cf.H.C.id(3,i)+Ns;
                    Pch_rmm(:,i) = (Pp_rmm(:,Si)+Pp_rmm(:,Di))/2;
                end
                %project channels onto cortex
                Pch_c1_rmm = nirs_coreg_optodes_unnormalized(Pch_rmm,...
                    Pvoid,p_cutoff,fc1,fSeg,1,dir_coreg);
                Pch_c1_wmm = Q*Pch_c1_rmm;
                for Pi = 1:NP
                    if ~Pvoid(Pi)
                        %inversion: unnormalized -> normalized
                        Pp_c1_wmm(:,Pi) = Q*Pp_c1_rmm(:,Pi);
                    end
                end;
            else
                if useKoreanCortexProjection
                    Pp_c1_wmm = projection_CS(Pp_wmm);%Using Korean template
                else
                    %This will no longer work because field
                    %job.coreg_choice.coreg_c1 doesn't exist anymore
                    Pp_c1_wmm = nirs_coreg_optodes(Pp_wmm,...
                        Pvoid,job.coreg_choice.coreg_c1,NIRS.Dt.ana.T1);
                    if isempty(Pp_c1_wmm)
                        disp('Coregistration of optodes onto wc1 failed. Now try to coregister onto template using projection_CS');
                        Pp_c1_wmm = projection_CS(Pp_wmm);
                    end
                end
                Pp_c1_rmm = zeros(4,NP);
                for Pi = 1:NP
                    if ~Pvoid(Pi)
                        %inversion: normalized -> unnormalized
                        Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi);
                    end
                end;
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
            
            % GENERATE TOPO DATA %
            %%%%%%%%%%%%%%%%%%%%%%
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
                ch_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0);
            
            %Additional projections
            %channels on skin
            rendered_MNI_skin = render_MNI_coordinates_new(ch_MNIw_vx_skin,...
                ch_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,1);
            %sources on cortex
            rendered_MNI_src = render_MNI_coordinates_new(src_MNIw_vx,...
                src_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0);
            %detectors on cortex
            rendered_MNI_det = render_MNI_coordinates_new(det_MNIw_vx,...
                det_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,0);
            
            %sources on skin
            rendered_MNI_src_skin = render_MNI_coordinates_new(src_MNIw_vx_skin,...
                src_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,1);
            %detectors on skin
            rendered_MNI_det_skin = render_MNI_coordinates_new(det_MNIw_vx_skin,...
                det_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,render_template,fSeg,dir_coreg,1);
            
            rend_file = fullfile(dir_coreg,'TopoData.mat');
            save(rend_file, 'rendered_MNI');
            NIRS.Dt.ana.rend = rend_file;
            dir_extra_coreg = fullfile(dir_coreg,'extra_coreg');
            if ~exist(dir_extra_coreg,'dir'), mkdir(dir_extra_coreg); end
            nirs_brain_project_2d(NIRS,dir_coreg,rendered_MNI,[],'r','','',[]);
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI,[],'r','','',[]); %just copy the files
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_skin,[],'r','','skin',[]);
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_src,rendered_MNI_det,'b','g','SD',Pvoid);
            nirs_brain_project_2d(NIRS,dir_extra_coreg,rendered_MNI_src_skin,rendered_MNI_det_skin,'b','g','SD_skin',Pvoid);           
            NIRS.flags.coregOK = 1;
        end
        save(newNIRSlocation,'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Coregistration failed for the ' int2str(iSubj) 'th subject. for ' job.NIRSmat{iSubj,1}]);
    end
end
out.NIRSmat = job.NIRSmat;