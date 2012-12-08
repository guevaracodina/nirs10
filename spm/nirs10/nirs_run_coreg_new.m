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
ForceReprocess = job.ForceReprocess;
% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
        if ~isfield(NIRS,'flags'), NIRS.flags = []; end
        job.NIRSmat{iSubj,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'coregOK') || job.force_redo)
            %Normalization of anatomical image
            anatT1 = NIRS.Dt.ana.T1;
            [dirT1, fil, ext] = fileparts(anatT1);
            fwT1 = [dirT1 filesep 'w' fil ext(1:4)];
            if spm_existfile(fwT1) && ~ForceReprocess
                disp(['Spatial normalisation already run in file ' fwT1 ' - skipping']);
            else
                nirs_batch_normalize(anatT1);
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
                    Pvoid(Pi,1) = 1;
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
            
            % OPTODE POSITIONS ON SKIN SURFACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(job.segT1_4fit{1,1})
                try
                    fsegT1_4fit = NIRS.Dt.ana.T1seg;
                catch
                    fsegT1_4fit = NIRS.Dt.ana.T1;
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
            %PP: I have not checked this function 
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
            
            % Save all changes to NIRS structure
            save(job.NIRSmat{iSubj},'NIRS');
            
            % GENERATE TOPO DATA %
            %%%%%%%%%%%%%%%%%%%%%%            
            [pth2,dummy,dummy2] = fileparts(job.NIRSmat{iSubj,1});
            [dir1,fil1,ext1] = fileparts(NIRS.Dt.ana.T1);
            fwT1 = fullfile(dir1,['w' fil1 ext1]);
            NIRS.Dt.ana.fwT1 = fwT1;
            try
                V = spm_vol(fwT1);
                wT1_info.mat = V.mat;
                wT1_info.dim = V.dim;
                %let's use positions of optodes on cortex
                %loop over channel ids
                Nch0 = size(NIRS.Cf.H.C.id,2)/2;
                ch_MNIw_vx = zeros(4,Nch0); 
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
                [rend, rendered_MNI] = render_MNI_coordinates(ch_MNIw_vx,...
                    ch_MNI_vx,wT1_info, NIRS.Dt.ana.wT1.VF,job.render_choice);
                %[rend, rendered_MNI] = render_MNI_coordinates(ch_MNIw_vx,ch_MNI_vx_skin,wT1_info, NIRS.Dt.ana.wT1.VF,job.render_choice);
                
                %%%% la notation serait Cp_rmv
                for kk = 1:6
                    rendered_MNI{kk}.ren = rend{kk}.ren;
                end
                rend_file = fullfile(pth2,'TopoData.mat');
                save(rend_file, 'rendered_MNI');
                NIRS.Dt.ana.rend = rend_file;
                
                Nch = size(rendered_MNI{1}.rchn,1);
                load Split
                for kk=1:6
                    fh0(kk) = figure;
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
                        [side_hemi spec_hemi] = nirs_get_brain_view(kk);
                        subj_str = '';
                        try
                            if isfield(NIRS.Dt.s,'subj_id')
                                subj_str = [NIRS.Dt.s.subj_id '_'];
                            end
                        end
                        
                        filen1 = fullfile(pth2,[subj_str spec_hemi '.fig']);
                        saveas(fh0(kk),filen1,'fig');
                        filen2 = fullfile(pth2,[subj_str spec_hemi '.tif']);
                        print(fh0(kk), '-dtiffn', filen2);
                        close(fh0(kk));                        
                    end
                end
            catch exception
                disp(exception.identifier);
                disp(exception.stack(1));
                disp('Could not create TopoData.mat file');
            end            
            NIRS.flags.coregOK = 1;            
        end
        save(job.NIRSmat{iSubj,1},'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Coregistration failed for the ' int2str(iSubj) 'th subject. for ' job.NIRSmat{iSubj,1}]);
    end
end
out.NIRSmat = job.NIRSmat;