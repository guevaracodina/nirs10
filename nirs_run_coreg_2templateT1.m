function out = nirs_run_coreg_2templateT1(job)
% Clement

% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    try
        NIRS = [];
        load(job.NIRSmat{iSubj,1});
        
        % SPATIAL NORMALIZATION OF ANATOMICAL IMAGE %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                [DirSPM,~,~] = fileparts(which('spm'));
                NIRS.Dt.ana.tT1 = fullfile(DirSPM,'templates','T1.nii');
            end
        end
        
        [dirT1, fil, ext] = fileparts(NIRS.Dt.ana.T1);
        fwT1 = [dirT1 filesep 'w' fil ext(1:4)];
        if spm_existfile(fwT1)
            disp('Spatial normalisation already run - skipping');
        else
            
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
            
        end
        
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
        
        
        % FIT ROM TO RMM POSITIONS, USING FIDUCIES %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Positions of fiducial points
        if job.fid_in_subject_MNI
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
        disp(['Error Value for subject ' int2str(iSubj) ': ' num2str(errVal)]);
        NIRS.Dt.pro.errValofCoreg_mm2 = errVal;
        
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
        
        %%%[1] CORRIGER LES POSITIONS DE BRAINSIGHT SI NECESSAIRE
        %%% (2)enregistrer les sous dites positions
        NIRS.Cf.H.P.r.m.mm.p(1,2) = -83;
        NIRS.Cf.H.P.r.m.mm.p(1,7) = -82;
        NIRS.Cf.H.P.r.m.mm.p(1,8) = -86;
        NIRS.Cf.H.P.r.m.mm.p(1,9) = -88;
        
        %%%%%% A CHANGER
        save('Y:\cours\GBM8307\casqueSPM\NIRS.mat','NIRS');
        clear NIRS
        load('Y:\cours\GBM8307\casqueSPM\NIRS.mat');
        
        %%% (3)appliquer la transformation inverse Q
        Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
        % unnormalized -> normalized, for optodes
        Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm
        Pp_c1_wmm = projection_CS(Pp_wmm);
        
        NP = size(Pp_wmm,2);
        NIRS.Cf.H.P.N = NP;
        NIRS.Cf.H.P.w.m.mm.p = Pp_wmm;        % Save Normalized coordinates of optodes
        NIRS.Cf.H.P.w.m.mm.c1.p = Pp_c1_wmm;
        NIRS.Cf.H.P.void = Pvoid;
        
        %%% (2bis) correction des optodes en dehors du volume de T1_segmented
        %%% puisque tres peu d air.......
        %%%%%% A CHANGER
        NIRS.Cf.H.P.w.m.mm.p(2,1) = 89;
        save('Y:\cours\GBM8307\casqueSPM\NIRS.mat','NIRS')
        clear NIRS
        load('Y:\cours\GBM8307\casqueSPM\NIRS.mat');
        
        
        Pp_t1rmm = NIRS.Cf.H.P.w.m.mm.p;
        Pp_c1_t1rmm = NIRS.Cf.H.P.w.m.mm.c1.p;
        
        %%% (4) a ce stade j'ai les bonnes positions en wmm donc sur l'image
        %%% normalisee du sujet............................................. qui
        %%% est ''fsegT1_4fit'' et donc sur lq te;plqtem puisau4on est dans le
        %%% ;[e;e r2f2rentiel et on peut approximer que c'est les memes images
        %%% Le TRUC : c'est que les positions en mm sont les memes sur la T1
        %%% template ou sur la normalisee du sujet.........
        % fsegT1_4fit = 'D:\Users\Cl�ment\cours\GBM8307\projet_casque\00044_segmented_T1.nii';
        fsegT1_4fit = 'Y:\cours\GBM8307\casqueSPM\00044_segmented_T1.nii';
        
        jobF.NP = NP;
        jobF.Pp_rmm = Pp_t1rmm(1:3,:);
        jobF.Pp_c1_rmm = Pp_c1_t1rmm(1:3,:);
        jobF.image_in = {fsegT1_4fit};
        jobF.lby = 'coreg';
        outF = nirs_fit_probe(jobF);
        Pfp_t1rmm = outF{1};
        
        % from MNI real space (mm) to MNI voxel space
        V_4fit = spm_vol(fsegT1_4fit);
        Pfp_t1rmv = [];
        for i=1:NP
            Pfp_t1rmv(:,i) = V_4fit.mat\[Pfp_t1rmm(:,i);1];
        end
        Pfp_t1rmv = Pfp_t1rmv(1:3,:);
        
        % fake Cs.temp...
        NIRS.Cf.dev.wl = [690 830];
        %%%
        NIRS.Cs.temp.segR = {fsegT1_4fit}; %ROI from temp
        NIRS.Cs.temp.Pfp_roi_rmv = Pfp_t1rmv;
        NIRS.Cs.temp.Pfp_roi_rmm = Pfp_t1rmm;
        NIRS.Cs.temp.Pp_roi_rmm = Pp_t1rmm(1:3,:);
        NIRS.Cs.temp.Pp_roi_c1_rmm = Pp_c1_t1rmm(1:3,:);
        NIRS.Cs.temp.Pkpt = [1 2 3 4 5 6 7 8 9 10 11 12];
        NIRS.Cs.temp.NSkpt = 3;
        NIRS.Cs.temp.NDkpt = 9;
        
%         save('D:\Users\Cl�ment\cours\GBM8307\casqueSPM\NIRS.mat','NIRS');
        save('Y:\cours\GBM8307\casqueSPM\NIRS.mat','NIRS');
        
        % GENERATE TOPO DATA %
        %%%%%%%%%%%%%%%%%%%%%%
        
        load(job.NIRSmat{iSubj});
        if job.GenDataTopo
            [pth2,~,~] = fileparts(job.NIRSmat{iSubj,1});
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
                ch_MNI_vx = []; %%%% la notation serait Cp_rmv
                %number of sources
                Ns = NIRS.Cf.H.S.N;
                for i=1:(size(NIRS.Cf.H.C.id,2)/2)
                    %indices of source and detector
                    Si = NIRS.Cf.H.C.id(2,i);
                    Di = NIRS.Cf.H.C.id(3,i)+Ns;
                    pos = V.mat\(Q*[(Pp_c1_rmm(:,Si)+Pp_c1_rmm(:,Di))/2;1]);
                    
                    ch_MNI_vx = [ch_MNI_vx pos]; %%%% la notation serait Cp_rmv
                end
                [rend, rendered_MNI] = render_MNI_coordinates(ch_MNI_vx, wT1_info);%%%% la notation serait Cp_rmv
                for kk = 1:6
                    rendered_MNI{kk}.ren = rend{kk}.ren;
                end
                rend_file = fullfile(pth2,'TopoData.mat');
                save(rend_file, 'rendered_MNI');
                NIRS.Dt.ana.rend = rend_file;
                
                %Viewer from NIRS_SPM
                viewer_ON = 0;
                if viewer_ON
                    %rendered_MNI = varargin{1};
                    Nch = size(rendered_MNI{1}.rchn,1);
                    load Split
                    for kk=2:6
                        figure;
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
                    end
                end
            catch
                disp('Could not create TopoData.mat file');
            end
        end
        save(job.NIRSmat{iSubj},'NIRS');
    catch
        disp(['Coregistration failed for the ' int2str(iSubj) 'th subject.']);
    end
end
out.NIRSmat = job.NIRSmat;%job.NIRSmat{iSubj};


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