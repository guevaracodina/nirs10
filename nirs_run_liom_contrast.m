function out = nirs_run_liom_contrast(job)
%Contrasts are generated in 3 different ways
%If no contrast is user-specified, automated contrasts are generated
%If contrasts are user-speficied, automated contrasts are not generated

%Contrasts need to be put in the SPM xCon structures
%If there is more than one session, contrasts can be run by session or
%over all sessions - the latter is the SPM standard way

%Views:
views_to_run = job.view;
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'

%Options
try
    GroupColorbars = job.GroupColorbars;
catch
    GroupColorbars = 0;
end
% try 
%     contrast_dir_name = job.contrast_dir_name;
% catch
%     contrast_dir_name = 'Contrast';
% end
try
    SmallFigures = job.SmallFigures;
catch
    SmallFigures = 0;
end
try
    write_neg_pos = job.write_neg_pos;
catch
    write_neg_pos = 0;
end
try
    GroupMultiSession = job.GroupMultiSession;
catch
    GroupMultiSession = 0;
end
try
    save_nifti_contrasts = job.save_nifti_contrasts;
catch
    save_nifti_contrasts = 1;
end
% try
%     Study_type = job.Study_type;
% catch
%     Study_type = 2;
% end

%When there are two or more sessions
%1: contrasts defined over more than 1 session are ignored
%0: contrasts defined over only 1 session are processed with full design
%matrix over all sessions
%2: both options 0 and 1 are run
try
    ProcessContrastsBySession = job.ProcessContrastsBySession;
catch
    ProcessContrastsBySession = 1;
end
try
    p_value = job.contrast_p_value;
catch
    p_value = 0.05;
end
%Colorbar - to allow specifying common min and max on all charts
try
    cbar.c_min = job.override_colorbar.colorbar_override.colorbar_min;
    cbar.c_max = job.override_colorbar.colorbar_override.colorbar_max;
    cbar.c_min2 = job.override_colorbar.colorbar_override.colorbar_min2;
    cbar.c_max2 = job.override_colorbar.colorbar_override.colorbar_max2;
    cbar.colorbar_override = 1;
catch
    cbar.colorbar_override = 0;
end
try
    output_unc = job.output_unc;
catch
    output_unc = 0;
end
%
try
    switch job.figures_visible
        case 1
            cbar.visible = 'on';
        case 0
            cbar.visible = 'off';
    end
catch
    cbar.visible = 'off';
end
%To generate inverted hemodynamic responses. Be careful when using this option
try
    GInv = job.GenerateInverted;
catch
    GInv = 1;
end
%Display feature - to group figures into subplots - careful, this may not
%group subplots in the correct order depending on which contrasts and other
%options you have selected
try
    GFIS = job.GroupFiguresIntoSubplots;
catch
    GFIS = 1;
end
%Booleans to choose which figures to write to disk, if any
switch job.contrast_figures
    case 0
        gen_fig = 0;
        gen_tiff = 0;
    case 1
        gen_fig = 1;
        gen_tiff = 1;
    case 2
        gen_fig = 1;
        gen_tiff = 0;
    case 3
        gen_fig = 0;
        gen_tiff = 1;
end
%SPM contrasts
if ~isfield(job,'consess')
    job.consess = [];
end
%automated contrasts if no user-specified contrast
if ~isempty(job.consess)
    automated_contrasts = 0;
else
    automated_contrasts = 1;
end
try 
    NonlinearEpilepsyOn = job.NonlinearEpilepsyOn;
catch
    NonlinearEpilepsyOn = 0;
end
%Gaussian spatial LPF
try
    radius = job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
    spatial_LPF = 1;
catch
    spatial_LPF = 0;
end
%New directory
try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        run_contrast_OK = 1;
        load(job.NIRSmat{Idx,1});
        NC = NIRS.Cf.H.C.N;
        %load topographic information (formerly known as preproc_info)
        try
            fname_ch = job.TopoData{1};
            if ~isempty(fname_ch)
                load(fname_ch);
                NIRS.Dt.ana.rend = fname_ch;
            else
                if isfield(NIRS.Dt.ana,'rend')
                    fname_ch = NIRS.Dt.ana.rend;
                    load(fname_ch);
                else
                    disp(['No TopoData structure for subject ' int2str(Idx) ...
                        '. Rerun coregistration or even first module and make sure that TopoData is generated.']);
                    run_contrast_OK = 0;
                end
            end
        catch  exception
            disp(exception.identifier);
            disp(exception.stack(1));
            run_contrast_OK = 0;
        end
        %load SPM - first GLM - might want to generalize
        dir1 = NIRS.SPM{1};
        load(fullfile(dir1,'SPM.mat'));
        %load TOPO (topographic maps) if already (partially) generated
        ftopo = fullfile(dir1,'TOPO.mat');
        TOPO = [];
        
        if NewDirCopyNIRS
            [dirN fil1 ext1] =fileparts(job.NIRSmat{Idx,1});
            dir2 = [dirN filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            job.NIRSmat{Idx,1} = newNIRSlocation;
            [dir fil1 ext1] = fileparts(ftopo);
            ftopo = fullfile(dir2, [fil1 ext1]);
            dir1 = dir2;
        else
            if exist(ftopo,'file'), load(ftopo); end
        end
        %Fill Z Structure, for passing most generic data
        Z = [];
        Z.gen_fig = gen_fig;
        Z.gen_tiff = gen_tiff;
        Z.p_value = p_value;
        Z.GroupColorbars = GroupColorbars;
        Z.dir1 = dir1;
        Z.cbar = cbar;
        Z.GInv = GInv;
        Z.GFIS = GFIS;
        Z.output_unc = output_unc;
        Z.SmallFigures = SmallFigures;
        Z.write_neg_pos = write_neg_pos;
        Z.save_nifti_contrasts = save_nifti_contrasts;
        %Z.contrast_dir_name = contrast_dir_name;
        %Handles for figures
        Pt = []; %positive responses, tube
        Pu = []; %positive responses, uncorrected
        Nt = []; %negative, tube
        Nu = []; %negative, uncorrected
        Ct = []; %positive and negative combined, tube
        Cu = []; %uncorrected
        %get contrasts
        [SPM xCon SSxCon] = nirs_get_contrasts(SPM,job,automated_contrasts,0,NonlinearEpilepsyOn);
        if run_contrast_OK
            %Big loop over views
            for v1=1:size(views_to_run,2)
                try
                    brain_view = views_to_run(v1);
                    [side_hemi spec_hemi] = nirs_get_brain_view(brain_view);
                    %Structure for passing GLM and interpolation data
                    clear W
                    % channel information
                    rchn = rendered_MNI{side_hemi}.rchn;
                    cchn = rendered_MNI{side_hemi}.cchn;
                    %rendering surface
                    brain = rendered_MNI{side_hemi}.ren;
                    %Fill W structure, less generic data than Z structure
                    W.brain = brain * 0.5;
                    W.s1 = size(brain, 1);
                    W.s2 = size(brain, 2);
                    %find channels which are visible from this projection view
                    W.index_ch = find(rchn ~= -1);
                    W.spatial_LPF = spatial_LPF;
                    if spatial_LPF
                        W.radius = radius;
                    end
                    if isempty(W.index_ch)
                        TOPO.v{side_hemi}.Warning = 'No channel found for this view';
                        disp(['No channel for view ' int2str(v1) ': Probable coregistration problem. Skipping this view']);
                    else
                        %split into HbO and HbR interpolations
                        wl = NIRS.Cf.dev.wl;
                        %for NIRS acquisition with 2 wavelengths only
                        
                        % MD: it seems to me that the channels corresponding to
                        % HbO/HbR should be fixed after the module
                        % nirs_run_ODtoHbOHbR ? Not sure this gives the right
                        % order....?!?
                        
                        %                         if wl(1) > 750
                        %                             %first wavelength is "HbO-like"
                        %                             W.ch_HbO = W.index_ch;
                        %                             W.ch_HbR = NC/2 + W.index_ch;
                        %                         else
                        %                             W.ch_HbR = W.index_ch;
                        %                             W.ch_HbO = NC/2 + W.index_ch;
                        %                         end
                        W.ch_HbO = W.index_ch;
                        W.ch_HbR = NC/2 + W.index_ch;
                        
                        try
                            W.ch_HbT = NC+W.index_ch;
                        catch exception
                            disp(exception.identifier);
                            disp(exception.stack(1));
                        end
                        %nch = length(index_ch);
                        W.rchn = rchn(W.index_ch);
                        W.cchn = cchn(W.index_ch);
                        W.spec_hemi = spec_hemi;
                        W.side_hemi = side_hemi;
                        
                        %                     switch Study_type
                        %                         case 0 %Single subject
                        %                             ProcessContrastsBySession = 2;
                        %                             GroupMultiSession = 0;
                        %                         case 1 %Group multi session
                        %                             ProcessContrastsBySession =
                        %                             GroupMultiSession = 1;
                        %                         case 2 %Group single session
                        %                             ProcessContrastsBySession =
                        %                             GroupMultiSession = 0;
                        %                     end
                        
                        if  ~(ProcessContrastsBySession == 1) %case 0 or 2
                            if ~GroupMultiSession
                                %REMOVE contrasts of wrong length
                                nC = length(xCon);
                                nCon = [];
                                for j1=1:nC
                                    tc = xCon(j1).c;
                                    if size(tc,1) == length(SPM.xXn)*nr
                                        %keep this contrast
                                        if isempty(nCon)
                                            nCon = xCon(j1);
                                        else
                                            nCon(end+1) = xCon(j1);
                                        end
                                    end
                                end
                            else
                                nCon = xCon;
                            end
                            %group of sessions
                            if GFIS
                                Pt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                    num2str(p_value) '_Pos'],'NumberTitle','off');
                                %subplot(fh0Pt,nC,3,1);
                                Pu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                    num2str(p_value) '_Pos'],'NumberTitle','off');
                                %subplot(fh0Pu,nC,3,1);
                                if GInv
                                    Nt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                        num2str(p_value) '_Neg'],'NumberTitle','off');
                                    %subplot(fh0Nt,nC,3,1);
                                    Nu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                        num2str(p_value) '_Neg'],'NumberTitle','off');
                                    Ct = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                        num2str(p_value)],'NumberTitle','off');
                                    Cu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                        num2str(p_value)],'NumberTitle','off');
                                end
                            end
                            beta_tmpO = [];
                            beta_tmpR = [];
                            beta_tmpT = [];
                            for f1 = 1:length(SPM.xXn)
                                beta_tmpO = [beta_tmpO; SPM.xXn{f1}.beta(:, W.ch_HbO)];
                                beta_tmpR = [beta_tmpR; SPM.xXn{f1}.beta(:, W.ch_HbR)];
                                try
                                    beta_tmpT = [beta_tmpT; SPM.xXn{f1}.beta(:, W.ch_HbT)];
                                end
                            end
                            W.var = SPM.xX.var; %careful, var can be a Matlab function,
                            %but instead we want W.var
                            W.beta_HbO = beta_tmpO(:); %taken as one vector
                            W.beta_HbR = beta_tmpR(:); %taken as one vector
                            W.mtx_var_HbO = diag(W.var(W.ch_HbO));
                            W.mtx_var_HbR = diag(W.var(W.ch_HbR));
                            try
                                W.beta_HbT = beta_tmpT(:); %taken as one vector
                                W.mtx_var_HbT = diag(W.var(W.ch_HbT));
                            catch exception
                                disp(exception.identifier);
                                disp(exception.stack(1));
                            end
                            
                            W.corr_beta = SPM.xX.corr_beta;
                            [TOPO] = constrasts_core(Z,W,TOPO,SPM.xX,nCon,0,Pt,Pu,Nt,Nu,Ct,Cu);
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if ~(ProcessContrastsBySession == 0) %case 1 or 2
                            %by session
                            %                     %REMOVE contrasts defined over more than one session
                            %                     nC = length(xCon);
                            %                     nCon = [];
                            %                     for j1=1:nC
                            %                         tc = xCon(j1).c;
                            %                         rk = 0;
                            %                         for f1=1:length(SPM.xXn)
                            %                             si = 1+nr*(f1-1); ei = nr*f1;
                            %                             if rank(tc(si:ei,:)) > 0
                            %                                 rk = rk+1;
                            %                             end
                            %                         end
                            %                         if rk == 1
                            %                             if isempty(nCon)
                            %                                 nCon = xCon(j1);
                            %                             else
                            %                                 nCon(end+1) = xCon(j1);
                            %                             end
                            %                         end
                            %                     end
                            %xCon = SSxCon;
                            %loop over sessions
                            for f1=1:length(SPM.xXn)
                                if GFIS
                                    Pt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                        num2str(p_value) '_S' int2str(f1) '_Pos'],'NumberTitle','off');
                                    %subplot(fh0Pt,nC,3,1);
                                    Pu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                        num2str(p_value) '_S' int2str(f1) '_Pos'],'NumberTitle','off');
                                    %subplot(fh0Pu,nC,3,1);
                                    if GInv
                                        Nt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                            num2str(p_value) '_S' int2str(f1) '_Neg'],'NumberTitle','off');
                                        %subplot(fh0Nt,nC,3,1);
                                        Nu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                            num2str(p_value) '_S' int2str(f1) '_Neg'],'NumberTitle','off');
                                        %Combined figures
                                        Ct = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                            num2str(p_value) '_S' int2str(f1)],'NumberTitle','off');
                                        Cu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                            num2str(p_value) '_S' int2str(f1)],'NumberTitle','off');
                                    end
                                end
                                try
                                    %for NIRS_SPM method
                                    W.var = SPM.xXn{f1}.ResSS./SPM.xXn{f1}.trRV;
                                    %covariance of beta estimates
                                    W.corr_beta = SPM.xXn{f1}.Bcov;
                                catch exception
                                    disp(exception.identifier);
                                    disp(exception.stack(1));
                                    %for WLS and BGLM methods
                                    W.corr_beta = SPM.xXn{f1}.Bvar;
                                    %will not work as we don't have var = ResSS/trRV
                                end
                                %GLM estimates - which beta though???
                                beta_tmp = SPM.xXn{f1}.beta(:, W.ch_HbO);
                                W.beta_HbO = beta_tmp(:); %taken as one vector
                                beta_tmp = SPM.xXn{f1}.beta(:, W.ch_HbR);
                                W.beta_HbR = beta_tmp(:); %taken as one vector
                                W.mtx_var_HbO = diag(W.var(W.ch_HbO));
                                W.mtx_var_HbR = diag(W.var(W.ch_HbR));
                                try
                                    beta_tmp = SPM.xXn{f1}.beta(:, W.ch_HbT);
                                    W.beta_HbT = beta_tmp(:); %taken as one vector
                                    W.mtx_var_HbT = diag(W.var(W.ch_HbT));
                                catch exception
                                    disp(exception.identifier);
                                    disp(exception.stack(1));
                                end
                                %special case when number of regressors
                                %varies between sessions
                                if automated_contrasts
                                     %get contrasts
                                    [SPM xCon SSxCon] = nirs_get_contrasts(SPM,job,automated_contrasts,f1,NonlinearEpilepsyOn);
                                end
                                [TOPO] = constrasts_core(Z,W,TOPO,SPM.xXn{f1},SSxCon,f1,Pt,Pu,Nt,Nu,Ct,Cu);
                                TOPO.SSxConS{f1} = SSxCon; %store contrasts by session
                            end %end for f1
                        end 
                    end
                    TOPO.v{side_hemi}.s1 = W.s1; %sizes of topographic projection
                    TOPO.v{side_hemi}.s2 = W.s2;
                    TOPO.v{side_hemi}.view = spec_hemi; %%% view of the brain
                catch exception
                    disp(exception.identifier);
                    disp(exception.stack(1));
                    disp(['Could not create contrasts for view ' spec_hemi ' for subject ' int2str(Idx)]);
                end
            end %end for v1
        end
        %TOPO.SPM = SPM; %save modified SPM - too big - not required
        
        try TOPO.xX = SPM.xX; end
        TOPO.xCon = xCon; %would not work if new contrasts are later added
        TOPO.SSxCon = SSxCon;
        save(ftopo,'TOPO','-v7.3'); %file can be large - really?
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not create contrasts for subject' int2str(Idx)]);
    end
    NIRS.TOPO = ftopo;
    if NewDirCopyNIRS
        save(newNIRSlocation,'NIRS');
    else
        save(job.NIRSmat{Idx,1},'NIRS');
    end
end
out.NIRSmat = job.NIRSmat;
end

function Q = interpolation_kernel(corr_beta,nch,mtx_var,W)
s1 = W.s1;
s2 = W.s2;
rchn = W.rchn;
cchn = W.cchn;
%identity over remaining channels
mtx_eye = eye(nch);
cov_beta = kron(mtx_var, corr_beta);
[U, S, V] = svd(cov_beta);
cov_beta_r = U*(S.^(0.5))*V';
%mesh of topographically projected brain size; note s1 and s2 inverted
[x, y] = meshgrid(1:s2, 1:s1);

B = zeros(s1, s2, nch);
Bx = zeros(s1, s2, nch);
By = zeros(s1, s2, nch);
for kk = 1:nch
    %a grid for channel positions?
    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic'); %could try nearest?
    if  kk == 1
        %mask the NaN
        mask = isnan(grid_eye);
        mask = 1- mask;
        [rmask{1}, cmask{1}] = find(mask == 1);
        index_mask0 = find(mask == 0);
    end
    grid_eye(index_mask0) = 0;
    [temp_Bx, temp_By] = gradient(grid_eye);
    temp_Bx(index_mask0) = 0;
    temp_By(index_mask0) = 0;
    
    B(:,:, kk) = grid_eye;
    Bx(:,:,kk) = temp_Bx;
    By(:,:,kk) = temp_By;
end
Q.cov_beta_r =cov_beta_r;
Q.B = B;
Q.Bx = Bx;
Q.By = By;
Q.rmask = rmask;
Q.cmask = cmask;
end

function [C xCon] = loop_contrasts(xCon,beta,corr_beta,mtx_var,Q,xX,W)
cov_beta_r = Q.cov_beta_r;
B = Q.B;
Bx = Q.Bx;
By = Q.By;
rmask = Q.rmask;
cmask = Q.cmask;
s1 = W.s1;
s2 = W.s2;
rmv = rmask{1};
cmv = cmask{1};
nC = size(xCon,2);
nm = length(rmv);
%preallocate
sum_kappa = zeros(nC,1);
kappa = zeros(nC,s1,s2);
c_interp_beta = zeros(nC,s1,s2); %analog of con0001.nii in SPM
c_cov_interp_beta = zeros(nC,s1,s2);
c_interp_T = zeros(nC,s1,s2); %T-stat
c_interp_ess = zeros(nC,s1,s2); %analog of ess0001.nii in SPM
c_interp_ess0 = zeros(s1,s2); %extra-sum of squares of full model (ResMS)
c_interp_F = zeros(nC,s1,s2); %F-stat
sz_xCon  = size(xCon(1).c,1);
%identity matrix of size number of regressors
tmp = eye(sz_xCon);
if isfield(xCon(1),'STAT')
    xConStructOK = 1;
else
    xConStructOK = 0;
end

for kk = 1:nm
    %this is different for HbO and HbR
    B2(:,1) = B(rmv(kk), cmv(kk), :);
    B2x(:,1) = Bx(rmv(kk), cmv(kk), :);
    B2y(:,1) = By(rmv(kk), cmv(kk), :);
    B3 = kron(B2, tmp);
    B3x = kron(B2x, tmp);
    B3y = kron(B2y, tmp);
    B3t = kron(B2', tmp);
    d = (B2'*mtx_var*B2); %ResSS/TrRV
    c_interp_ess0(rmv(kk), cmv(kk)) = d; %ResSS/TrRV -- recall normalization of ResSS by TrRV is included here (while it is a scale factor in the SPM nifti)
    for c1 = 1:nC
        c = xCon(c1).c;
        
        if (xConStructOK && xCon(c1).STAT == 'T') || ...
                ~xConStructOK %then must be a T-stat
            
            %this is the same for HbO and HbR
            c_corr_beta = c' * corr_beta * c;
            %this is different for HbO and HbR
            P = cov_beta_r * (B3 * c);
            Px = cov_beta_r * (B3x * c);
            Py = cov_beta_r * (B3y * c);
            tmp_1 = P'*P; tmp_2 = tmp_1^(-1/2); tmp_3 = tmp_2^3;
            u_derx = Px*tmp_2 - (P*(P'*Px))*tmp_3;
            u_dery = Py*tmp_2 - (P*(P'*Py))*tmp_3;
            %For each contrast and each channel, we get kappa, c_interp_beta
            %and c_cov_interp_beta
            
            %Positive contrasts
            kappa(c1,rmv(kk), cmv(kk)) =  ...
                sqrt(abs(det([u_derx'*u_derx u_derx'*u_dery; ...
                u_dery'*u_derx u_dery'*u_dery])));
            
            c_interp_beta(c1,rmv(kk), cmv(kk)) = (c' * B3t) * beta;
            
            c_cov_interp_beta(c1,rmv(kk), cmv(kk)) = (B2'*mtx_var*B2) * c_corr_beta;
            c_interp_T(c1,rmv(kk), cmv(kk)) = c_interp_beta(c1,rmv(kk), cmv(kk))/(c_cov_interp_beta(c1,rmv(kk), cmv(kk)))^0.5;
        else
            %This is an F-stat; if xCon(c1).STAT == 'F'
            %Do a GLM for the reduced model for each F contrast
            hsqr = xCon(c1).h * B3t;
            %Numerator of F-test
            c_ResSS =  (beta'*hsqr')*(hsqr*beta);
            %Interpolate
            c_interp_F(c1,rmv(kk), cmv(kk)) = c_ResSS/(d*xCon(c1).trRV);
            c_interp_ess(c1,rmv(kk), cmv(kk)) = c_ResSS/xCon(c1).trRV; %note that we also normalize here by eidf
        end
    end
end
for c1 = 1:nC %not for F contrasts...
    tm = kappa(c1,:,:);
    sum_kappa(c1) = sum(tm(:));
end
if W.spatial_LPF %does not work properly -- do not use
    K.k1 = s1;
    K.k2 = s2;
    K.radius = W.radius;
    K = spatial_LPF('set',K);
    c_interp_beta = spatial_LPF('lpf',K,c_interp_beta);
    c_cov_interp_beta = spatial_LPF('lpf',K,c_cov_interp_beta);
    %should interpolate all the rest too...
end
C.sum_kappa = sum_kappa;
C.c_interp_beta = c_interp_beta;
C.c_cov_interp_beta = c_cov_interp_beta;
C.c_interp_T = c_interp_T;
C.c_interp_F = c_interp_F;
C.c_interp_ess0 = c_interp_ess0;
C.c_interp_ess = c_interp_ess;
end

function [Pt,Pu,Nt,Nu,Ct,Cu] = interpolated_maps(Z,W,C,xCon,f1,erdf,hb,Pt,Pu,Nt,Nu,Ct,Cu)
%
sum_kappa = C.sum_kappa;
c_interp_beta = C.c_interp_beta;
c_cov_interp_beta = C.c_cov_interp_beta;
c_interp_T = C.c_interp_T;
c_interp_F = C.c_interp_F;
c_interp_ess = C.c_interp_ess;
c_interp_ess0 = C.c_interp_ess0;
s1 = W.s1;
s2 = W.s2;
spec_hemi = W.spec_hemi;
p_value = Z.p_value;
GInv = Z.GInv;
GFIS = Z.GFIS;

nC = size(xCon,2);
if isfield(xCon(1),'STAT')
    xConStructOK = 1;
else
    xConStructOK = 0;
end
%Structure F
load Split
F.split = split;
F.pathn = Z.dir1;
F.erdf = erdf;
if f1 > 0
    strf = ['_S' int2str(f1)];
    strf_fig = [' S' int2str(f1)];
else
    strf = '';
    strf_fig = '';
end
filestr = [num2str(p_value) '_' spec_hemi '_' hb strf];
filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb strf_fig];

%CF: copy figure structure
CF.GInv = GInv;
CF.split = split;
CF.nC = nC;

%loop over contrasts
try
    for c1=1:nC
        if (xConStructOK && xCon(c1).STAT == 'T') || ...
                ~xConStructOK %then must be a T-stat
            impose_T_bound = 1;
            index_mask = find(squeeze(c_cov_interp_beta(c1,:,:)) ~= 0);
            T_map = zeros(s1, s2);
            if impose_T_bound
                %put a bound on low variance -- as in SPM (see ResMS bound)
                tmp = squeeze(c_cov_interp_beta(c1,index_mask));
                bound = (1e-3)^0.5*max(tmp(isfinite(tmp)));
                tmp2 = bound*ones(size(tmp));
                tmp2(tmp > bound) = tmp(tmp > bound);
                T_map(index_mask) = squeeze(c_interp_beta(c1,index_mask))./ ...
                    sqrt(tmp2); %c_interp_T(c1,index_mask);
            else
                T_map(index_mask) = squeeze(c_interp_beta(c1,index_mask))./ ...
                    sqrt(squeeze(c_cov_interp_beta(c1,index_mask))); 
            end
            tstr = 'T';
        else
            %F-stats
            index_mask = find(squeeze(c_interp_F(c1,:,:)) ~= 0);
            %still use variable T_map, though it is now an F_map
            T_map = zeros(s1, s2);
            impose_F_bound = 1;
            if impose_F_bound
                tmp = squeeze(c_interp_ess0(index_mask));%should calculate just once
                bound = 1e-3*max(tmp(isfinite(tmp)));
                tmp2 = bound*ones(size(tmp));
                tmp2(tmp > bound) = tmp(tmp > bound);
                T_map(index_mask) = squeeze(c_interp_ess(c1,index_mask))./tmp2';                    
            else
                T_map(index_mask) = squeeze(c_interp_F(c1,index_mask));
            end
            tstr = 'F';
        end
        %Positive responses
        %names
        F.contrast_info = [filestr '_Pos' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Pos' xCon(c1).name];
        F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
        F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
        
        F.T_map = T_map;
        F.eidf = xCon(c1).eidf;
        F.sum_kappa = sum_kappa(c1);
        F.tstr = tstr;
        F.hb = hb;
        %contrast - what will be saved as nifti if requested
        if strcmp(tstr,'T')
            %F.con = zeros(s1,s2);
            %F.con = squeeze(c_interp_beta(c1,index_mask));
            F.con = squeeze(c_interp_beta(c1,:,:));
            F.ess = [];
        else
            F.con = [];
            %F.ess = zeros(s1,s2);
            %F.ess = squeeze(c_interp_F(c1,index_mask));
            F.ess = squeeze(c_interp_ess(c1,:,:));
        end
        %     if c1 == 3
        %         a=1;
        %     end
        %uncorrected - if we generate inverted responses, or HbO or HbT
        if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z); %DF: draw figure structure: handles to figure, axes and colorbar
                %copy figure structure
                if GFIS, [Pu,Nu,Cu] = nirs_copy_figure(Pu,Nu,Cu,DF,CF,c1,hb,1,tstr); end
            end
            %tube formula corrected for Tstats - or Bonferroni for Fstats
            %if tstr == 'T' %only for Tstats for now
            DF = nirs_draw_figure(3,F,W,Z);
            if GFIS, [Pt,Nt,Ct] = nirs_copy_figure(Pt,Nt,Ct,DF,CF,c1,hb,1,tstr); end
            %end
        end
        if  tstr == 'T' %for F-stat, do not invert the map - always positive
            %repeat for negative contrasts
            F.T_map = -T_map; %only used for generating Neg contrast figures
        end
        
        F.contrast_info = [filestr '_Neg' xCon(c1).name];
        F.contrast_info_for_fig = [filestr_fig ' Neg' xCon(c1).name];
        %contrast - what will be saved as nifti if requested
        if strcmp(tstr,'T')
            F.con = squeeze(c_interp_beta(c1,:,:));
            F.ess = [];
        else
            F.con = [];
            F.ess = squeeze(c_interp_ess(c1,:,:));
        end
        if GInv || strcmp(hb,'HbR')
            if Z.output_unc
                DF = nirs_draw_figure(2,F,W,Z);
                if GFIS, [Pu,Nu,Cu] = nirs_copy_figure(Pu,Nu,Cu,DF,CF,c1,hb,0,tstr); end
            end
            DF = nirs_draw_figure(3,F,W,Z);
            if GFIS, [Pt,Nt,Ct] = nirs_copy_figure(Pt,Nt,Ct,DF,CF,c1,hb,0,tstr); end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end


function [TOPO] = constrasts_core(Z,W,TOPO,xX,xCon,f1,Pt,Pu,Nt,Nu,Ct,Cu)
%Note that even for contrasts pertaining to individual
%sessions, the number of degrees of freedom it taken
%to be that for all the sessions - this is what is
%done in SPM - this looks like an error, however, it
%should be small under situations of interest, as the
%number of dof in each session is fairly large, even
%more so for fNIRS than for fMRI

side_hemi = W.side_hemi;
GFIS = Z.GFIS;
GInv = Z.GInv;
W.erdf = xX.erdf;
gen_fig = Z.gen_fig;
gen_tiff = Z.gen_tiff;
nC = size(xCon,2);
try
    wX = xX.xKXs.X;
catch
    wX = xX.xKXs; %quick fix for GroupMultiSessions - more properly, full xKXs
    %structure needs to be generated at the subject level from the
    %multisessions, in order to treat F contrasts
end
%Complete xCon structure for F contrasts (fields h and trRV)
for c1 = 1:nC
    if xCon(c1).STAT == 'F'
        xCon(c1).h = spm_FcUtil('Hsqr',xCon(c1), wX); %Need filtered design matrix
        switch xX.K.LParam.type
            case {'hrf', 'Gaussian'}
                S = xX.K.KL;
            case 'none'
                S = speye(nScan);
        end
        switch xX.K.HParam.type
            case 'DCT'
                S = S - xX.K.X0 * (xX.K.X0' * S);
                %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
        end
        %Calculate modified nubmer of degrees of freedom due to filtering
        %for the sum of squares difference between the full and reduced models
        trRV2 = approx_trRV(xX.xKXs.X,xX.pKX,S,xCon(c1).c);
        xCon(c1).trRV = trRV2;
    end
end

%HbO
try
    %Note that var(ch_HbO) depends on HbO vs HbR
    Q = interpolation_kernel(W.corr_beta,length(W.ch_HbO),W.mtx_var_HbO,W);
    [C xCon] = loop_contrasts(xCon,W.beta_HbO,W.corr_beta,W.mtx_var_HbO,Q,xX,W);
    if gen_fig || gen_tiff
        [Pt,Pu,Nt,Nu,Ct,Cu] = interpolated_maps(Z,W,C,xCon,f1,W.erdf,'HbO',Pt,Pu,Nt,Nu,Ct,Cu);
    end
    %fill TOPO
    TOPO = fill_TOPO(TOPO,C,side_hemi,f1,'HbO');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

%HbR
try
    Q = interpolation_kernel(W.corr_beta,length(W.ch_HbR),W.mtx_var_HbR,W);
    
    [C xCon] = loop_contrasts(xCon,W.beta_HbR,W.corr_beta,W.mtx_var_HbR,Q,xX,W);
    
    if gen_fig || gen_tiff
        [Pt,Pu,Nt,Nu,Ct,Cu] = interpolated_maps(Z,W,C,xCon,f1,W.erdf,'HbR',Pt,Pu,Nt,Nu,Ct,Cu);
    end
    TOPO = fill_TOPO(TOPO,C,side_hemi,f1,'HbR');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

try
    %HbT
    Q = interpolation_kernel(W.corr_beta,length(W.ch_HbT),W.mtx_var_HbT,W);
    
    [C xCon] = loop_contrasts(xCon,W.beta_HbT,W.corr_beta,W.mtx_var_HbT,Q,xX,W);
    
    if gen_fig || gen_tiff
        [Pt,Pu,Nt,Nu,Ct,Cu] = interpolated_maps(Z,W,C,xCon,f1,W.erdf,'HbT',Pt,Pu,Nt,Nu,Ct,Cu);
    end
    TOPO = fill_TOPO(TOPO,C,side_hemi,f1,'HbT');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

try %works for both group and single sessions
    if GFIS
        if Z.write_neg_pos || ~GInv
            save_assembled_figures(Z,W,Pt,'Pos','tube',f1);
            if Z.output_unc
                save_assembled_figures(Z,W,Pu,'Pos','unc',f1);
            else
                try close(Pu); end
            end
        else
            try close(Pt); end
            try close(Pu); end
        end
        if GInv
            if Z.write_neg_pos
                save_assembled_figures(Z,W,Nt,'Neg','tube',f1);
                if Z.output_unc
                    save_assembled_figures(Z,W,Nu,'Neg','unc',f1);
                else
                    try close(Nu); end
                end
            else
                try close(Nt); end
                try close(Nu); end
            end
            save_assembled_figures(Z,W,Ct,'','tube',f1);
            if Z.output_unc
                save_assembled_figures(Z,W,Cu,'','unc',f1);
            else
                try close(Cu); end
            end
        end
    end
catch exception
    disp(exception);
end
end

function inv_hb = inv_get_chromophore(h1)
switch h1
    case 'HbO'
        inv_hb = 1;
    case 'HbR'
        inv_hb = 2;
    case 'HbT'
        inv_hb = 3;
end
end

function TOPO = fill_TOPO(TOPO,C,side_hemi,f1,hb)
%Fills TOPO structure with interpolated beta, cov_beta, and F-contrast
hbi = inv_get_chromophore(hb);
%CONVENTION:
%'g': group of sessions
%'group': group analysis over subjects (see liom_group)
%'s': individual sessions
if f1 == 0
    TOPO.v{side_hemi}.g.hb{hbi}.sum_kappa = C.sum_kappa;
    TOPO.v{side_hemi}.g.hb{hbi}.c_interp_beta = C.c_interp_beta;
    TOPO.v{side_hemi}.g.hb{hbi}.c_cov_interp_beta = C.c_cov_interp_beta;
    TOPO.v{side_hemi}.g.hb{hbi}.c_interp_F = C.c_interp_F;
    TOPO.v{side_hemi}.g.hb{hbi}.c_interp_T = C.c_interp_T;
    TOPO.v{side_hemi}.g.hb{hbi}.c_interp_ess = C.c_interp_ess;
    TOPO.v{side_hemi}.g.hb{hbi}.c_interp_ess0 = C.c_interp_ess0;
    TOPO.v{side_hemi}.g.hb{hbi}.hb = hb;
else
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.sum_kappa = C.sum_kappa;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_interp_beta = C.c_interp_beta;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_cov_interp_beta = C.c_cov_interp_beta;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_interp_F = C.c_interp_F;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_interp_T = C.c_interp_T;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_interp_ess = C.c_interp_ess;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.c_interp_ess0 = C.c_interp_ess0;
    TOPO.v{side_hemi}.s{f1}.hb{hbi}.hb = hb;
end
end

function [SPM xCon SSxCon] = nirs_get_contrasts(SPM,job,automated_contrasts,s1,NonlinearEpilepsyOn)
%Construct the full design matrix over all sessions
SS_SPM = SPM;
SS_SPM.xCon = [];
if s1
    SS_SPM.xX = SPM.xXn{s1};
else
    SS_SPM.xX = SPM.xXn{1};
end

%spm contrasts
try
    tmp1 = []; tmp2 = []; tmp3 = 0; tmp4 = 0; tmp5 = 0; tmp6 = 0; tmp7 = [];
    for i1=1:length(SPM.xXn)
        %filtered design matrix
        tmp1 = blkdiag(tmp1,SPM.xXn{i1}.xKXs.X);
        %design matrix
        tmp2 = blkdiag(tmp2,SPM.xXn{i1}.X);
        %residual sum of squares of full model
        tmp3 = tmp3 + SPM.xXn{i1}.ResSS;
        tmp5 = tmp5 + SPM.xXn{i1}.trRV;
        tmp6 = tmp6 + SPM.xXn{i1}.trRVRV;
        tmp7 = blkdiag(tmp7,SPM.xXn{i1}.Bcov);
    end
    SPM.xX.xKXs = tmp1;
    SPM.xX.X = tmp2;
    SPM.xX.ResSS = tmp3;
    %Not clear from Satterthwaite approximation whether
    %to use tmp5/tmp6 or tmp4 for # of degrees of freedom
    %for multiple sessions
    SPM.xX.trRV = tmp5;
    SPM.xX.trRVRV = tmp6;
    SPM.xX.erdf = tmp5^2/tmp6;
    SPM.xX.erdf_check = tmp4;
    SPM.xX.corr_beta = tmp7;
    SPM.xX.var = SPM.xX.ResSS/SPM.xX.trRV;
    SPM.xCon = [];
    %Generate the SPM-type contrasts over all sessions
    SPM = nirs_spm_run_con(job,SPM);
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
%Generate automated contrasts if required
%number of regressors for one session
nr = size(SS_SPM.xX.X,2);
%clear contrast contrast_name xCon
%xCon = [];
contrastT = {}; contrastF = {};
contrastT_name = {}; contrastF_name = {};
%negative contrasts can be treated later as to avoid a
%duplication of long calculations
%NonlinearEpilepsyOn = 1; %Remember to revert back to 0, used for most studies
try
    if automated_contrasts
        if NonlinearEpilepsyOn
            try
                %number of confounds
                ncf = NIRSconfounds.NumChConfoundsActual;
            catch
                ncf = 0;
            end
            %Need to check if GLM was run with or without derivs
            switch SPM.xBF.name
                case {'hrf','Gamma functions'}
                    
                    %canonical HRF - no derivatives - no F contrasts
                    %1st Volterra kernel
                    contrastT{1} = [1 zeros(1,nr-1)];
                    contrastT_name{1} = 'T1';
                    %2nd Volterra kernel
                    if SPM.job.volt > 1
                        if nr > 5+ncf %Careful, this might not be the correct number
                            %if there are more confounding regressors
                            %assume 2 stimuli - only take the first one
                            contrastT{2} = [0 0 1 zeros(1,nr-3)];
                            contrastT_name{2} = 'T2';
                        else
                            %assume only 1 stimulus
                            contrastT{2} = [0 1 zeros(1,nr-2)];
                            contrastT_name{2} = 'T2';
                        end
                    end
                    if SPM.job.volt == 3
                        if nr > 5+ncf
                            %assume 2 stimuli - only take the first one
                            contrastT{3} = [0 0 0 0 0 0 1 zeros(1,nr-7)];
                            contrastT_name{3} = 'T3';
                        else
                            %assume only 1 stimulus
                            contrastT{3} = [0 0 1 zeros(1,nr-3)];
                            contrastT_name{3} = 'T3';
                        end
                    end
                case 'hrf (with time derivative)'
                    %Automated T contrasts
                    %1st Volterra kernel
                    contrastT{1} = [1 zeros(1,nr-1)];
                    contrastT_name{1} = 'T1';
                    %2nd Volterra kernel
                    if SPM.job.volt > 1
                        if nr > 10+ncf %Careful, this might not be the correct number
                            %if there are more confounding regressors
                            %assume 2 stimuli - only take the first one
                            contrastT{2} = [0 0 0 0 1 zeros(1,nr-7)];
                            contrastT_name{2} = 'T2';
                        else
                            %assume only 1 stimulus
                            contrastT{2} = [0 0 1 zeros(1,nr-4)];
                            contrastT_name{2} = 'T2';
                        end
                    end
                    %Automated F contrasts - careful, F contrasts need to be a cell
                    contrastF{1} = {[eye(2) zeros(2,nr-2)]};
                    contrastF_name{1} = 'F1';
                    if SPM.job.volt > 1
                        %temp_mat = [ zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
                        mat_Volt = [eye(3)]; %[eye(3) zeros(3,6); zeros(3,3) temp_mat];
                        if nr > 10+ncf
                            contrastF{2} = {[zeros(3,4) mat_Volt zeros(3,nr-7)]};
                            contrastF_name{2} = 'F2';
                        else
                            contrastF{2} = {[zeros(3,2) mat_Volt zeros(3,nr-5)]};
                            contrastF_name{2} = 'F2';
                        end
                    end
                    
                    %Automated F contrasts
                case 'hrf (with time and dispersion derivatives)'
                    %Automated T contrasts
                    %1st Volterra kernel
                    contrastT{1} = [1 zeros(1,nr-1)];
                    contrastT_name{1} = 'T1';
                    %2nd Volterra kernel
                    if SPM.job.volt > 1
                        if nr > 12+ncf %Careful, this might not be the correct number
                            %if there are more confounding regressors
                            %assume 2 stimuli - only take the first one
                            contrastT{2} = [0 0 0 0 0 0 1 zeros(1,nr-7)];
                            contrastT_name{2} = 'T2';
                        else
                            %assume only 1 stimulus
                            contrastT{2} = [0 0 0 1 zeros(1,nr-4)];
                            contrastT_name{2} = 'T2';
                        end
                    end
                    %Automated F contrasts - careful, F contrasts need to be a cell
                    contrastF{1} = {[eye(3) zeros(3,nr-3)]};
                    contrastF_name{1} = 'F1';
                    if SPM.job.volt > 1
                        %temp_mat = [ zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
                        mat_Volt = [eye(6)]; %[eye(3) zeros(3,6); zeros(3,3) temp_mat];
                        if nr > 12+ncf
                            contrastF{2} = {[zeros(6,6) mat_Volt zeros(6,nr-12)]};
                            contrastF_name{2} = 'F2';
                        else
                            contrastF{2} = {[zeros(6,3) mat_Volt zeros(6,nr-9)]};
                            contrastF_name{2} = 'F2';
                        end
                    end
            end
        else
            for i0=1:nr-1
                contrastT{i0} = [zeros(1,i0-1) 1 zeros(1,nr-i0)];
                try  
                    if s1
                        contrastT_name{i0} = ['_' validate_name(SS_SPM.Sess(s1).U(i0).name{1})];
                    else
                        contrastT_name{i0} = ['_' validate_name(SS_SPM.Sess(1).U(i0).name{1})];
                    end
                catch
                    contrastT_name{i0} = ['C' int2str(i0)];
                end
            end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
%Find out if some of the SPM contrasts could be used for
%single-session contrasts
try
    for i1=1:length(job.consess)
        if isfield(job.consess{i1},'tcon')
            l1 = length(job.consess{i1}.tcon.convec);
            if l1 <= nr
                %suitable for single session contrast
                %pad by zeros
                contrastT{end+1} = [job.consess{i1}.tcon.convec zeros(1,nr-l1)];
                contrastT_name{end+1} = job.consess{i1}.tcon.name;
            end
        else
            if isfield(job.consess{i1},'fcon')
                l1 = size(job.consess{i1}.fcon.convec{1},2);
                l2 = size(job.consess{i1}.fcon.convec{1},1);
                if l1 <= nr
                    %suitable for single session contrast
                    %pad by zeros
                    contrastF{end+1} = {[job.consess{i1}.fcon.convec{1} zeros(l2,nr-l1)]};
                    contrastF_name{end+1} = job.consess{i1}.fcon.name;
                end
            end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
%Generate single-session contrasts, using either automated
%contrasts or user-specified contrasts
try
    %store contrasts into SSxCon %Single-Session xCon
    if iscell(contrastT) == 1
        % for more than one contrast
        for kk = 1:length(contrastT)
            tmp_c = contrastT{kk};
            %store contrasts into xCon
            SSxCon(kk).name = contrastT_name{kk};
            SSxCon(kk).STAT = 'T';
            SSxCon(kk).c = tmp_c(:);
            %fill other fields
            SSxCon(kk).X0 = [];
            SSxCon(kk).iX0 = 'c'; %not used
            SSxCon(kk).X1o = [];
            SSxCon(kk).eidf = 1; %not used
            SSxCon(kk).Vcon = [];
            SSxCon(kk).Vspm = [];
        end
    end
    if iscell(contrastF) == 1
        nc1=length(SSxCon);
        % for more than one contrast
        for kk = 1:length(contrastF)
            tmp_c = contrastF{kk};
            %store contrasts into xCon
            SSxCon(kk+nc1).name = contrastF_name{kk};
            SSxCon(kk+nc1).STAT = 'F';
            SSxCon(kk+nc1).c = tmp_c{1}';
            %fill other fields
            SSxCon(kk+nc1).X0 = [];
            SSxCon(kk+nc1).iX0 = 'c'; %not used
            SSxCon(kk+nc1).X1o = [];
            SSxCon(kk+nc1).eidf = rank(tmp_c{1});
            SSxCon(kk+nc1).Vcon = [];
            SSxCon(kk+nc1).Vspm = [];
        end
    end
end
try
    jobt = [];
    for j1=1:length(contrastT)
        jobt.consess{j1}.tcon.name = contrastT_name{j1};
        jobt.consess{j1}.tcon.convec = contrastT{j1};
        jobt.consess{j1}.tcon.sessrep = 'both'; %'none'; %'both';
    end
    %contrasts get added to previously generated contrasts stored
    %in SPM structure
    SPM = nirs_spm_run_con(jobt,SPM);
end
try
    jobf = [];
    for j1=1:length(contrastF)
        jobf.consess{j1}.fcon.name = contrastF_name{j1};
        jobf.consess{j1}.fcon.convec = contrastF{j1};
        jobf.consess{j1}.fcon.sessrep = 'both'; %'none'; %'both';
    end
    SPM = nirs_spm_run_con(jobf,SPM);
end
%Contrasts over all sessions
xCon = SPM.xCon;

%replicate for SSxCon
try
    jobt = [];
    for j1=1:length(contrastT)
        jobt.consess{j1}.tcon.name = contrastT_name{j1};
        jobt.consess{j1}.tcon.convec = contrastT{j1};
        jobt.consess{j1}.tcon.sessrep = 'none';
    end
    SS_SPM = nirs_spm_run_con(jobt,SS_SPM);
end
try
    jobf = [];
    for j1=1:length(contrastF)
        jobf.consess{j1}.fcon.name = contrastF_name{j1};
        jobf.consess{j1}.fcon.convec = contrastF{j1};
        jobf.consess{j1}.fcon.sessrep = 'none';
    end
    SS_SPM = nirs_spm_run_con(jobf,SS_SPM);
end
%Single session contrasts
SSxCon = SS_SPM.xCon;
end

function name = validate_name(name)
%check if there are unallowed characters in file name part and remove them
    name = regexprep(name,'>','g');
    name = regexprep(name,'<','l');
    name = regexprep(name,' ','_');
    name = regexprep(name,'?','u');
    name = regexprep(name,'!','e');
    name = regexprep(name,'|','v');
    name = regexprep(name,'#','n');
    name = regexprep(name,'%','p');
    name = regexprep(name,'&','a');
    name = regexprep(name,'+','s');
    name = regexprep(name,'@','t');
    name = regexprep(name,'$','d');
    name = regexprep(name,'^','c');
    name = regexprep(name,'(','_');
    name = regexprep(name,')','_');
end