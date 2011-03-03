function out = nirs_run_liom_contrast(job)
%Views: 
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'

views_to_run = job.view;
contrast_data = job.contrast_data;

%Options
try
    p_value = job.contrast_p_value;
catch
    p_value = 0.05;
end
try
    cbar.c_min = job.override_colorbar.colorbar_override.colorbar_min;
    cbar.c_max = job.override_colorbar.colorbar_override.colorbar_max;
    cbar.colorbar_override = 1;
catch
    cbar.colorbar_override = 0;
end
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
try 
    GInv = job.GenerateInverted;
catch
    GInv = 1;
end
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
if ~isempty(contrast_data)
    automated_contrasts = 0;
else
    automated_contrasts = 1;
end
%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        NC = NIRS.Cf.H.C.N;
        %load topographic information (formerly known as preproc_info)
        try
            fname_ch = job.TopoData{1};
            load(fname_ch);
            NIRS.Dt.ana.rend = fname_ch;
        catch
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
        end
            
        %load SPM - first GLM - might want to generalize 
        dir1 = NIRS.SPM{1};
        load(fullfile(dir1,'SPM.mat'));
        ftopo = fullfile(dir1,'TOPO.mat');
        TOPO = [];
        try load(ftopo); end
        %Big loop over views 
        for v1=1:size(views_to_run,2)
            try
                brain_view = views_to_run(v1);
                switch brain_view
                    case 1 % 'ventral'
                        spec_hemi = 'ventral';
                        side_hemi = 1;
                    case 2 % 'dorsal'
                        spec_hemi = 'dorsal';
                        side_hemi = 2;
                    case 3 %'right_lateral'
                        spec_hemi = 'right';
                        side_hemi = 3;
                    case 4 %'left_lateral'
                        spec_hemi = 'left';
                        side_hemi = 4;
                    case 5 %'frontal'
                        spec_hemi = 'frontal';
                        side_hemi = 5;
                    case 6 %'occipital'
                        spec_hemi = 'occipital';
                        side_hemi = 6;
                end
                % channel information
                rchn = rendered_MNI{side_hemi}.rchn;
                cchn = rendered_MNI{side_hemi}.cchn;
                %rendering surface
                brain = rendered_MNI{side_hemi}.ren;
                brain = brain * 0.5;
                s1 = size(brain, 1);
                s2 = size(brain, 2);
                %find channels which are visible from this projection view
                index_ch = find(rchn ~= -1);
                %nch = length(index_ch);
                rchn = rchn(index_ch);
                cchn = cchn(index_ch);
                %number of regressors
                nr = size(SPM.xXn{1}.X,2);
                clear contrast contrast_name
                %negative contrasts can be treated later as to avoid a
                %duplication of long calculations
                if automated_contrasts
                    %1st Volterra kernel
                    contrast{1} = [1 zeros(1,nr-1)];
                    contrast_name{1} = '1';
                    %contrast{2} = [-1 zeros(1,nr-1)];
                    %contrast_name{2} = 'n1'; 
                    %2nd Volterra kernel
                    if SPM.job.volt > 1
                        if nr > 5 %Careful, this might not be the correct number
                            %if there are more confounding regressors
                            %assume 2 stimuli - only take the first one
                            contrast{2} = [0 0 1 zeros(1,nr-3)];
                            contrast_name{2} = '2';
                        else
                            %assume only 1 stimulus
                            contrast{2} = [0 1 zeros(1,nr-2)];
                            contrast_name{2} = '2';
                        end
                    end
                    if SPM.job.volt == 3
                        if nr > 5
                            %assume 2 stimuli - only take the first one
                            contrast{3} = [0 0 0 0 0 0 1 zeros(1,nr-7)];
                            contrast_name{3} = '3';
                        else
                            %assume only 1 stimulus
                            contrast{3} = [0 0 1 zeros(1,nr-3)];
                            contrast_name{3} = '3';
                        end 
                    end
                else
                    %user specified contrast specification
                    %simple T contrasts
                    for ic1=1:size(contrast_data,2)
                        tmp_c = contrast_data(ic1).contrast_c;
                        %pad with zeros
                        contrast{ic1} = [tmp_c zeros(1,nr-size(tmp_c,2))];
                        contrast_name{ic1} = contrast_data(ic1).contrast_name;
                    end
                end 
                %store contrasts into xCon
                if iscell(contrast) == 1
                    % for more than one contrast
                    for kk = 1:length(contrast)
                        tmp_c = contrast{kk};
                        %store contrasts into xCon
                        xCon(kk).c = tmp_c(:);
                        xCon(kk).n = contrast_name{kk};
                    end
                else
                    xCon(1).c = contrast(:);
                end

                %loop over sessions        
                for f1=1:length(SPM.xXn)
                    %split into HbO and HbR interpolations 
                    wl = NIRS.Cf.dev.wl;
                    %for NIRS acquisition with 2 wavelengths only
                    if wl(1) > 750
                        %first wavelength is "HbO-like"
                        ch_HbO = index_ch;
                        ch_HbR = NC/2 + index_ch;
                    else
                        ch_HbR = index_ch;
                        ch_HbO = NC/2 + index_ch; 
                    end
                    try
                        ch_HbT = NC+index_ch;
                    end
                    nC = length(xCon);
                    if GFIS                      
                        fh0Pt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                            num2str(p_value) '_S' int2str(f1) '_Pos'],'NumberTitle','off');
                        %subplot(fh0Pt,nC,3,1);
                        fh0Pu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                            num2str(p_value) '_S' int2str(f1) '_Pos'],'NumberTitle','off');
                        %subplot(fh0Pu,nC,3,1);
                        if GInv
                            fh0Nt = figure('Visible',cbar.visible,'Name',['A_tube_' ...
                                num2str(p_value) '_S' int2str(f1) '_Neg'],'NumberTitle','off');
                            %subplot(fh0Nt,nC,3,1);
                            fh0Nu = figure('Visible',cbar.visible,'Name',['A_unc_' ...
                                num2str(p_value) '_S' int2str(f1) '_Neg'],'NumberTitle','off');
                            %subplot(fh0Nu,nC,3,1);
                        else
                            fh0Nt = [];
                            fh0Nu = [];
                        end
                    else %not used
                        fh0Pt = [];
                        fh0Pu = [];
                        fh0Nt = [];
                        fh0Nu = [];
                    end
                    try
                        %for NIRS_SPM method
                        var = SPM.xXn{f1}.ResSS./SPM.xXn{f1}.trRV; 
                        %covariance of beta estimates
                        corr_beta = SPM.xXn{f1}.Bcov;
                    catch
                        %for WLS and BGLM methods
                        corr_beta = SPM.xXn{f1}.Bvar;
                        %will not work as we don't have var = ResSS/trRV
                    end                                        
                    %GLM estimates - which beta though???
                    beta_tmp = SPM.xXn{f1}.beta(:, ch_HbO);
                    beta_HbO = beta_tmp(:); %taken as one vector
                    beta_tmp = SPM.xXn{f1}.beta(:, ch_HbR);
                    beta_HbR = beta_tmp(:); %taken as one vector
                    mtx_var_HbO = diag(var(ch_HbO));
                    mtx_var_HbR = diag(var(ch_HbR));
                    try 
                        beta_tmp = SPM.xXn{f1}.beta(:, ch_HbT);
                        beta_HbT = beta_tmp(:); %taken as one vector
                        mtx_var_HbT = diag(var(ch_HbT)); 
                    end
                    erdf = SPM.xXn{f1}.erdf;
                    %HbO
                    %Note that var(ch_HbO) depends on HbO vs HbR
                    [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(...
                        corr_beta,length(ch_HbO),mtx_var_HbO,s1,s2,rchn,cchn);

                    [sum_kappa c_interp_beta c_cov_interp_beta] = ...
                        loop_contrasts(xCon,beta_HbO,corr_beta,cov_beta_r,...
                        mtx_var_HbO,B,Bx,By,rmask,cmask,s1,s2);
                    if gen_fig || gen_tiff
                        [fh0Pt,fh0Pu,fh0Nt,fh0Nu] = interpolated_maps(...
                            xCon,sum_kappa,c_interp_beta,c_cov_interp_beta,...
                            s1,s2,brain,spec_hemi,f1,dir1,erdf,p_value,'HbO',gen_fig,gen_tiff,cbar,GInv,...
                            GFIS,fh0Pt,fh0Pu,fh0Nt,fh0Nu);
                    end
                    TOPO.v{side_hemi}.s{f1}.hb{1}.sum_kappa = sum_kappa;
                    TOPO.v{side_hemi}.s{f1}.hb{1}.c_interp_beta = c_interp_beta;
                    TOPO.v{side_hemi}.s{f1}.hb{1}.c_cov_interp_beta = c_cov_interp_beta;

                    %HbR
                    [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(...
                        corr_beta,length(ch_HbR),mtx_var_HbR,s1,s2,rchn,cchn);              

                    [sum_kappa c_interp_beta c_cov_interp_beta] = ...
                        loop_contrasts(xCon,beta_HbR,corr_beta,cov_beta_r,...
                        mtx_var_HbR,B,Bx,By,rmask,cmask,s1,s2);
                    if gen_fig || gen_tiff
                        [fh0Pt,fh0Pu,fh0Nt,fh0Nu] = interpolated_maps(...
                            xCon,sum_kappa,c_interp_beta,c_cov_interp_beta,...
                            s1,s2,brain,spec_hemi,f1,dir1,erdf,p_value,'HbR',gen_fig,gen_tiff,cbar,GInv,...
                            GFIS,fh0Pt,fh0Pu,fh0Nt,fh0Nu);
                    end
                    TOPO.v{side_hemi}.s{f1}.hb{2}.sum_kappa = sum_kappa;
                    TOPO.v{side_hemi}.s{f1}.hb{2}.c_interp_beta = c_interp_beta;
                    TOPO.v{side_hemi}.s{f1}.hb{2}.c_cov_interp_beta = c_cov_interp_beta;              

                    try
                        %HbT
                        [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(...
                            corr_beta,length(ch_HbT),mtx_var_HbT,s1,s2,rchn,cchn);              

                        [sum_kappa c_interp_beta c_cov_interp_beta] = ...
                            loop_contrasts(xCon,beta_HbT,corr_beta,cov_beta_r,...
                            mtx_var_HbT,B,Bx,By,rmask,cmask,s1,s2);
                        if gen_fig || gen_tiff
                            [fh0Pt,fh0Pu,fh0Nt,fh0Nu] = interpolated_maps(...
                                xCon,sum_kappa,c_interp_beta,c_cov_interp_beta,...
                                s1,s2,brain,spec_hemi,f1,dir1,erdf,p_value,'HbT',gen_fig,gen_tiff,cbar,GInv,...
                                GFIS,fh0Pt,fh0Pu,fh0Nt,fh0Nu);
                        end
                        TOPO.v{side_hemi}.s{f1}.hb{3}.sum_kappa = sum_kappa;
                        TOPO.v{side_hemi}.s{f1}.hb{3}.c_interp_beta = c_interp_beta;
                        TOPO.v{side_hemi}.s{f1}.hb{3}.c_cov_interp_beta = c_cov_interp_beta;   
                    end
                    if GFIS 
                        save_session_figures(fh0Pt,gen_fig,gen_tiff,p_value,'Pos',dir1,spec_hemi,cbar,'tube',f1);
                        save_session_figures(fh0Pu,gen_fig,gen_tiff,p_value,'Pos',dir1,spec_hemi,cbar,'unc',f1);
                        if GInv
                            save_session_figures(fh0Nt,gen_fig,gen_tiff,p_value,'Neg',dir1,spec_hemi,cbar,'tube',f1);
                            save_session_figures(fh0Nu,gen_fig,gen_tiff,p_value,'Neg',dir1,spec_hemi,cbar,'unc',f1);
                        end
                    end
                end %end for f1
                TOPO.v{side_hemi}.s1 = s1; %sizes of topographic projection
                TOPO.v{side_hemi}.s2 = s2;
                TOPO.v{side_hemi}.view = spec_hemi; %%% view of the brain
                %TOPO.v{side_hemi}.side_hemi = side_hemi;
            catch
                disp(['Could not create contrasts for view ' spec_hemi ' for subject ' int2str(Idx)]);
            end
        end %end for v1
        TOPO.xCon = xCon; %would not work if new contrasts are later added        
        save(ftopo,'TOPO');
    catch
        disp(['Could not create contrasts for subject' int2str(Idx)]);
    end
    NIRS.TOPO = ftopo;
    save(job.NIRSmat{Idx,1},'NIRS');
end
out.NIRSmat = job.NIRSmat;
end

function [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(corr_beta,nch,...
    mtx_var,s1,s2,rchn,cchn)
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
    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
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
end

function [sum_kappa c_interp_beta c_cov_interp_beta] =...
    loop_contrasts(xCon,beta,corr_beta,cov_beta_r,mtx_var,B,Bx,By,rmask,cmask,s1,s2)
%PP Big loop over individual T-contrasts -- how to generalize to
%F-contrasts?
rmask_vector = rmask{1};
cmask_vector = cmask{1};
nC = size(xCon,2);
nm = length(rmask_vector);
%preallocate
sum_kappa = zeros(nC,1);
kappa = zeros(nC,s1,s2);
c_interp_beta = zeros(nC,s1,s2);
c_cov_interp_beta = zeros(nC,s1,s2);
sz_xCon  = size(xCon(1).c,1);
% B2 = zeros(nch, 1);
% B2x = zeros(nch, 1);
% B2y = zeros(nch, 1);
%identity matrix of size number of regressors
tmp = eye(sz_xCon);
for kk = 1:nm
    %this is different for HbO and HbR     
    B2(:,1) = B(rmask_vector(kk), cmask_vector(kk), :);
    B2x(:,1) = Bx(rmask_vector(kk), cmask_vector(kk), :);
    B2y(:,1) = By(rmask_vector(kk), cmask_vector(kk), :); 
    B3 = kron(B2, tmp);
    B3x = kron(B2x, tmp);
    B3y = kron(B2y, tmp);
    B3t = kron(B2', tmp);
    for c1 = 1:nC
        %this is the same for HbO and HbR
        c_corr_beta = xCon(c1).c' * corr_beta * xCon(c1).c; 
        %this is different for HbO and HbR
        P = cov_beta_r * (B3* xCon(c1).c);
        Px = cov_beta_r * (B3x * xCon(c1).c);
        Py = cov_beta_r * (B3y * xCon(c1).c);
        tmp_1 = P'*P; tmp_2 = tmp_1^(-1/2); tmp_3 = tmp_2^3; 
        u_derx = Px*tmp_2 - (P*(P'*Px))*tmp_3;
        u_dery = Py*tmp_2 - (P*(P'*Py))*tmp_3;
        %For each contrast and each channel, we get kappa, c_interp_beta
        %and c_cov_interp_beta
        
        %Positive contrasts
        kappa(c1,rmask_vector(kk), cmask_vector(kk)) =  ...
           sqrt(abs(det([u_derx'*u_derx u_derx'*u_dery; u_dery'*u_derx u_dery'*u_dery])));
        c_interp_beta(c1,rmask_vector(kk), cmask_vector(kk)) = ...
                (xCon(c1).c' * B3t) * beta;
        c_cov_interp_beta(c1,rmask_vector(kk), cmask_vector(kk)) = ...
                (B2'*mtx_var*B2) * c_corr_beta;
           
        %%%%%%INSTEAD, JUST RECALL THAT ONLY c_interp_beta flips sign%%%%%    
        %Negative contrasts  
        %c_corr_beta is unchanged, P, Px, Py flip sign
        %tmp_1,2,3,4 are unchanged, u_derx, u_dery flip sign
        %kappa is unchanged, c_interp_beta flips sign
        %c_cov_interp_beta is unchanged
% %         kappa(2*c1,rmask_vector(kk), cmask_vector(kk)) = ...
% %             kappa(2*c1-1,rmask_vector(kk), cmask_vector(kk));
% %         c_interp_beta(2*c1,rmask_vector(kk), cmask_vector(kk)) = ...
% %             -c_interp_beta(2*c1-1,rmask_vector(kk), cmask_vector(kk));     
% %         c_cov_interp_beta(2*c1,rmask_vector(kk), cmask_vector(kk)) = ...
% %             c_cov_interp_beta(2*c1-1,rmask_vector(kk), cmask_vector(kk));    
    end
end
for c1 = 1:nC
    tm = kappa(c1,:,:);
    sum_kappa(c1) = sum(tm(:)); 
end
end    


function [fh0Pt,fh0Pu,fh0Nt,fh0Nu] = interpolated_maps(xCon,sum_kappa,c_interp_beta,c_cov_interp_beta,...
    s1,s2,brain,spec_hemi,f1,pathn,erdf,p_value,hb,gen_fig,gen_tiff,cbar,GInv,...
    GFIS,fh0Pt,fh0Pu,fh0Nt,fh0Nu)
%loop over contrasts
nC = size(xCon,2);
for c1=1:nC
    index_mask = find(squeeze(c_cov_interp_beta(c1,:,:)) ~= 0);
    T_map = zeros(s1, s2);
    T_map(index_mask) = squeeze(c_interp_beta(c1,index_mask))./ ...
        sqrt(squeeze(c_cov_interp_beta(c1,index_mask)));
    %names
    contrast_info = [num2str(p_value) '_' spec_hemi '_' hb '_S' int2str(f1) '_Pos' xCon(c1).n];
    contrast_info_for_fig = [num2str(p_value) ' ' spec_hemi ' ' hb '_S' int2str(f1) ' Pos' xCon(c1).n];
    
    load Split
%     %no threshold
%     nirs_draw_figure(1,brain,T_map,contrast_info,...
%         contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar);
    %uncorrected
    if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
        [fh1 ax1 hc1] = nirs_draw_figure(2,brain,T_map,contrast_info,...
            contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar); 
        if GFIS, [fh0Pu, fh0Nu] = nirs_copy_figure(fh0Pu,fh0Nu,c1,nC,hb,GInv,fh1,ax1,hc1,1,split); end
        %tube
        [fh1 ax1 hc1] = nirs_draw_figure(3,brain,T_map,contrast_info,...
        contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar); 
        if GFIS, [fh0Pt, fh0Nt] = nirs_copy_figure(fh0Pt,fh0Nt,c1,nC,hb,GInv,fh1,ax1,hc1,1,split); end
    end
    %repeat for negative contrasts
    T_map = - T_map;
    contrast_info = [num2str(p_value) '_' spec_hemi '_' hb '_S' int2str(f1) '_Neg' xCon(c1).n];
    contrast_info_for_fig = [num2str(p_value) ' ' spec_hemi ' ' hb '_S' int2str(f1) ' Neg' xCon(c1).n];

%     %no threshold
%     nirs_draw_figure(1,brain,T_map,contrast_info,...
%         contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar);
    if GInv || strcmp(hb,'HbR')
        %uncorrected
        [fh1 ax1 hc1] = nirs_draw_figure(2,brain,T_map,contrast_info,...
            contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar);  
        if GFIS, [fh0Pu, fh0Nu] = nirs_copy_figure(fh0Pu,fh0Nu,c1,nC,hb,GInv,fh1,ax1,hc1,0,split); end
        %tube
        [fh1 ax1 hc1] = nirs_draw_figure(3,brain,T_map,contrast_info,...
            contrast_info_for_fig,split,pathn,erdf,sum_kappa(c1),p_value,gen_fig,gen_tiff,cbar); 
        if GFIS, [fh0Pt, fh0Nt] = nirs_copy_figure(fh0Pt,fh0Nt,c1,nC,hb,GInv,fh1,ax1,hc1,0,split); end
    end
end
end

% figure(1);
% ax1=axes('Parent',1);
% plot(ax1,rand(100,1));
% title(ax1,'some title');
% xlabel(ax1,'time');
% ylabel(ax1,'amplitude');
% 
% figure(2);
% ax2=subplot(211);
% copyobj(allchild(ax1),ax2);

function save_session_figures(fh0,gen_fig,gen_tiff,p_value,Inv,dir1,spec_hemi,cbar,str_cor,f1) 
filen1 = fullfile(dir1,[str_cor '_' num2str(p_value) '_' spec_hemi '_S' int2str(f1) '_' Inv '.fig']);
if gen_fig
    saveas(fh0,filen1,'fig');
end
if gen_tiff
    filen2 = fullfile(dir1,['A_' str_cor '_' num2str(p_value) '_' spec_hemi '_S' int2str(f1) '_' Inv '.tiff']);
    print(fh0, '-dtiffn', filen2);
end
if strcmp(cbar.visible, 'off')
    try close(fh0); end
end
end