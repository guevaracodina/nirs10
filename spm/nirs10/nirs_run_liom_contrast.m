function out = nirs_run_liom_contrast(job)
%Views: 
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'

%views_to_run = [4 3]; %[4 3];
views_to_run = job.view;
contrast_data = job.contrast_data;
automated_contrasts = 1;
%Options
p_value = 0.05;
%flag_correction = 0;
flag_figure = 1;

%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        %load topographic information (formerly known as preproc_info)
        fname_ch = NIRS.Dt.ana.rend;
        load(fname_ch);
        %load SPM
        load(NIRS.SPM);
        
        %Big loop over views 
        for v1=1:size(views_to_run,2)
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
            nch = length(index_ch);
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
                contrast_name{1} = 'p1';
                %contrast{2} = [-1 zeros(1,nr-1)];
                %contrast_name{2} = 'n1'; 
                %2nd Volterra kernel
                if SPM.job.volt > 1
                    if nr > 4
                        %assume 2 stimuli - only take the first one
                        contrast{2} = [0 0 1 zeros(1,nr-3)];
                        contrast_name{2} = 'p2';
                        %contrast{4} = [0 0 -1 zeros(1,nr-3)];
                        %contrast_name{4} = 'n2'; 
                    else
                        %assume only 1 stimulus
                        contrast{2} = [0 1 zeros(1,nr-2)];
                        contrast_name{2} = 'p2';
                        %contrast{4} = [0 -1 zeros(1,nr-2)];
                        %contrast_name{4} = 'n2'; 
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
                               
                %HbO
                %Note that var(ch_HbO) depends on HbO vs HbR
                sz_xCon  = size(xCon(1).c,1);
                [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(...
                    corr_beta,length(ch_HbO),var(ch_HbO),...
                    sz_xCon,s1,s2,rchn,cchn);
                
                [sum_kappa c_interp_beta c_cov_interp_beta] = ...
                    loop_contrasts(xCon,corr_beta,cov_beta_r,...
                    B,Bx,By,rmask,cmask,s1,s2);
                
                TOPO.v{v1}.s{f1}.HbO.sum_kappa = sum_kappa;
                TOPO.v{v1}.s{f1}.HbO.c_interp_beta = c_interp_beta;
                TOPO.v{v1}.s{f1}.HbO.c_cov_interp_beta = c_cov_interp_beta;
                
                %HbR
                [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(...
                    corr_beta,length(ch_HbR),var(ch_HbR),...
                    sz_xCon,s1,s2,rchn,cchn);              

                [sum_kappa c_interp_beta c_cov_interp_beta] = ...
                    loop_contrasts(xCon,corr_beta,cov_beta_r,...
                    B,Bx,By,rmask,cmask,s1,s2);
                
                TOPO.v{v1}.s{f1}.HbR.sum_kappa = sum_kappa;
                TOPO.v{v1}.s{f1}.HbR.c_interp_beta = c_interp_beta;
                TOPO.v{v1}.s{f1}.HbR.c_cov_interp_beta = c_cov_interp_beta;              
  
           end %end for f1
           TOPO.v{v1}.s1 = s1; %sizes of topographic projection
           TOPO.v{v1}.s2 = s2;
           TOPO.v{v1}.view = spec_hemi; %%% view of the brain
           TOPO.v{v1}.side_hemi = side_hemi;
        end %end for v1
        TOPO.xCon = xCon;
        [dir1, ~, ~] = fullfile(NIRS.SPM);
        ftopo = fullfile(dir1,'TOPO.mat');
        save(ftopo,'TOPO');
    catch
        disp(['Could not create contrasts for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;
end

function [cov_beta_r B Bx By rmask cmask] = interpolation_kernel(corr_beta,nch,...
    var,sz_xCon,s1,s2,rchn,cchn)
%identity over remaining channels
mtx_eye = eye(nch);
mtx_var = diag(var);
cov_beta = kron(mtx_var, corr_beta);
[U, S, V] = svd(cov_beta);
cov_beta_r = U*(S.^(0.5))*V';
%identity matrix of size 1??
tmp = eye(sz_xCon);
%mesh of topographically projected brain size; note s1 and
%s1 inverted
[x, y] = meshgrid(1:s2, 1:s1);

B2 = zeros(nch, 1);
B2x = zeros(nch, 1);
B2y = zeros(nch, 1);
%?
B = zeros(s1, s2, nch);
Bx = zeros(s1, s2, nch);
By = zeros(s1, s2, nch);
%masks
rmask{1} = [];
cmask{1} = [];

disp('Extracting interpolation kernels...');
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
    loop_contrasts(xCon,corr_beta,cov_beta_r,B,Bx,By,rmask,cmask,s1,s2)
%PP Big loop over individual T-contrasts -- how to generalize to
%F-contrasts?
rmask_vector = rmask{1};
cmask_vector = cmask{1};
nCon = size(xCon,2);
nm = length(rmask_vector);
%preallocate
sum_kappa = zeros(nCon,1);
kappa = zeros(nCon,s1,s2);
c_interp_beta = zeros(nCon,s1,s2);
c_cov_interp_beta = zeros(nCon,s1,s2);
for kk = 1:nm
    %this is different for HbO and HbR     
    B2(:,1) = B(rmask_vector(kk), cmask_vector(kk), :);
    B2x(:,1) = Bx(rmask_vector(kk), cmask_vector(kk), :);
    B2y(:,1) = By(rmask_vector(kk), cmask_vector(kk), :); 
    for Ic = 1:nCon
        %this is the same for HbO and HbR
        c_corr_beta = xCon(Ic).c' * corr_beta * xCon(Ic).c; 
        %this is different for HbO and HbR
        P = cov_beta_r * kron(B2, tmp)* xCon(Ic).c;
        Px = cov_beta_r * kron(B2x, tmp) * xCon(Ic).c;
        Py = cov_beta_r * kron(B2y, tmp) * xCon(Ic).c;
        tmp_1 = P'*P; tmp_2 = tmp_1^(-1/2); tmp_3 = tmp_2^3; tmp_4 = P*P';
        u_derx = Px.*tmp_2 - (tmp_4*Px).*tmp_3;
        u_dery = Py.*tmp_2 - (tmp_4*Py).*tmp_3;
        %For each contrast and each channel, we get kappa, c_interp_beta
        %and c_cov_interp_beta
        
        %Positive contrasts
        kappa(Ic,rmask_vector(kk), cmask_vector(kk)) =  ...
           sqrt(abs(det([u_derx'*u_derx u_derx'*u_dery; u_dery'*u_derx u_dery'*u_dery])));
        c_interp_beta(Ic,rmask_vector(kk), cmask_vector(kk)) = ...
                xCon(Ic).c' * kron(B2', eye(size(xCon(Ic).c,1))) * beta;
        c_cov_interp_beta(Ic,rmask_vector(kk), cmask_vector(kk)) = ...
                (B2'*mtx_var*B2) * c_corr_beta;
           
        %%%%%%INSTEAD, JUST RECALL THAT ONLY c_interp_beta flips sign%%%%%    
        %Negative contrasts  
        %c_corr_beta is unchanged, P, Px, Py flip sign
        %tmp_1,2,3,4 are unchanged, u_derx, u_dery flip sign
        %kappa is unchanged, c_interp_beta flips sign
        %c_cov_interp_beta is unchanged
% %         kappa(2*Ic,rmask_vector(kk), cmask_vector(kk)) = ...
% %             kappa(2*Ic-1,rmask_vector(kk), cmask_vector(kk));
% %         c_interp_beta(2*Ic,rmask_vector(kk), cmask_vector(kk)) = ...
% %             -c_interp_beta(2*Ic-1,rmask_vector(kk), cmask_vector(kk));     
% %         c_cov_interp_beta(2*Ic,rmask_vector(kk), cmask_vector(kk)) = ...
% %             c_cov_interp_beta(2*Ic-1,rmask_vector(kk), cmask_vector(kk));    
    end
end
for Ic = 1:nCon
    tm = kappa(Ic,:,:);
    sum_kappa(Ic) = sum(tm(:)); 
end
end    


function interpolated_maps(xCon,sum_kappa,c_interp_beta,c_cov_interp_beta,s1,s2)
%loop over contrasts
index_mask = find(c_cov_interp_beta ~= 0);
T_map = zeros(s1, s2);
T_map(index_mask) = c_interp_beta(index_mask)./sqrt(c_cov_interp_beta(index_mask));
min_T = min(T_map(index_mask));
max_T = max(T_map(index_mask));

smin_T = max_T - ((max_T - min_T)./63) * 127;
sbar = linspace(smin_T, max_T, 128);
T_brain = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
T_brain(index_mask) = T_map(index_mask);

contrast_info = ['cinterp_SPM_nirs_' name '_Contrast_' contrast_name{Ic} '_View_' spec_hemi];
contrast_info_for_fig = ['cinterp SPM nirs ' name ' Contrast ' contrast_name{Ic} ' View ' spec_hemi];
%if flag_correction
    str_cor1 = 'tube';
%else
    str_cor2 = 'unc';
%end
if flag_figure == 1
    fh1 = figure('Name',[contrast_info],'NumberTitle','off');
    title(contrast_info_for_fig);
    imagesc(T_brain); %,'cdatamapping','scaled');
    load Split
    colormap(split)
    axis off
    axis image
    hc = colorbar;
    set(hc, 'YLim', [sbar(65) sbar(128)]);
    y_tick = linspace(sbar(65), sbar(128), 5)';
    set(hc, 'YTick', y_tick);
    set(hc, 'FontSize', 8);
end


filen_interp_SPM = fullfile(pathn,[contrast_info '.mat']);
save(filen_interp_SPM, 'cinterp_SPM_nirs');

erdf = SPM_nirs.xX.erdf;
%if flag_correction == 1 % tube formula correction
    z_value = 1:0.0001:7;
    p_value_tube = ((sum_kappa * gamma(3/2))/(2*(pi^(3/2))))*(1-gammainc((z_value(:).^2)/2, 3/2));
    index_z = [];
    ini_ran = 10^(-10);
    n = 0;
    while isempty(index_z) == 1
        ran = ini_ran * (10^n);
        n = n+1;
        index_z = find(p_value_tube > p_value - ran & p_value_tube < p_value + ran);
    end
    index_z = index_z(end);
    th_z = z_value(index_z);

fh2 = nirs_SPM_NIRS_draw_figure(th_z,brain,contrast_info,contrast_info_for_fig,T_map,flag_figure,str_cor1,split);
%if flag_correction == 0 % no correction
    th_z = spm_invTcdf(1-p_value, erdf);
%end
fh3 = nirs_SPM_NIRS_draw_figure(th_z,brain,contrast_info,contrast_info_for_fig,T_map,flag_figure,str_cor2,split);

filen1 = fullfile(pathn,[contrast_info '_noT.fig']);
filen2 = fullfile(pathn,[contrast_info '_' str_cor1 '.fig']);
filen3 = fullfile(pathn,[contrast_info '_' str_cor2 '.fig']);
saveas(fh1,filen1,'fig');
saveas(fh2,filen2,'fig');
saveas(fh3,filen3,'fig');

print(fh1,'-dpsc2','-append', ResultsFile);
print(fh2,'-dpsc2','-append', ResultsFile);
print(fh3,'-dpsc2','-append', ResultsFile);
try close(fh1); end
try close(fh2); end
try close(fh3); end
end