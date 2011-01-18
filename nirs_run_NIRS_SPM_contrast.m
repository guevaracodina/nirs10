function out = nirs_run_NIRS_SPM_contrast(job)
%Load NIRS.mat information
%clear NIRS
%load(job.NIRSmat{1,1});
fname_ch = job.NIRS_SPM_Coregistration_Channels{1,1};
views_to_run = job.view;
contrast_data = job.contrast_data;
try 
    load(fname_ch);
    preproc_info.ch_config.nholder = '1set';
catch
    disp('Could not load info on channel coregistration - aborting');
    out = [];
    return;
end
Dmx_files = job.Dmx_files;

%Options
p_value = 0.05;
%flag_correction = 0;
flag_figure = 1;

%Views: 
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'

%views_to_run = [4 3]; %[4 3];
%brain_view = 'left_lateral';

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
    
    brain = preproc_info.rend_ch_pos{side_hemi}.ren;
    brain = brain * 0.5;
    s1 = size(brain, 1);
    s2 = size(brain, 2);

    % channel information
    rchn = preproc_info.rend_ch_pos{side_hemi}.rchn;
    cchn = preproc_info.rend_ch_pos{side_hemi}.cchn;

    index_ch = find(rchn ~= -1);
    if strcmp(preproc_info.ch_config.nholder, '2sets') == 1
        ch1 = preproc_info.ch_config.ch_holder1;
        ch2 = preproc_info.ch_config.ch_holder2;
        index_tmp1 = find(rchn(ch1) ~= -1);
        index_tmp2 = find(rchn(ch2) ~= -1);
        index_side1_ch = ch1(index_tmp1);
        index_side2_ch = ch2(index_tmp2);
    end
    nch = length(index_ch);
    rchn = rchn(index_ch);
    cchn = cchn(index_ch);
        
     %Loop over contrasts: done already in activation_map_batch for
        %T-contrasts

        %Loop over views

        %Loop over methods of false positive corrections
%          contrast = [1 0 0];
%          contrast_name{1} = 'Positive';

        %contrast{1} = [1 0 0 ];
        %contrast_name{1} = 'Positive';
        %contrast{1} = [-1 0 0];
        %contrast_name{1} = 'Negative';

    %     contrast{1} = [1 0 0 0 0 0];
    %     contrast_name{1} = 'Can';
    %     contrast{2} = [1 0 -1 0 0 0];
    %     contrast_name{2} = 'CanXCan';

        %contrast{2} = [-1 0 0 0];

        
    %Loop over files %might want to do that later - loop over views first
    for f1=1:size(Dmx_files,1)
        try
        fname_SPM = Dmx_files{f1};
        clear SPM_nirs
        load(fname_SPM);
        try
            [pathn, name, ext, versn] = fileparts(fname_SPM);
            pathn = [pathn filesep];
        catch
            index = find(fname_SPM == filesep);
            pathn = fname_SPM(1:pathn(end));
        end

        ResultsFile = fullfile(pathn,'ResultsAll.ps');

        %complete the contrast specification
        nbeta = size(SPM_nirs.xX.X,2);
        %simple T contrasts
        for ic1=1:size(contrast_data,2)
            tmp_c = contrast_data(ic1).contrast_c;
            %pad with zeros
            contrast{ic1} = [tmp_c zeros(1,nbeta-size(tmp_c,2))];
            contrast_name{ic1} = contrast_data(ic1).contrast_name;
        end
        
        %[T_map, T_brain, T_brain_over] = activation_map_batch(fname_SPM, fname_ch, contrast, brain_view, p_value, flag_correction, flag_figure);


        % for more than one contrast
        if iscell(contrast) == 1
            for kk = 1:length(contrast)
                tmp_c = contrast{kk};
                xCon(kk).c = tmp_c(:);
            end
        else
            xCon(1).c = contrast(:);
        end

        var = SPM_nirs.nirs.ResSS./SPM_nirs.xX.trRV; %PP note trRV depends only on design matrix and is independent of HbO, HbR, HbT
        corr_beta = SPM_nirs.xX.Bcov;
        mtx_eye = eye(nch);
        beta_tmp = SPM_nirs.nirs.beta(:, index_ch);
        beta = beta_tmp(:);
        mtx_var = zeros(nch,nch);
        for kk = 1:nch
            mtx_var(kk,kk) = var(index_ch(kk));
        end
        cov_beta = kron(mtx_var, corr_beta);
        [U, S, V] = svd(cov_beta);
        cov_beta_r = U*(S.^(0.5))*V';
        tmp = eye(size(xCon(1).c,1));

        %[Ic, xCon] = nirs_spm_conman(SPM_nirs, 'T', Inf, 'Select contrasts...', 'for conjunction', 1);
        %SPM_nirs.xCon = xCon;

        [x, y] = meshgrid(1:s2, 1:s1);
        c_interp_beta = zeros(s1, s2);
        c_cov_interp_beta = zeros(s1, s2);
        kappa = zeros(s1, s2);
        B2 = zeros(nch, 1);
        B2x = zeros(nch, 1);
        B2y = zeros(nch, 1);
        B = zeros(s1, s2, nch);
        Bx = zeros(s1, s2, nch);
        By = zeros(s1, s2, nch);
        rmask{1} = [];
        rmask{2} = [];

        switch preproc_info.ch_config.nholder
            case '1set'
                disp('Extracting interpolation kernels...');
                for kk = 1:nch
                    grid_eye = griddata(cchn, rchn, (mtx_eye(:,kk))', x, y, 'cubic');
                    if  kk == 1
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
            case '2sets'
                if isempty(index_side1_ch) ~= 1
                    nch_side1 = length(index_side1_ch);
                    disp('Extracting interpolation kernels...(1set)');
                    eye_tmp = eye(nch_side1);
                    rchn = preproc_info.rend_ch_pos{side_hemi}.rchn;
                    cchn = preproc_info.rend_ch_pos{side_hemi}.cchn;

                    for kk = 1:nch_side1
                        grid_eye = griddata(cchn(index_side1_ch), rchn(index_side1_ch), (eye_tmp(1:nch_side1, kk))', x, y, 'cubic');
                        if kk == 1
                            mask = isnan(grid_eye);
                            mask = 1 - mask;
                            [rmask{1}, cmask{1}] = find(mask == 1);
                            index_mask0 = find(mask == 0);
                        end
                        grid_eye(index_mask0) = 0;
                        [temp_Bx, temp_By] = gradient(grid_eye);
                        temp_Bx(index_mask0) = 0;
                        temp_By(index_mask0) = 0;
                        B(:, :, kk) = grid_eye;
                        Bx(:, :, kk) = temp_Bx;
                        By(:, :, kk) = temp_By;
                    end
                else
                    nch_side1 = 0;
                end
                if isempty(index_side2_ch) ~= 1
                    nch_side2 = length(index_side2_ch);
                    disp('Extracting interpolation kernels...(2set)');
                    eye_tmp = eye(nch_side2);
                    rchn = preproc_info.rend_ch_pos{side_hemi}.rchn;
                    cchn = preproc_info.rend_ch_pos{side_hemi}.cchn;

                    for kk = 1:nch_side2
                        grid_eye = griddata(cchn(index_side2_ch), rchn(index_side2_ch), (eye_tmp(1:nch_side2, kk))', x, y, 'cubic');
                        if kk == 1
                            mask = isnan(grid_eye);
                            mask = 1 - mask;
                            [rmask{2}, cmask{2}] = find(mask == 1);
                            index_mask0 = find(mask == 0);
                        end
                        grid_eye(index_mask0) = 0;
                        [temp_Bx, temp_By] = gradient(grid_eye);
                        temp_Bx(index_mask0) = 0;
                        temp_By(index_mask0) = 0;
                        B(:,:, kk+nch_side1) = grid_eye;
                        Bx(:,:,kk+nch_side1) = temp_Bx;
                        By(:,:,kk+nch_side1) = temp_By;
                    end
                end
        end  %%% end of case '2sets'
        clear index_mask0;
        clear temp_Bx
        clear temp_By
        clear grid_eye
        clear x
        clear y
        %PP Big loop over individual T-contrasts -- how to generalize to
        %F-contrasts?
        for Ic = 1:size(xCon,2)
            if Ic == 1
                disp('for the 1st contrast vector,');
            elseif Ic == 2
                disp('for the 2nd contrast vector,');
            elseif Ic == 3
                disp('for the 3rd contrast vector,');
            elseif Ic > 3
                disp(['for the ' num2str(Ic) 'th contrast vector,']);
            end
            c_corr_beta = xCon(Ic).c' * corr_beta * xCon(Ic).c;
            try
                disp('generating interpolated response of HbX and its covariance from the 1st channel set...');
                rmask_vector = rmask{1};
                cmask_vector = cmask{1};
                for kk = 1:length(rmask_vector)
                    B2(:,1) = B(rmask_vector(kk), cmask_vector(kk), :);
                    B2x(:,1) = Bx(rmask_vector(kk), cmask_vector(kk), :);
                    B2y(:,1) = By(rmask_vector(kk), cmask_vector(kk), :);                    
                    P = cov_beta_r * kron(B2, tmp)* xCon(Ic).c;
                    Px = cov_beta_r * kron(B2x, tmp) * xCon(Ic).c;
                    Py = cov_beta_r * kron(B2y, tmp) * xCon(Ic).c;
                    %PP
                    tmp_1 = P'*P; tmp_2 = tmp_1^(-1/2); tmp_3 = tmp_2^3; tmp_4 = P*P';
%                     u_derx = Px./((P'*P)^(1/2)) - (P*P'*Px)./((P'*P)^(3/2));
%                     u_dery = Py./((P'*P)^(1/2)) - (P*P'*Py)./((P'*P)^(3/2));
                    u_derx = Px.*tmp_2 - (tmp_4*Px).*tmp_3;
                    u_dery = Py.*tmp_2 - (tmp_4*Py).*tmp_3;
                    kappa(rmask_vector(kk), cmask_vector(kk)) =  abs(det([u_derx'*u_derx u_derx'*u_dery; u_dery'*u_derx u_dery'*u_dery]));
                    c_interp_beta(rmask_vector(kk), cmask_vector(kk)) = xCon(Ic).c' * kron(B2', eye(size(xCon(Ic).c,1))) * beta;
                    c_cov_interp_beta(rmask_vector(kk), cmask_vector(kk)) = (B2'*mtx_var*B2) * c_corr_beta;
                end
            end
            if strcmp(preproc_info.ch_config.nholder, '2sets') == 1 && isempty(rmask{2}) ~= 1
                disp('generating interpolated response of HbX and its covariance from the 2nd channel set...');
                rmask_vector = rmask{2};
                cmask_vector = cmask{2};
                for kk = 1:length(rmask_vector)
                    B2(:,1) = B(rmask_vector(kk), cmask_vector(kk), :);
                    B2x(:,1) = Bx(rmask_vector(kk), cmask_vector(kk), :);
                    B2y(:,1) = By(rmask_vector(kk), cmask_vector(kk), :);
                    P = cov_beta_r * kron(B2, tmp)* xCon(Ic).c;
                    Px = cov_beta_r * kron(B2x, tmp) * xCon(Ic).c;
                    Py = cov_beta_r * kron(B2y, tmp) * xCon(Ic).c;
                    u_derx = Px./((P'*P)^(1/2)) - (P*P'*Px)./((P'*P)^(3/2));
                    u_dery = Py./((P'*P)^(1/2)) - (P*P'*Py)./((P'*P)^(3/2));
                    kappa(rmask_vector(kk), cmask_vector(kk)) =  abs(det([u_derx'*u_derx u_derx'*u_dery; u_dery'*u_derx u_dery'*u_dery]));
                    c_interp_beta(rmask_vector(kk), cmask_vector(kk)) = xCon(Ic).c' * kron(B2', eye(size(xCon(Ic).c,1))) * beta;
                    c_cov_interp_beta(rmask_vector(kk), cmask_vector(kk)) = (B2'*mtx_var*B2) * c_corr_beta;
                end
            end
            index_mask = find(c_cov_interp_beta ~= 0);
            kappa = sqrt(kappa);
            sum_kappa = sum(kappa(:));
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
            %contrast_info = ['cinterp_SPM_nirs_' name '_' SPM_nirs.nirs.Hb '_Contrast_' contrast_name{Ic} '_View_' spec_hemi];
            %contrast_info_for_fig = ['cinterp SPM nirs ' name ' ' SPM_nirs.nirs.Hb ' Contrast ' contrast_name{Ic} ' View ' spec_hemi];
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

            cinterp_SPM_nirs.cbeta = c_interp_beta(:)';
            cinterp_SPM_nirs.ccov_beta = c_cov_interp_beta(:)';
            cinterp_SPM_nirs.s1 = s1;
            cinterp_SPM_nirs.s2 = s2;
            cinterp_SPM_nirs.sum_kappa = sum_kappa;
            cinterp_SPM_nirs.xCon = xCon;
            cinterp_SPM_nirs.Ic = Ic;
            cinterp_SPM_nirs.fname_ch = fname_ch; %%% file name (channel)
            cinterp_SPM_nirs.fname_SPM = fname_SPM; %%% file name (SPM)
            cinterp_SPM_nirs.viewBrain = spec_hemi; %%% view of the brain

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
            %end
            %fh1 = imcapture(gcf, 'all', 150);
            
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

            %out_T_map{Ic} =  T_map;
            %out_T_brain{Ic} = T_brain;
            %out_T_brain_over{Ic} = T_brain_over;
            %disp('done.');
            clear T_map T_brain T_brain_over kappa sum_kappa;
            clear cinterp_SPM_nirs cinterp_beta cinterp_cov_beta;  
        end %end for Ic %PP some of these quantities might not need to be recalculated for each contrast!!!
        catch
        end
    end %end for f1
end %end for v1
%save(job.NIRSmat{1,1},'NIRS');
%out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
out = [];
end


    
