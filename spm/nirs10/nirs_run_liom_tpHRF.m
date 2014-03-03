function out = nirs_run_liom_tpHRF(job)
%==========================================================================
%Ke Peng, LIOM, Polytechnique
%2013-12-17, version 0.1, Function created
%==========================================================================

for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
        Topodatafile = NIRS.Dt.ana.rend;
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        sess = job.session_select;
        onsets = NIRS.Dt.fir.Sess(sess).U;
        [new_path,name,ext] = fileparts(newNIRSlocation);
        %******************************************************************
        %Time-course
        %Data_dir = NIRS.Dt.s.p;% 
        %Data_dir = [NIRS.Dt.s.p '\Stat'];
        Data_dir = [NIRS.Dt.s.p '\dataSPM\coreg\Stat'];
        SessionTimeCourse = fopen_NIR(fullfile(Data_dir, ['Sess' int2str(sess) '.nir']), NC*3/2); %Load all time course of the session
        chromophore = job.chromophore_select;
        hb = get_chromophore(chromophore);
        if isfield(job.HRFmethod, 'simple_averaging')
            time_before_spk = job.HRFmethod.simple_averaging.HRFstart;
            time_after_spk = job.HRFmethod.simple_averaging.HRFend;
            EvtofInterest = job.HRFmethod.simple_averaging.EvtInterest;
            EvtofInterest_name = job.HRFmethod.simple_averaging.EvtInterest_name;
            if isfield(job.HRFmethod.simple_averaging.TCavoidance, 'avoidance_off')
                av.enabled = 0;
            else
                av.enabled = 1;
                if isfield(job.HRFmethod.simple_averaging.TCavoidance, 'avoidance_nonintere_only')
                    av.type = 1;
                    av.EvtofNonInterest_sz = job.HRFmethod.simple_averaging.TCavoidance.avoidance_nonintere_only.EvtNonInterest_sz;
                    av.EvtofNonInterest_spk = job.HRFmethod.simple_averaging.TCavoidance.avoidance_nonintere_only.EvtNonInterest_spk;
                elseif isfield(job.HRFmethod.simple_averaging.TCavoidance, 'avoidance_intern_only')
                    av.type = 2;
                elseif isfield(job.HRFmethod.simple_averaging.TCavoidance, 'avoidance_all')
                    av.type = 3;
                    av.EvtofNonInterest_sz = job.HRFmethod.simple_averaging.TCavoidance.avoidance_all.EvtNonInterest_sz; 
                    av.EvtofNonInterest_spk = job.HRFmethod.simple_averaging.TCavoidance.avoidance_all.EvtNonInterest_spk;
                end
            end
            [estiHRF estiHRF_std] = average_HRF_timecourse(SessionTimeCourse, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, EvtofInterest_name, chromophore, av, job.HRFmethod.simple_averaging.HRFsubtraction);
        elseif isfield(job.HRFmethod, 'deconvolution')
            time_before_spk = job.HRFmethod.deconvolution.HRFstart;
            time_after_spk = job.HRFmethod.deconvolution.HRFend;
            EvtofInterest = job.HRFmethod.deconvolution.EvtInterest;
            EvtofInterest_name = job.HRFmethod.deconvolution.EvtInterest_name;
            decov.poly_order = job.HRFmethod.deconvolution.poly_order;
            if isfield(job.HRFmethod.deconvolution.deriv_kernel, 'ols')
                decov.deriv_kernel = 'ols';
                decov.kl = 1;
            elseif isfield(job.HRFmethod.deconvolution.deriv_kernel, 'ridge')
                decov.deriv_kernel = 'ridge';
                decov.kl = 2;
                decov.k_val = job.HRFmethod.deconvolution.deriv_kernel.ridge.k_val;
                decov.rv_file = fullfile(new_path, ['S' int2str(sess) '_' hb '_trRV.mat']);
            end
            [estiHRF dm] = deconvolve_HRF_timecourse(SessionTimeCourse, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, EvtofInterest_name, chromophore, decov);
        end
        %******************************************************************
        
        %******************************************************************
        %Output
        if isfield(job.Outputformat, 'channel_view')
            sample_interval = job.Outputformat.channel_view.sample_interval;
            channel_no = job.Outputformat.channel_view.channel_no;
            estiHRF_l = size(estiHRF,2);
            [t0 t0_p] = sample_idx(sample_interval, estiHRF_l, fs);
            y = estiHRF(channel_no,t0_p);
            if isfield(job.HRFmethod, 'simple_averaging')
                figure_name = ['Avg_S' int2str(sess) '_C' int2str(channel_no) '_' hb];
                ebar_id = estiHRF_std(((channel_no-1)*estiHRF_l+1) : (channel_no*estiHRF_l));
                ebar = ebar_id(t0_p)';
            elseif isfield(job.HRFmethod, 'deconvolution')
                figure_name = ['Dec_S' int2str(sess) '_C' int2str(channel_no) '_' hb];
                row = 1 : size(dm.Z,2);
                rv_file = fullfile(new_path, ['S' int2str(sess) '_' hb '_trRV.mat']);
                ci = confidence_interval(1, dm.X, dm.ResSS, fs, row, channel_no, rv_file);
                ebar = ci(t0_p)';
            end
            fh = figure('Visible','on','Name',figure_name,'NumberTitle','off');
            t = t0 - abs(time_before_spk);
            errorbar(t,y,ebar,'b');
            xlabel('Time (s)');
            ylabel('\DeltaC (a.u.)');
            set(gca,'xlim',[time_before_spk,time_after_spk]);
            filen = fullfile(new_path,[figure_name '.fig']);
            print(fh, '-dpng', filen,'-r300');
            d = [t;y;ebar];
            fidt = fullfile(new_path,[figure_name '.mat']);
            save(fidt, 'd', '-mat');
            filen2 = fullfile(new_path,[figure_name '.png']);
            print(fh, '-dpng', filen2,'-r300');
            pause(1);
            close(fh);
        elseif isfield(job.Outputformat, 'intepolated_view')
            %******************************************************************
            %Interpolation kernel
            load(Topodatafile);
            brain_view = job.Outputformat.intepolated_view.brain_view;
            [side_hemi spec_hemi] = nirs_get_brain_view(brain_view);
            rchn = rendered_MNI{brain_view}.rchn;
            cchn = rendered_MNI{brain_view}.cchn;
            %Two options
            W.AllowExtrapolation = job.Outputformat.intepolated_view.AllowExtrapolation; %Option: 0: do not extrapolate
            W.no_interpolation = job.Outputformat.intepolated_view.no_interpolation; %Option: 1: do not interpolate
            %find channels which are visible from this projection view
            W.index_ch = find(rchn ~= -1);
            %rendering surface
            brain = rendered_MNI{brain_view}.ren;
            W.s1 = size(brain, 1);
            W.s2 = size(brain, 2);
            %split into HbO and HbR interpolations
            W.ch = W.index_ch;
            W.rchn = rchn(W.index_ch);
            W.cchn = cchn(W.index_ch);
            %Generate interpolation matrix
            Q = interpolation_kernel_cine_simplified(W);
            estiHRF = estiHRF(W.ch,:);
            %******************************************************************
            %Rendering
            sample_interval = job.Outputformat.intepolated_view.sample_interval;
            estiHRF_l = size(estiHRF,2);
            [t0 t0_p] = sample_idx(sample_interval, estiHRF_l, fs);
            thz = job.Outputformat.intepolated_view.threshold_value;
            if isfield(job.HRFmethod, 'simple_averaging')
                figure_name = ['Avg_S' int2str(sess) '_' hb '_' spec_hemi];
            elseif isfield(job.HRFmethod, 'deconvolution')
                figure_name = ['Dec_S' int2str(sess) '_' hb '_' spec_hemi];
            end
            [W, Z, F] = render_config(estiHRF, W, rendered_MNI, side_hemi, spec_hemi, new_path, brain_view, thz, figure_name);
            newsavepath = [F.pathn '\sess' int2str(sess)];
            if ~exist(newsavepath, 'dir')
                 mkdir(F.pathn,['sess' int2str(sess)]);
            end
            F.pathns = newsavepath;
            newsavepath_fig = [newsavepath '\fig'];
            if ~exist(newsavepath_fig, 'dir')
                 mkdir(newsavepath,'fig');
            end
            F.pathnsfig = newsavepath_fig;
            for i0 = 1 : length(t0_p)
                cp = t0_p(i0);
                ct = t0(i0);
                map = estiHRF(:,cp)';
                map_it = map * Q.B;
                map_it_reshaped = zeros(W.s1,W.s2);
                map_it_reshaped(Q.index_mask) = map_it(1,:);
                ct = ct-abs(time_before_spk);
                F.contrast_info_both = [int2str(ct) 's'];
                F.contrast_info_for_fig = F.contrast_info_both;
                F.s_map = map_it_reshaped;
                DF = nirs_render_map(F,W,Z);
                %Save figure
                if isfield(DF, 'fh2')
                    set(DF.fh2, 'visible', 'on');
                    text(200, 450, F.contrast_info_both, 'Color', 'r', 'FontSize',18);
                    DF.fh2 = save_figure(DF.fh2, Z, F); 
                    M(i0) = getframe(DF.fh2);
                    pause(0.5);
                    set(DF.fh2, 'visible', 'off');
                elseif isfield(DF, 'fh1')
                    set(DF.fh1, 'visible', 'on');
                    text(200, 450, F.contrast_info_both, 'Color', 'r', 'FontSize',18);
                    DF.fh1 = save_figure(DF.fh1, Z, F);
                    M(i0) = getframe(DF.fh1);
                    pause(0.5);
                    set(DF.fh1, 'visible', 'off');                   
                end
            end
            forename = fullfile(F.pathn, figure_name);
            movie2avi(M,[forename,'.avi'],'FPS', 1/sample_interval, 'compression', 'None');
            %movie2avi(mov, 'avifile.avi', 'compression', 'Cinepak');
        end       
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1));
        disp(['Could not estimate HRF for subject' int2str(iSubj) ' for ' job.NIRSmat{iSubj,1}]);
    end
end
outNIRSmat = NIRS;
out.NIRSmat = outNIRSmat;

function [W, Z, F] = render_config(estiHRF, W, rendered_MNI0, side_hemi, spec_hemi, new_path, brain_view, thz, figure_name)
W.side_hemi = side_hemi;
W.spec_hemi = spec_hemi;
%View dependent info for figures
brain = rendered_MNI0{brain_view}.ren;
msk = brain>1;brain(msk)=1;
msk = brain<0;brain(msk)=0;
brain = brain * 0.5;
W.brain = brain;
%For single subject group of sessions analysis
if isfield(rendered_MNI0{W.side_hemi},'view_mask_2d')
    W.brain_view_mask_2d = rendered_MNI0{W.side_hemi}.view_mask_2d;
end
W.thz = thz;
F.pathn = new_path;
load Split
F.split = split;
F.tstr = 'T';
F.contrast_info_both = 'test4';
F.contrast_info_both_for_fig = 'test5';
Z.GroupColorbars = 0;
F.str_cor = figure_name;
F.contrast_info = figure_name;
F.contrast_info_for_fig = 'test3';
Z.gen_fig = 0;
Z.gen_tiff = 0;
map_p_idx = find(estiHRF > 0);
map_n_idx = find(estiHRF < 0);
map_c_min = min(estiHRF(map_p_idx));
map_c_max = max(estiHRF(map_p_idx));
map_c_min2 = min(estiHRF(map_n_idx));
map_c_max2 = max(estiHRF(map_n_idx));
c_min = max(map_c_min, thz);
c_max = map_c_max;
c_min2 = map_c_min2;
c_max2 = min(map_c_max2, -thz);
Z.cbar.c_min = floor(c_min*10)/10;
Z.cbar.c_max = ceil(c_max*10)/10;
Z.cbar.c_min2 = floor(c_min2*10)/10;
Z.cbar.c_max2 = ceil(c_max2*10)/10;
Z.cbar.colorbar_override = 1;
Z.cbar.visible = 'off';
Z.write_neg_pos = 0;

function [t0 t0_p] = sample_idx(sample_interval, estiHRF_l, fs)
t = floor(estiHRF_l/fs);
t0 = 0:sample_interval:t;
t0_p = round(t0*fs) + 1; %The first point is 1 instead of 0
t0_p = t0_p(find(t0_p <= estiHRF_l));
%estiHRF_rep = repmat(estiHRF_c,t0_l,estiHRF_l);
%t0_rep = repmat(t0_p',t0_l,estiHRF_l);
%[idx vmin] = min(abs(estiHRF_rep-t0_rep),[],2);

function fh = save_figure(fh, Z, F)
if ~Z.gen_fig %If not saved before
    filen1 = fullfile(F.pathnsfig,[F.contrast_info '_' F.contrast_info_both '.fig']);
    saveas(fh,filen1,'fig');
end
if ~Z.gen_tiff
    filen2 = fullfile(F.pathns,[F.contrast_info '_' F.contrast_info_both '.png']);
    print(fh, '-dpng', filen2,'-r300');
end

function ci = confidence_interval(pre_coloring, X, ResSS, fs, row, ch, rv_filen)
%calculate confidence intervals for each beta
alpha = 0.05; %level of significance
if pre_coloring
    %calculate standard deviation of beta
    try
        load(rv_filen);
    catch
        HParam.type = 'none';
        LParam.type = 'hrf';
        K = struct( 'HParam', HParam,...
            'RT', 1/fs ,...
            'row', row,...
            'LParam', LParam);
        K = spm_filter_HPF_LPF_WMDL(K);
        S = K.KL; %Tempoal autocorrelation kernel
        clear K
        pKX = pinv(X'*X) * X'; %pKX = pinv(X'X)*X'
        Bcov = (pKX * S);
        Bcov = Bcov * Bcov';%Covariance matrix of beta
        %calculate effective DOF
        Var1 = S * S'; %Var1 = V;
        clear S
        Var2 = Var1 - X * (pKX * Var1); %R = I - X*pinv(X'X)X'; Var2 = RV;
        Var3 = Var1 * Var2; %Var3 = VRV;
        clear Var1
        Var4 = Var3 - X * (pKX * Var3);  %Var4 = RVRV;
        clear Var3
        trRV = sum(diag(Var2));
        trRVRV = sum(diag(Var4));
        rv.Bcov = Bcov;
        rv.trRV = trRV;
        rv.trRVRV = trRVRV;
        save(rv_filen, '-struct', 'rv');
    end
    % Confidence interval
    try
        erdf = trRV^2/trRVRV; %effective degrees of freedom
        tval = tinv(alpha/2, erdf); %corresponding t-value
        ervar = ResSS/trRV;%Usual estimator of variance of error of betas in a least-square scheme
        Bse = sqrt(diag(Bcov) .* ervar(ch));%Standard deviation of beta of one channel, (correction to the usual estimator, Ye et al. 2009)
        ci = abs(tval*Bse); %t_value * SE(beta). Formula to calulcate the CI: beta +/- t_value*SE(beta);
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp('Problem calculating degrees of freedom');
    end
end

%         sample_interval = job.sample_interval;
%         tp_interval = round(sample_interval*fs);
%         cp = 1;
%         i0 = 0;
%         while cp <= size(estiHRF,2)
%             map = estiHRF(:,cp)';
%             map_it = map * Q.B;
%             map_it_reshaped = zeros(W.s1,W.s2);
%             %amplitude_i_reshaped = reshape(amplitude_i,W.s1,W.s2);
%             map_it_reshaped(Q.index_mask) = map_it(1,:);
%             i0 = i0+1;
%             if 0
%                 h = figure;
%                 c = max(max(abs(estiHRF)));
%                 imagesc(amplitude_i_reshaped);
%                 colorbar
%                 caxis([-c c]);
%                 %Write markers
%                 if isfield(job.HRFmethod, 'simple_averaging') || isfield(job.HRFmethod, 'deconvolution')
%                     mk = [int2str(time_before_spk + (i0-1)*sample_interval) ' s'];
%                     text(300, 450, mk, 'Color', 'b', 'FontSize',18);
%                 end
%                 %save
%                 switch chromophore
%                     case 0
%                         filestr = [new_path '\S' int2str(sess) '_' spec_hemi '_HbO_' int2str(i0)];
%                     case 1
%                         filestr = [new_path '\S' int2str(sess) '_' spec_hemi '_HbR_' int2str(i0)];
%                     case 2
%                         filestr = [new_path '\S' int2str(sess) '_' spec_hemi '_HbT_' int2str(i0)];
%                 end
%                 saveas(gcf,filestr,'png');
%                 pause(1);
%                 close(h);
%             else
%                 [W, Z, F] = render_config(estiHRF, W, rendered_MNI0, side_hemi, spec_hemi, new_path, brain_view);
%                 F.s_map = map_it_reshaped;
%                 DF = nirs_render_map(F,W,Z);
%             end
%             cp = cp + tp_interval;
%         end
%         %******************************************************************