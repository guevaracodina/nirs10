function out = nirs_run_extract_map_data(job)
%Views: 
views_to_run = job.view;
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'
try
    job.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_max_HbR;
    select_mode = 1;
catch
    try
        job.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_max_all;
        select_mode = 2;
    catch
        try
            job.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_coordinates;
            coord_min = job.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_coordinates.extract_coord_min;
            coord_max = job.extract_select_mode.extract_auto_mode.extract_select_auto_mode.extract_coordinates.extract_coord_max;
            select_mode = 3;
        catch
            try
                job.extract_select_mode.extract_manual_mode;
                select_mode = 0;
            catch                
                select_mode = 0;
            end
        end
    end
end
try
    extract_struct_name = job.extract_struct_name;
catch
    extract_struct_name = 'ED';
end
try 
    extract_base_contrast = job.extract_base_contrast;
catch
    extract_base_contrast = 1;
end
try
    extract_contrast = job.extract_contrast;
catch
    extract_contrast = 1;
end
try 
    %job.extract_average_mode.extract_radius;
    radius = job.extract_average_mode.extract_radius.extract_radius_val;
    extract_mode = 0;
catch
    try
        %job.extract_average_mode.extract_threshold;
        threshold = job.extract_average_mode.extract_threshold.extract_threshold_val;
        radius = job.extract_average_mode.extract_threshold.extract_radius_val;
        extract_mode = 1;
    catch
        extract_mode = 0;
    end
end
    
%Store all extracted data in ED structure
ED.job = job; 
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
        %Load TOPO
        try 
            ftopo = NIRS.TOPO; 
        catch
            ftopo = fullfile(dir1,'TOPO.mat');
        end
        TOPO = [];
        load(ftopo); 

        %Big loop over views 
        for v1=1:size(views_to_run,2)
            try
                brain_view = views_to_run(v1);
                try 
                    TOPO.v{brain_view};
                    view_estimated = 1;
                catch
                    view_estimated = 0;
                end
                if view_estimated
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
                    %Structure for passing GLM and interpolation data
                    W = [];
                    % channel information
                    rchn = rendered_MNI{side_hemi}.rchn;
                    cchn = rendered_MNI{side_hemi}.cchn;
                    %rendering surface
                    brain = rendered_MNI{side_hemi}.ren;
                    W.brain = brain * 0.5;
                    W.s1 = size(brain, 1);
                    W.s2 = size(brain, 2);
                    %find channels which are visible from this projection view
                    W.index_ch = find(rchn ~= -1);
                    if isempty(W.index_ch)
                        TOPO.v{side_hemi}.Warning = 'No channel found for this view';
                        disp(['No channel for view ' int2str(side_hemi) ': Probable coregistration problem. Skipping this view']);
                    else
                        %split into HbO and HbR interpolations
                        wl = NIRS.Cf.dev.wl;
                        %for NIRS acquisition with 2 wavelengths only
                        
                        % MD: it seems to me that the channels corresponding to
                        % HbO/HbR should be fixed after the module
                        % nirs_run_ODtoHbOHbR ? Not sure this gives the right
                        % order....?!?
                        
                        if wl(1) > 750
                            %first wavelength is "HbO-like"
                            W.ch_HbO = W.index_ch;
                            W.ch_HbR = NC/2 + W.index_ch;
                        else
                            W.ch_HbR = W.index_ch;
                            W.ch_HbO = NC/2 + W.index_ch;
                        end
                        try
                            W.ch_HbT = NC+W.index_ch;
                        end
                        %nch = length(index_ch);
                        W.rchn = rchn(W.index_ch);
                        W.cchn = cchn(W.index_ch);
                        W.spec_hemi = spec_hemi;
                        W.side_hemi = side_hemi;
                        
                        
                        beta_tmpO = [];
                        beta_tmpR = [];
                        beta_tmpT = [];
                        cHb{1} = W.ch_HbO;
                        cHb{2} = W.ch_HbR;
                        cHb{3} = W.ch_HbT;
                        for f1 = 1:length(SPM.xXn)
                            %beta for all sessions, with regressors in rows and channels in columns
                            beta_tmpO = [beta_tmpO; SPM.xXn{f1}.beta(:, cHb{1})];
                            beta_tmpR = [beta_tmpR; SPM.xXn{f1}.beta(:, cHb{2})];
                            try
                                beta_tmpT = [beta_tmpT; SPM.xXn{f1}.beta(:, cHb{3})];
                            end
                        end
                        
                        Q = interpolation_kernel_short(length(cHb{1}),W);
                        %TO DO
                        %Big loop over selected contrasts
                        base_con = extract_base_contrast;
                        for c1 = 1:length(extract_contrast)
                            try
                            con = extract_contrast(c1);
                            for h1=1:3 %careful if HbT not generated
                                %Select point(s) on map
                                switch select_mode
                                    case 0 %manual
                                        if h1 == 1 && c1 == 1 %only go through once                                           
                                            m = TOPO.v{side_hemi}.group.hb{2}.c{base_con}.Tmap;
                                            disp(['Displayed HbR map for contrast ' int2str(base_con)]);
                                            M = get_min_max(m);
                                            figure; imagesc(m);
                                            M.rmin = spm_input('Enter rmin ',1);
                                            M.cmin = spm_input('Enter cmin ',1);
                                            M.rmax = spm_input('Enter rmax ',1);
                                            M.cmax = spm_input('Enter cmax ',1);
                                        end
                                    case 1 %extract_max_HbR - might be recalculating
                                        m = TOPO.v{side_hemi}.group.hb{2}.c{base_con}.Tmap;
                                        M = get_min_max(m);
                                    case 2 %extract_max_all
                                        %might want base_con here?
                                        m = TOPO.v{side_hemi}.group.hb{h1}.c{con}.Tmap;
                                        M = get_min_max(m); %
                                    case 3 %extract_coordinates
                                        M.rmin = coord_min(1); M.cmin = coord_min(2);
                                        M.rmax = coord_max(1); M.cmax = coord_max(2);
                                    otherwise
                                end
                                lmin = []; %list of pixels for min and max clusters
                                lmax = [];
                                %find pixels to average
                                switch extract_mode
                                    case 0 %radius
                                        %double loop x & y to span a circle
                                        for x=-radius:radius
                                            for y=-radius:radius
                                                if x*x+y*y <= radius*radius
                                                    %add pixel to list of pixels to average
                                                    %but also check that the
                                                    %pixel is on the map
                                                    if ~(brain(M.rmax+x,M.cmax+y)==0)
                                                        %lmax = [lmax;[x y]];
                                                        lmax = [lmax;[M.rmax+x M.cmax+y]];
                                                    end
                                                    if ~(brain(M.rmin+x,M.cmin+y)==0)
                                                        %lmin = [lmin;[x y]];
                                                        lmin = [lmin;[M.rmin+x M.cmin+y]];
                                                    end
                                                end
                                            end
                                        end
                                    case 1 %threshold, within a set radius
                                        %double loop x & y to span a circle
                                        for x=-radius:radius
                                            for y=-radius:radius
                                                if x*x+y*y <= radius*radius
                                                    %add pixel to list of pixels to average
                                                    %but also check that the
                                                    %pixel is on the map
                                                    if ~(brain(M.rmax+x,M.cmax+y)==0) && m(M.rmax+x,M.cmax+y) > threshold
                                                        %lmax = [lmax;[x y]];
                                                        lmax = [lmax;[M.rmax+x M.cmax+y]];
                                                    end
                                                    if ~(brain(M.rmin+x,M.cmin+y)==0) && m(M.rmin+x,M.cmin+y) < -threshold
                                                        %lmin = [lmin;[x y]];
                                                        lmin = [lmin;[M.rmin+x M.cmin+y]];
                                                    end
                                                end
                                            end
                                        end
                                    otherwise
                                end
                                ED.v{side_hemi}.hb{h1}.lmax = lmax;
                                ED.v{side_hemi}.hb{h1}.lmin = lmin;
                                ED.v{side_hemi}.hb{h1}.M = M;
                                %Extract time series at selected points
                                %from SPM.xY.P -- does that exist for group
                                %studies?
                                if c1 == 1
                                    for f1=1:length(SPM.xY.P) %should not need to reopen so often
                                        try
                                            %filtered data
                                            f0 = SPM.xY.Pf{f1};
                                            Y = fopen_NIR(f0,SPM.xY.Cf)';
                                            filteredOK = 1;
                                            %HbTrun = 1;
                                        catch
                                            f0 = SPM.xY.P{f1};
                                            Y = fopen_NIR(f0,NIRS.Cf.H.C.N)';
                                            filteredOK = 0;
                                            %HbTrun = 0;
                                        end
                                        try %for when HbT was not run
                                            ED.v{side_hemi}.s{f1}.hb{h1}.Ymin = interp_series(Y(:,cHb{h1}),lmin,Q);
                                            ED.v{side_hemi}.s{f1}.hb{h1}.Ymax = interp_series(Y(:,cHb{h1}),lmax,Q);
                                            ED.v{side_hemi}.s{f1}.hb{h1}.Ymin_Sigma = std(ED.v{side_hemi}.s{f1}.hb{h1}.Ymin);
                                            ED.v{side_hemi}.s{f1}.hb{h1}.Ymax_Sigma = std(ED.v{side_hemi}.s{f1}.hb{h1}.Ymax);
                                            ED.v{side_hemi}.s{f1}.hb{h1}.filteredOK = filteredOK;
                                        end
                                    end
                                    %group Std
                                    try
                                        Ymin = [];
                                        Ymax = [];
                                        for f1=1:length(SPM.xY.P)
                                            Ymin = [Ymin; ED.v{side_hemi}.s{f1}.hb{h1}.Ymin];
                                            Ymax = [Ymax; ED.v{side_hemi}.s{f1}.hb{h1}.Ymax];
                                        end
                                        if length(SPM.xY.P) > 1
                                            ED.v{side_hemi}.g.hb{h1}.Ymin = Ymin;
                                            ED.v{side_hemi}.g.hb{h1}.Ymax = Ymax;
                                            ED.v{side_hemi}.g.hb{h1}.Ymin_Sigma = std(Ymin);
                                            ED.v{side_hemi}.g.hb{h1}.Ymax_Sigma = std(Ymax);
                                            ED.v{side_hemi}.g.hb{h1}.filteredOK = filteredOK;
                                        end
                                    end
                                end
                                %Extract statistics data from maps
                                
                                for f1=1:length(SPM.xXn)
                                    try
                                    if c1 ==1
                                        ED.v{side_hemi}.s{f1}.hb{h1}.bmin = interp_series(SPM.xXn{f1}.beta(:,cHb{h1}),lmin,Q);
                                        ED.v{side_hemi}.s{f1}.hb{h1}.bmax = interp_series(SPM.xXn{f1}.beta(:,cHb{h1}),lmax,Q);
                                        %normalized by standard deviation
                                        ED.v{side_hemi}.s{f1}.hb{h1}.bnmin = ED.v{side_hemi}.s{f1}.hb{h1}.bmin/ED.v{side_hemi}.s{f1}.hb{h1}.Ymin_Sigma;
                                        ED.v{side_hemi}.s{f1}.hb{h1}.bnmax = ED.v{side_hemi}.s{f1}.hb{h1}.bmax/ED.v{side_hemi}.s{f1}.hb{h1}.Ymax_Sigma;
                                    end
                                    %for each contrast
                                    c0 = TOPO.SSxCon(c1).c;
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bmin = c0'*ED.v{side_hemi}.s{f1}.hb{h1}.bmin;
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bmax = c0'*ED.v{side_hemi}.s{f1}.hb{h1}.bmax;
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bnmin = c0'*ED.v{side_hemi}.s{f1}.hb{h1}.bnmin;
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bnmax = c0'*ED.v{side_hemi}.s{f1}.hb{h1}.bnmax;
                                    %should be the same as:
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bcmin = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_interp_beta(c1,:,:)),lmin,[]);
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bcmax = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_interp_beta(c1,:,:)),lmax,[]);
                                    %normalized by standard deviation
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bNmin = ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bcmin/ED.v{side_hemi}.s{f1}.hb{h1}.Ymin_Sigma;
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bNmax = ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bcmax/ED.v{side_hemi}.s{f1}.hb{h1}.Ymax_Sigma;
                                    %interpolated covariance
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.covmin = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_cov_interp_beta(c1,:,:)),lmin,[]);
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.covmax = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_cov_interp_beta(c1,:,:)),lmax,[]);
                                    %interpolated F
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.Fmin = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_interp_F(c1,:,:)),lmin,[]);
                                    ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.Fmax = interp_series(squeeze(TOPO.v{side_hemi}.s{f1}.hb{h1}.c_interp_F(c1,:,:)),lmax,[]);
                                    end
                                    
                                    %                                 W.var = SPM.xX.var; %careful, var can be a Matlab function,
                                    %                                 %but instead we want W.var
                                    %                                 W.beta_HbO = beta_tmpO(:); %taken as one vector
                                    %                                 W.beta_HbR = beta_tmpR(:); %taken as one vector
                                    %                                 W.mtx_var_HbO = diag(W.var(W.ch_HbO));
                                    %                                 W.mtx_var_HbR = diag(W.var(W.ch_HbR));
                                    %                                 try
                                    %                                     W.beta_HbT = beta_tmpT(:); %taken as one vector
                                    %                                     W.mtx_var_HbT = diag(W.var(W.ch_HbT));
                                    %                                 end
                                    %
                                    %                                 W.corr_beta = SPM.xX.corr_beta;
                                    %                                 [TOPO] = extract_data_core(Z,W,TOPO,SPM.xXn{f1},nCon,f1);
                                    %
                                end %end for f1
                                %Group of sessions
                                if isfield(TOPO.v{side_hemi},'g')
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.tmin = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.tmap_group(:,:)),lmin,[]);
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.tmax = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.tmap_group(:,:)),lmax,[]);
                                    
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmin = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.beta_group(:,:)),lmin,[]);
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmax = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.beta_group(:,:)),lmax,[]);
                                    %normalized by standard deviation
                                    %ED.v{side_hemi}.g.hb{h1}.c{c1}.bNmin = ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmin/ED.v{side_hemi}.g.hb{h1}.Ymin_Sigma;
                                    %ED.v{side_hemi}.g.hb{h1}.c{c1}.bNmax = ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmax/ED.v{side_hemi}.g.hb{h1}.Ymax_Sigma;
                                    %interpolated covariance
                                    try
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.stdmin = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.std_group(:,:)),lmin,[]);
                                    ED.v{side_hemi}.g.hb{h1}.c{c1}.stdmax = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c{2*c1-1}.std_group(:,:)),lmax,[]);
                                    end
                                    %interpolated F
                                    %ED.v{side_hemi}.g.hb{h1}.c{c1}.Fmin = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c_interp_F(c1,:,:)),lmin,[]);
                                    %ED.v{side_hemi}.g.hb{h1}.c{c1}.Fmax = interp_series(squeeze(TOPO.v{side_hemi}.g.hb{h1}.c_interp_F(c1,:,:)),lmax,[]);

                                else
                                    if isfield(TOPO.v{side_hemi},'group')
                                        %group study proper
                                        ED.v{side_hemi}.group.hb{h1}.c{c1}.bcmin = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c{2*c1-1}.beta_group(:,:)),lmin,[]);
                                        ED.v{side_hemi}.group.hb{h1}.c{c1}.bcmax = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c{2*c1-1}.beta_group(:,:)),lmax,[]);
                                        %normalized by standard deviation
                                        %ED.v{side_hemi}.group.hb{h1}.c{c1}.bNmin = ED.v{side_hemi}.group.hb{h1}.c{c1}.bcmin/ED.v{side_hemi}.group.hb{h1}.Ymin_Sigma;
                                        %ED.v{side_hemi}.group.hb{h1}.c{c1}.bNmax = ED.v{side_hemi}.group.hb{h1}.c{c1}.bcmax/ED.v{side_hemi}.group.hb{h1}.Ymax_Sigma;
                                        %interpolated covariance
                                        try
                                            ED.v{side_hemi}.group.hb{h1}.c{c1}.stdmin = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c{2*c1-1}.std_group(:,:)),lmin,[]);
                                            ED.v{side_hemi}.group.hb{h1}.c{c1}.stdmax = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c{2*c1-1}.std_group(:,:)),lmax,[]);
                                        end
                                        %interpolated F
                                        %ED.v{side_hemi}.group.hb{h1}.c{c1}.Fmin = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c_interp_F(c1,:,:)),lmin,[]);
                                        %ED.v{side_hemi}.group.hb{h1}.c{c1}.Fmax = interp_series(squeeze(TOPO.v{side_hemi}.group.hb{h1}.c_interp_F(c1,:,:)),lmax,[]);
                                    end
                                end
                            end %end for h1
                        end %end for c1
                        end
                        %Organize data in a more convenient way
                        try                           
                            for h1=1:3
                                c_min = [];
                                c_max = [];
                                for c1=1:length(ED.v{side_hemi}.s{1}.hb{1}.c)
                                    tmp_min = [];
                                    tmp_max = [];
                                    tmp_std_min = [];
                                    tmp_std_max = [];
                                    tmin = [];
                                    tmax = [];
                                    for f1=1:length(SPM.xXn)
                                        tmp_min = [tmp_min ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bNmin];
                                        tmp_max = [tmp_max ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.bNmax];
                                        try
                                            tmp_std_min = [tmp_std_min (ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.covmin).^0.5];
                                            tmp_std_max = [tmp_std_max (ED.v{side_hemi}.s{f1}.hb{h1}.c{c1}.covmax).^0.5];
                                        end
                                    end
%                                    try
                                        %add group result to list
                                        if isfield(TOPO.v{side_hemi},'g')
                                            tmp_min = [tmp_min sum(tmp_min)/length(tmp_min) ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmin];
                                            tmp_max = [tmp_max sum(tmp_max)/length(tmp_max) ED.v{side_hemi}.g.hb{h1}.c{c1}.bcmax];
                                            tmp_std_min = [tmp_std_min std(tmp_min)/length(tmp_min)^0.5 ED.v{side_hemi}.g.hb{h1}.c{c1}.stdmin];
                                            tmp_std_max = [tmp_std_max std(tmp_max)/length(tmp_max)^0.5 ED.v{side_hemi}.g.hb{h1}.c{c1}.stdmax];
                                            tmin = ED.v{side_hemi}.g.hb{h1}.c{c1}.tmin;
                                            tmax = ED.v{side_hemi}.g.hb{h1}.c{c1}.tmax;
                                        else
                                            if isfield(TOPO.v{side_hemi},'group')
                                                tmp_min = [tmp_min  sum(tmp_min)/length(tmp_min) ED.v{side_hemi}.group.hb{h1}.c{c1}.stdmin];
                                                tmp_max = [tmp_max  sum(tmp_max)/length(tmp_max) ED.v{side_hemi}.group.hb{h1}.c{c1}.stdmax];
                                                tmin = ED.v{side_hemi}.group.hb{h1}.c{c1}.tmin;
                                                tmax = ED.v{side_hemi}.group.hb{h1}.c{c1}.tmax;
                                            end
                                        end
 %                                   end
                                    c_min = [c_min;tmp_min];
                                    c_max = [c_max;tmp_max];
                                    epilepsy = 1;
                                    if epilepsy
                                        %add ratio 
                                        cmin(end+1,:) = c_min(1,:)./c_min(2,:);
                                        cmax(end+1,:) = c_max(1,:)./c_max(2,:);
                                    end
                                end
                                ED.v{side_hemi}.hb{h1}.bNmin = c_min;
                                ED.v{side_hemi}.hb{h1}.bNmax = c_max; 
                                ED.v{side_hemi}.hb{h1}.tmin = tmin;
                                ED.v{side_hemi}.hb{h1}.tmax = tmax;
                                ED.v{side_hemi}.hb{h1}.note = 'bNmin: normalized amplitudes by session or subject, then unweighted average and precision-weighted group result, for each contrast and perhaps ratio of first 2 contrasts'; 
                            end
                        end
                    end
                end
            catch exception
                disp(exception);
                disp(['Could not extract data for view ' spec_hemi ' for subject ' int2str(Idx)]);
            end
        end %end for v1
    catch exception
        disp(exception);
        disp(['Could not create extract data']);
    end
    outfile = fullfile(dir1,[extract_struct_name '.mat']);
    NIRS.ED = outfile;
    save(outfile,'ED');
    save(job.NIRSmat{Idx,1},'NIRS');
end
out.NIRSmat = job.NIRSmat;
end
 
function M = get_min_max(m)
    %min
    [M.min M.pmin] = min(m(:));
    [M.rmin M.cmin] = ind2sub(size(m),M.pmin);
    
    %min
    [M.max M.pmax] = max(m(:));
    [M.rmax M.cmax] = ind2sub(size(m),M.pmax);
end