function out = nirs_run_liom_intrasubject_average(job)
filter_vasomotion = 1;
%physiological confounds
if isfield(job.NIRSchannelsConfound,'NIRSconfounds')
    NIRSconfounds = job.NIRSchannelsConfound.NIRSconfounds;
    NIRSconfoundsOn = 1;
    NumChConfounds = NIRSconfounds.NumChConfounds;
    MinChDist = NIRSconfounds.MinChDist;
    MaxChDist = NIRSconfounds.MaxChDist;
else
    NIRSconfoundsOn = 0;
end
if isfield(job.vasomotion_choice,'vasomotion_on')
    vasomotion_on = 1;
    select_chromophore = job.vasomotion_choice.vasomotion_on.select_chromophore;
else
    vasomotion_on = 0;
end
%PCA - the functionality of this has not been tested
PCA = job.channel_pca;

%HPF - Butterworth infinite impulse response filter
if isfield(job.hpf_butter,'hpf_butter_On')
    hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
    HPFbutter = 1;
else
    if isfield(job.hpf_butter,'remove_linear')
        HPFbutter = 2;
        hpf_butter_freq = 0; %not used
        hpf_butter_order = 3; %not used
    else
        if isfield(job.hpf_butter,'GLM_remove_linear')
            HPFbutter = 3; %no high pass filter
            hpf_butter_freq = 0; %not used
            hpf_butter_order = 3; %not used
        else
            if isfield(job.hpf_butter,'SPM_cosine_filter')
                HPFbutter = 4; %no high pass filter
                hpf_butter_freq = 1/300; %1/240; %4 minute
                hpf_butter_order = 3; %not used
            else
                HPFbutter = 0; %no high pass filter
                hpf_butter_freq = 0; %not used
                hpf_butter_order = 3; %not used
            end
        end
    end
end

AvgFilters = job.AvgFilters;
if isfield(job.averaging_choice,'average_all_data')
    need_onsets = 0;
    U = [];
else
    need_onsets = 1;
end

%HPF - filter from NIRS_SPM - note that Butterworth HPF can be used with it
if isfield(AvgFilters.nirs_hpf,'hpf_dct')
    HPF = ['DCT, ' int2str(AvgFilters.nirs_hpf.hpf_dct.hpf_dct_cutoff)];
else
    if isfield(AvgFilters.nirs_hpf,'hpf_wavelet')
        HPF = ['wavelet,' int2str(AvgFilters.nirs_hpf.hpf_wavelet.hpf_wavelet_iter)];
    else
        if isfield(AvgFilters.nirs_hpf,'hpf_none')
            HPF = 'none';
        else
            disp('Unrecognized high pass filter');
        end
    end
end

%LPF - filter from NIRS_SPM
if isfield(AvgFilters.nirs_lpf,'lpf_gauss')
    FWHM = AvgFilters.nirs_lpf.lpf_gauss.fwhm1;
    LPF = 'gaussian';
else
    if isfield(AvgFilters.nirs_lpf,'lpf_hrf')
        LPF = 'hrf';
    else
        if isfield(AvgFilters.nirs_lpf,'lpf_none')
            LPF = 'none';
        else
            disp('Unrecognized low pass filter');
        end
    end
end

%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        %Objective is to fill SPM structure
        SPM = [];
        SPM.Idx = Idx;
        %always store SPM analysis in some directory
        if ~isfield(job.NIRSmatCopyChoice,'NIRSmatCopy')
            job.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = 'Avg';
        end
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'Avg_OK') || job.force_redo)
            [spm_dir dummy] = fileparts(newNIRSlocation);
            NIRS.spm_dir = spm_dir;
            %use last step of preprocessing
            lst = length(NIRS.Dt.fir.pp);
            rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
            NC = NIRS.Cf.H.C.N;
            fs = NIRS.Cf.dev.fs;
            % MICHÈLE 21 sept. 2011
            %nsess = size(rDtp,1);
            nsessAll = size(rDtp,1);
            % Sessions to use (assumed same for all subjects)
            try
                idx_sess = job.sessions;
                if max(idx_sess)>nsessAll || min(idx_sess)<1 ||...
                        length(idx_sess)>nsessAll || ... % wrong entry...
                        isempty(idx_sess) % or default value!
                    goToCatch;
                end
            catch
                idx_sess = 1:nsessAll;
            end
            nsess = length(idx_sess);
            if need_onsets
                %Find onsets
                try
                    for f=1:nsess
                        iSess = idx_sess(f);
                        %PLEASE DO NOT MODIFY THE NEXT TWO LINES!!!!!!!!!!!!!!!!
                        NIRS.Dt.fir.Sess(1).U(1).name;
                        if ~isempty(NIRS.Dt.fir.Sess(iSess).U(1).name)
                            SPM.Sess(f) = NIRS.Dt.fir.Sess(iSess);
                        else % no onsets
                            SPM.Sess(f).U = [];
                            SPM.Sess(f).C.C = [];
                            SPM.Sess(f).C.name = cell(1,0);
                        end
                    end
                catch
                    %Ignore parametric modulations - cf spm_run_fmri_design.m
                    P.name = 'none';
                    P.h    = 0;
                    for f=1:nsess
                        iSess = idx_sess(f);
                        try
                            %load onset file
                            clear names onsets durations
                            load(job.subj(1,1).input_onsets{f}); %careful, must have same onsets for all subjects
                            for kk = 1:size(names, 2)
                                SPM.Sess(f).U(kk).name = names(kk);
                                SPM.Sess(f).U(kk).ons = onsets{kk};
                                SPM.Sess(f).U(kk).dur = durations{kk};
                                SPM.Sess(f).U(kk).P = P;
                            end
                        catch
                            %Could not load onset
                            disp(['Could not find onsets - assuming baseline scan (no stimuli) on session ' int2str(iSess) '.']);
                            % MICHÈLE 21 sept. 2011 - for resting state scans one must
                            % be allowed to include 0 conditions in the design matrix
                            % (only other regressors).
                            SPM.Sess(f).U = [];
                            SPM.Sess(f).C.C = [];
                            SPM.Sess(f).C.name = cell(1,0);
                            % This way, spm_get_ons will not prompt the user to
                            % manually enter conditions as it does when the U field
                            % does not exist.
                        end
                    end
                end
            end
            %Adding confound regressors
            for f=1:nsess
                iSess = idx_sess(f);
                C = [];
                Cname = {};
                try
                    if job.GLM_include_cardiac
                        %heart rate regressor
                        C = NIRS.Dt.fir.Sess(iSess).fR{1};
                        C = C - repmat(mean(C),[length(C),1]);
                        Cname = {'H'};
                    end
                    if job.GLM_include_Mayer
                        %Mayer wave regressor
                        C = [C NIRS.Dt.fir.Sess(iSess).mR{1}];
                        Cname = [Cname {'M'}];
                    end
                    %wl = NIRS.Cf.dev.wl;
                    %HbO_like = [];
                    %for i=1:length(wl)
                    %    if wl(i) > 750 %in nanometer
                    %       %found a wavelength that is HbO-like
                    %       HbO_like = [HbO_like i];
                    %   end
                    %end
                    %HbO channels
                    %chHbO = NIRS.Cf.H.C.wl== HbO_like;
                    chHbO = 1:NC/2;
                    fullHbO = NIRS.Cf.H.C.gp(chHbO);
                    if vasomotion_on
                        %get data for that session
                        d = fopen_NIR(rDtp{iSess,1},NC);
                        switch select_chromophore
                            case 1 %HbT
                                tmpC = 2*mean(d,1)';
                            case 2 %HbR
                                tmpC = mean(d(logical(1-chHbO),:),1)';
                            case 3 %HbO
                                tmpC = mean(d(chHbO,:),1)';
                            case 4 %Hbo&HbR
                                tmpC = [mean(d(chHbO,:),1)' mean(d(logical(1-chHbO),:),1)'];
                        end
                        if filter_vasomotion
                            cutoff=0.04; %in Hz, 25 s
                            FilterOrder=3; %Is this too weak?
                            tmpC = ButterHPF(fs,cutoff,FilterOrder,tmpC);
                            cutoff=0.125; %PP 8 s
                            tmpC = ButterLPF(fs,cutoff,FilterOrder,tmpC);
                        end
                        switch select_chromophore
                            case {1 2 3}
                                C = [C tmpC];
                                Cname = [Cname {'V'}];
                            case 4
                                C = [C tmpC];
                                Cname = [Cname {'V_R'} {'V_O'}];
                        end
                    end
                    if ~isempty(job.subj.multi_reg)
                        %"Multiple regressors" file
                        %nb = 0;
                        %for iFile = 1:length(job.subj.multi_reg) % if more
                        %than one file per session... not implemented
                        %nb = nb+1;
                        try
                            [dir fil ext] = fileparts(job.subj.multi_reg{f});
                            if strcmp(ext,'.mat')
                                regressors = load(job.subj.multi_reg{f});
                                % not going to work, because I don't know
                                % the name of the variable saved
                                % (regressors.variable?)...
                            elseif strcmp(ext,'.txt')
                                regressors = load(job.subj.multi_reg{f},'-ascii');
                            end
                            for iReg = 1:size(regressors,2)
                                C = [C regressors(:,iReg)];
                                Cname = [Cname {['Other' int2str(iReg)]}];
                            end
                        catch
                            % Possibly no "Multiple Regressors" files for some
                            % sessions (currently only last ones can be
                            % omitted)
                        end
                    end
                    if NIRSconfoundsOn
                        wl = NIRS.Cf.dev.wl;
                        %HbO_like = [];
                        %for i=1:length(wl)
                        %    if wl(i) > 750 %in nanometer
                        %        %found a wavelength that is HbO-like
                        %        HbO_like = [HbO_like i];
                        %   end
                        %end
                        %HbO channels
                        %chHbO = NIRS.Cf.H.C.wl== HbO_like;
                        chHbO = 1:NC/2;
                        fullHbO = NIRS.Cf.H.C.gp(chHbO);
                        
                        [B_HbO IX] = sort(fullHbO); %sort in ascending order
                        %Impose minimum and maximum bounds
                        IXmax = IX(B_HbO <= MaxChDist);
                        IXmin = IX(B_HbO >= MinChDist);
                        IX2 = intersect(IXmin,IXmax);
                        %need to sort again
                        tHbO = fullHbO(IX2);
                        [dummy, IX3] = sort(tHbO);
                        %Keep up to NumChConfounds
                        try
                            IX3 = IX3(1:NumChConfounds);
                        end
                        %Find IX2 in full channel list
                        HbOIX = NIRS.Cf.H.C.id(:,chHbO);
                        HbOIX2 = HbOIX(:,IX2);
                        HbOIX3 = HbOIX2(1,IX3); %raw data channel
                        HbOIX4 = [];
                        for j1 = 1:length(HbOIX3) %channel in processed data rDtp(lst)
                            HbOIX4 = [HbOIX4 find(HbOIX3(j1)==NIRS.Cf.H.C.id(1,:))];
                        end
                        %get data for that session
                        d = fopen_NIR(rDtp{iSess,1},NC);
                        d_conf = d(HbOIX4,:);
                        
                        for j1=1:length(HbOIX4)
                            C = [C d_conf(j1,:)'];
                            Cname = [Cname {['C' int2str(j1)]}];
                        end
                        NIRSconfounds.NumChConfoundsActual = length(HbOIX4);
                        NIRSconfounds.Ch_removed = HbOIX4;
                        %Create and save a new data set excluding these
                        %channels for HbO and HbR
                        %generate list of kept channels
                        ch_keep = 1:NC;
                        kept_ch = ones(1,NC);
                        kept_ch(HbOIX4) = 0;
                        %if HbO_like == 1
                        kept_ch(HbOIX4+NC/2) = 0;
                        %else
                        %    kept_ch(HbOIX4-NC/2) = 0;
                        %end
                        ch_keep = ch_keep(logical(kept_ch));
                        [dir3, fil3, ext3] = fileparts(rDtp{iSess,1});
                        new_name = fullfile(spm_dir,[fil3 ext3]);
                        d_kept = d(ch_keep,:);
                        fwrite_NIR(new_name,d_kept);
                        %add outfile name to NIRS
                        if f == 1
                            lst = lst+1;
                            NIRS.Dt.fir.pp(lst).pre = 'Stat - Removed Channels';
                            NIRS.Dt.fir.pp(lst).job = job;
                        end
                        NIRS.Dt.fir.pp(lst).p{f,1} = new_name;
                        NIRS.Dt.fir.pp(lst).kept{f,1} = ch_keep; %kept channels
                    end
                catch
                    try
                        if job.GLM_include_cardiac
                            C = NIRS.Dt.fir.Sess(iSess).cR{1};
                            Cname = {'H'};
                        end
                    catch
                    end
                end
                
                SPM.Sess(f).C.C    = C;
                SPM.Sess(f).C.name = Cname;
            end
            
            %Number of datapoints for each session
            nscan = [];
            for f=1:nsess
                iSess = idx_sess(f);
                try
                    %only use of data for design specification. By storing
                    %size(d,2) in NIRS, we would avoid loading all the data!
                    d = fopen_NIR(rDtp{iSess,1},NC);
                catch
                    disp(['Aborting. Could not load data file for session ' int2str(iSess)]);
                    %return
                end
                nscan = [nscan size(d,2)];
            end
            clear d
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CODE from NIRS_SPM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            SPM.nscan = nscan;
            SPM.xY.RT = 1/fs;
            SPM.xBF.T = job.time_res;
            SPM.xBF.T0 = 1;
            SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
            if  job.units == 0
                SPM.xBF.UNITS = 'scans';
            elseif  job.units == 1
                SPM.xBF.UNITS = 'secs';
            end
            SPM.xBF.name = 'hrf';
            SPM.xBF = spm_get_bf(SPM.xBF);
            bf      = SPM.xBF.bf;
            % separate specifications for non-replicated sessions
            %--------------------------------------------------------------------------
            rep     = 0;
            
            Xx    = [];
            Xb    = [];
            Xname = {};
            Bname = {};
            
            for s = 1:length(SPM.nscan)
                % number of scans for this session
                %----------------------------------------------------
                k   = SPM.nscan(s);
                
                if (s == 1) || ~rep %always true
                    if need_onsets
                    % Get inputs, neuronal causes or stimulus functions U
                    %------------------------------------------------------------------
                    U = spm_get_ons(SPM,s);
                    
                    %Will need this to generate contrasts using SPM later
                    % Convolve stimulus functions with basis functions
                    %------------------------------------------------------------------
                    [X,Xn,Fc] = nirs_spm_Volterra(U,bf,1);
                    else
                        U = [];
                        X = [];
                        Xn{1} = [];
                        Fc = [];
                    end
                    % Resample regressors at acquisition times (32 bin offset)
                    %-------------------------------------------------
                    try
                        X = X((0:(k - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
                    end
                    
                    % and orthogonalise (within trial type)
                    %--------------------------------------
                    for i = 1:length(Fc)
                        X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
                    end
                    % %                     %Start with an empty design matrix; don't convolve with
                    % %                     %basis functions
                    % %                     X = []; %only used for regressing out confounds
                    % %                     Xn{1} = [];
                    % %                     Fc = [];
                    % get user specified regressors
                    %================================
                    try
                        C     = SPM.Sess(s).C.C;
                        Cname = SPM.Sess(s).C.name;
                    catch
                        % covariates - C
                        %----------------
                        str   = sprintf('Session %d',s);
                        spm_input('Other regressors',1,'d',str)
                        C     = [];
                        c = spm_input('user specified','+1','w1',0);
                        while size(C,2) < c
                            str = sprintf('regressor %i',size(C,2) + 1);
                            C  = [C spm_input(str,2,'e',[],[k Inf])];
                        end
                        % and their names - Cnames
                        %--------------------------
                        Cname = {};
                        for i = 1:size(C,2)
                            str      = sprintf('regressor %i',i);
                            Cname{i} = spm_input('name of','+0','s',str);
                        end
                    end
                    
                    % append mean-corrected regressors and names
                    %-------------------------------------------
                    reg_rows = size(C,1);
                    if (reg_rows > 0) && ~(reg_rows== k)
                        str1='Error in nirs_run_liom_intrasubject_average.m:';
                        str2=sprintf('Session %d has %d scans but regressors have %d entries', s,k,reg_rows);
                        str3='These numbers should match';
                        warndlg({str1; str2; str3});
                        return
                    end
                    X      = [X spm_detrend(C)];
                    Xn     = {Xn{:}   Cname{:}};
                    
                    % Confounds: Session effects
                    %===========================
                    B      = ones(k,1);
                    Bn     = {'constant'};
                    
                end
                
                % Session structure array
                %-----------------------------------------
                SPM.Sess(s).U      = U;
                SPM.Sess(s).C.C    = C;
                SPM.Sess(s).C.name = Cname;
                SPM.Sess(s).row    = size(Xx,1) + (1:k);
                SPM.Sess(s).col    = size(Xx,2) + (1:size(X,2));
                SPM.Sess(s).Fc     = Fc;
                
                % Append names
                %---------------------------------------------------------------
                for i = 1:length(Xn)
                    Xname{end + 1} = [sprintf('Sn(%i) ',s) Xn{i}];
                end
                for i = 1:length(Bn)
                    Bname{end + 1} = [sprintf('Sn(%i) ',s) Bn{i}];
                end
                
                % append into Xx and Xb
                %===============================================================
                Xx    = blkdiag(Xx,X);
                Xb    = blkdiag(Xb,B);
                
            end %- for s
            
            % finished
            %-----------------------------------------------------------------------
            SPM.xX.X      = [Xx Xb];
            SPM.xX.iH     = [];
            SPM.xX.iC     = 1:size(Xx,2);
            SPM.xX.iB     = (1:size(Xb,2)) + size(Xx,2);
            SPM.xX.iG     = [];
            SPM.xX.name   = {Xname{:} Bname{:}};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nscan = SPM.nscan;
            nsess = length(nscan);
            
            %Principal component removal
            SPM.xX.PCA = PCA;
            
            %Butterworth high pass filter
            SPM.xX.HPFbutter = HPFbutter;
            SPM.xX.hpf_butter_freq = hpf_butter_freq;
            SPM.xX.hpf_butter_order = hpf_butter_order;
            try
                SPM.GenerateHbT = job.GenerateHbT;
            catch
                SPM.GenerateHbT = 0;
            end
            try
                SPM.NIRSconfounds = NIRSconfounds;
            end
            %Add model specification
            try SPM.xX.opt = Opt; catch; end
            
            %%% updated for wavelet-MDL detrending 2009-03-19
            
            str = 'Detrending?';
            if isempty(strfind(HPF, 'wavelet')) == 0 % wavelet-MDL
                index_NT = find(HPF == ',');
                if isempty(index_NT) == 1
                    NT = 4;
                else
                    NT = str2num(HPF(index_NT+1:end));
                end
                SPM.xX.K.HParam.type = 'Wavelet-MDL';
                SPM.xX.K.HParam.M = NT;
            elseif isempty(strfind(HPF, 'DCT')) == 0 % DCT
                index_cutoff = find(HPF == ',');
                if isempty(index_cutoff) == 1
                    cutoff = 128;
                else
                    cutoff = str2num(HPF(index_cutoff+1:end));
                end
                SPM.xX.K.HParam.type = 'DCT';
                SPM.xX.K.HParam.M = cutoff;
            elseif isempty(strfind(HPF, 'none')) == 0 %no filter
                SPM.xX.K.HParam.type = 'none';
            end
            
            if isempty(strfind(LPF, 'hrf')) == 0 % hrf smoothing
                SPM.xX.K.LParam.type = 'hrf';
            elseif isempty(strfind(LPF, 'gaussian')) == 0 % Gaussian smoothing
                SPM.xX.K.LParam.FWHM = FWHM;
                SPM.xX.K.LParam.type = 'Gaussian';
            else
                SPM.xX.K.LParam.type = 'none';
            end
            
            
            %location of the HbO and HbR (combined) files - note that
            %for the purpose of Averaging, we do not care if a channel
            %is HbO or HbR, so we can loop over all the channels
            for f=1:nsess
                iSess = idx_sess(f);
                SPM.xY.P{f} = NIRS.Dt.fir.pp(lst).p{iSess};
            end
            SPM.job = job;
            spm_file = fullfile(spm_dir,'SPM.mat');
            save(spm_file,'SPM');
            NIRS.SPM{1} = spm_file;
            
            %NIRS is now modified - it includes a link to the GLM
            if NIRSconfoundsOn
                %update NIRS matrix
                NIRS.Cf.H.C.N = length(ch_keep);
                try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(ch_keep); end
                try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,ch_keep); end
                try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(ch_keep); end
                try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(ch_keep); end
                try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(ch_keep); end
                save(newNIRSlocation,'NIRS'); %This is essential
                
                %May need to generate a new topodata
                clear matlabbatch
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmat = {newNIRSlocation};
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.force_redo = 1;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.NIRSmatCopyChoice.NIRSmatOverwrite = struct([]);
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.template_mode = 0;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.anatT1 = {''};
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.segT1_4fit = {''};
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.anatT1_template = {'W:\spm8\templates\T1.nii'};
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.fid_in_subject_MNI = 0;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.nasion_wMNI = [0 84 -48];
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.AL_wMNI = [-83 -19 -38];
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.AR_wMNI = [83 -19 -38];
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.GenDataTopo = 1;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.render_choice.render_template = struct([]);
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.View6Projections = 0;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.Save6Projections = 1;
                matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg1.ForceReprocess = 0;
                spm_jobman('run',matlabbatch);
                [dir0 fil0] = fileparts(newNIRSlocation);
                NIRS.Dt.ana.rend = fullfile(dir0,'TopoData.mat');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Filtering and removing confounds; averaging
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [SPM NIRS] = nirs_liom_average(NIRS,SPM);
            SPM.FlagAvg = 1; %Flag to indicate that we are in averaging mode
            save(spm_file,'SPM');
            
            NIRS.flags.Avg_OK = 1;
            save(newNIRSlocation,'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not do averaging for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;