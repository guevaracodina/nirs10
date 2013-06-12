function out = nirs_run_liom_GLM_specify(job)
filter_vasomotion = 1;
filter_design_matrix = 1; % Michèle Aug 5, 2012 : completely eliminated this option 
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
%to generate display of design matrix
flag_window = job.flag_window;
%Option to skip generation of trRV. TrRV is required for statistics
generate_trRV = job.generate_trRV;

if isfield(job.vasomotion_choice,'vasomotion_on')
    vasomotion_on = 1;
    select_chromophore = job.vasomotion_choice.vasomotion_on.select_chromophore;
else
    vasomotion_on = 0;
end
%Currently, only NIRS_SPM method works well
%Specify WLS, BGLM or NIRS_SPM parameters
meth0=job.wls_or_bglm;
if isfield(meth0,'NIRS_SPM')
    meth0.NIRS_SPM;
    meth1 = 3;
else
    if isfield(meth0,'WLS')
        meth0.WLS;
        meth1 = 1;
    else
        if isfield(meth0,'BGLM')
            meth0.BGLM;
            meth1 = 2;
        end
    end
end

switch meth1
    case 1
        Opt.meth = 'WLS';
        Opt.Design.L0=job.wls_or_bglm.WLS.WLS_L0; % maximum scale for signal decomposition=J-L0
        Opt.Design.J0=job.wls_or_bglm.WLS.WLS_J0;     % minimum scale to model physiology
        Opt.Design.threshold_drift=job.wls_or_bglm.WLS.WLS_threshold_drift; % Threshold for correlation analysis
    case 2
        Opt.meth = 'BGLM';
        Opt.Design.fmax=job.wls_or_bglm.BGLM.BGLM_fmax;   % maximum frequency for cosinusoidal drifts
        Opt.Design.degre=job.wls_or_bglm.BGLM.BGLM_degre; % maximum degree for polynomial drifts
        Opt.Design.threshold_drift=job.wls_or_bglm.BGLM.BGLM_threshold_drift; % Threshold for correlation analysis
    case 3
        Opt.meth = 'NIRS_SPM';
    otherwise
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

%HPF - filter from NIRS_SPM - note that Butterworth HPF can be used with it
if meth1 ==3
    if isfield(meth0.NIRS_SPM.nirs_hpf,'hpf_dct')
        HPF = ['DCT, ' int2str(meth0.NIRS_SPM.nirs_hpf.hpf_dct.hpf_dct_cutoff)];
    else
        if isfield(meth0.NIRS_SPM.nirs_hpf,'hpf_wavelet')
            HPF = ['wavelet,' int2str(meth0.NIRS_SPM.nirs_hpf.hpf_wavelet.hpf_wavelet_iter)];
            %wavelet_depth = job.nirs_hpf.hpf_wavelet.hpf_wavelet_depth;
        else
            if isfield(meth0.NIRS_SPM.nirs_hpf,'hpf_none')
                HPF = 'none';
            else
                disp('Unrecognized high pass filter');
            end
        end
    end
end

%LPF - filter from NIRS_SPM
if meth1 == 3
    if isfield(meth0.NIRS_SPM.nirs_lpf,'lpf_gauss')
        FWHM = meth0.NIRS_SPM.nirs_lpf.lpf_gauss.fwhm1;
        LPF = 'gaussian';
    else
        if isfield(meth0.NIRS_SPM.nirs_lpf,'lpf_hrf')
            LPF = 'hrf';
        else
            if isfield(meth0.NIRS_SPM.nirs_lpf,'lpf_none')
                LPF = 'none';
            else
                disp('Unrecognized low pass filter');
            end
        end
    end
end

%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        %Objective is to fill SPM structure
        SPM = [];
        %always store SPM analysis in some directory
        if ~isfield(job.NIRSmatCopyChoice,'NIRSmatCopy')
            job.NIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir = 'Stat';
        end            
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS,'flags') || ~isfield(NIRS.flags,'GLMspec_OK') || job.force_redo)
            [spm_dir dummy] = fileparts(newNIRSlocation);
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
            
            %Adding confound regressors
            
            for f=1:nsess
                iSess = idx_sess(f);
                C = [];
                Cname = {};
                try
                    if job.GLM_include_cardiac
                        %heart rate regressor
                        C = NIRS.Dt.fir.Sess(iSess).fR{1};
                        %remove NaNs
                        Cnan = find(isnan(NIRS.Dt.fir.Sess(iSess).fR{1}));
                        Ctemp = C;
                        Ctemp(Cnan) = [];
                        Cmean = mean(Ctemp);
                        C(Cnan) = Cmean;
                        C = C - repmat(Cmean,[length(C),1]);
                        Cname = {'H'};
                        if any(isnan(C))
                             C = [];
                             Cname = {};
                        end
                    end
                    if job.GLM_include_Mayer
                        %Mayer wave regressor
                        C = [C NIRS.Dt.fir.Sess(iSess).mR{1}];
                        Cname = [Cname {'M'}];
                    end
                    %wl = NIRS.Cf.dev.wl;
%                     HbO_like = [];
%                     for i=1:length(wl)
%                         if wl(i) > 750 %in nanometer
%                             %found a wavelength that is HbO-like
%                             HbO_like = [HbO_like i];
%                         end
%                     end
                     
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
%                         HbO_like = [];
%                         for i=1:length(wl)
%                             if wl(i) > 750 %in nanometer
%                                 %found a wavelength that is HbO-like
%                                 HbO_like = [HbO_like i];
%                             end
%                         end
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
                 
            
            % separate specifications for non-replicated sessions
            %--------------------------------------------------------------------------
            rep     = 0;
            %use field job.bases only to see if gamma was selected, otherwise
            %use old job.derivs field
            
            SPM.xBF = nirs_specify_bases(job);
            SPM.xBF.T = job.time_res;
            SPM.xBF.T0 = 1;
            SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
       
            %the rest is very close to spm_fmri_design.m
            
            % get basis functions
            %--------------------------------------------------------------------------
            try
                bf      = SPM.xBF.bf;
            catch
                SPM.xBF = nirs_get_bf(SPM.xBF);
                if size(SPM.xBF.bf,1) == 1 || size(SPM.xBF.bf,2) == 1
                    SPM.xBF.bf = SPM.xBF.bf/sum(SPM.xBF.bf); %normalize
                end
                bf      = SPM.xBF.bf;
            end
            
            V = job.volt;
            SPM.xBF.Volterra = V; % model interactions (Volterra)
            
            % 1st or 2nd order Volterra expansion?
            %--------------------------------------------------------------------------
            try
                V   = SPM.xBF.Volterra;
            catch
                V   = spm_input('model interactions (Volterra)','+1','y/n',[2 1]);
                SPM.xBF.Volterra  = V;
            end
            
            Xx    = [];
            Xb    = [];
            Xname = {};
            Xname_short = {};
            Bname = {};
            
            for s = 1:length(SPM.nscan)
                % number of scans for this session
                %----------------------------------------------------
                k   = SPM.nscan(s);
                
                if (s == 1) || ~rep %always true
                    % create convolved stimulus functions or inputs
                    %==================================================================
                    
                    % Get inputs, neuronal causes or stimulus functions U
                    %------------------------------------------------------------------
                    U = spm_get_ons(SPM,s);
                    
                    % Convolve stimulus functions with basis functions
                    %------------------------------------------------------------------
                    [X,Xn,Fc] = nirs_spm_Volterra(U,bf,V);
                    
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
                    %To orthogonalize 2nd Volterra:
                    %tX = spm_orth(X(:,1:2));
                    %X = [tX X(:,3:end)];
                    
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
                        str1='Error in nirs_run_NIRS_SPMspecify_batch.m:';
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
                    Xname_short{end+1} = Xn{i};
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
            SPM.xX.name_short   = {Xname_short{:} Bname{:}};
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
            
            
            
            %**********************************************************
            %Temporarily added by Ke Peng
            %**********************************************************

            if meth1 == 1
                SPM.xX.K.HParam.type = 'none';
                SPM.xX.K.LParam.type = 'none'; %Need to modify if using HPFs or LPFs
            end
            
            %**********************************************************
            
            
            %%% updated for wavelet-MDL detrending 2009-03-19
            if meth1 == 3
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
            end
         
            % related spm m-file : spm_fmri_spm_ui.m
            if meth1 == 3
                method_cor = job.wls_or_bglm.NIRS_SPM.nirs_noise;
            else
                method_cor = 0;
            end
            if method_cor == 0
                cVi = 'none';
            elseif method_cor == 1
                cVi = 'AR(1)';
            end
            
            if ~ischar(cVi)	% AR coeficient[s] specified
                SPM.xVi.Vi = spm_Ce(nscan,cVi(1:3));
                cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];
                
            else
                switch lower(cVi)
                    case 'none'		%  xVi.V is i.i.d
                        %---------------------------------------------------------------
                        SPM.xVi.V  = speye(sum(nscan));
                        cVi        = 'i.i.d';
                    otherwise		% otherwise assume AR(0.2) in xVi.Vi
                        %---------------------------------------------------------------
                        SPM.xVi.Vi = spm_Ce(nscan,0.2);
                        cVi        = 'AR(0.2)';
                end
            end
            SPM.xVi.form = cVi;
            
            SPM.xsDes = struct('Basis_functions', SPM.xBF.name, 'Sampling_period_sec', num2str(SPM.xY.RT), 'Total_number_of_samples', num2str(SPM.nscan));
            if flag_window == 1
                spm_DesRep('DesMtx',SPM.xX,[],SPM.xsDes)
            end
            
            %location of the HbO and HbR (combined) files - note that
            %for the purpose of GLM estimation, we do not care if a channel
            %is HbO or HbR, so we can loop over all the channels
            for f=1:nsess
                iSess = idx_sess(f);
                SPM.xY.P{f} = NIRS.Dt.fir.pp(lst).p{iSess};
            end
            if isfield(job,'TrRVRVexact')
                SPM.TrRVRVexact = job.TrRVRVexact;
            else
                SPM.TrRVRVexact = 0; %approximate
            end
            SPM.generate_trRV = generate_trRV;
            SPM.filter_design_matrix = filter_design_matrix;
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
                %Important: need to save NIRS before running batch coreg,
                %so that the channels are updated! 
                save(newNIRSlocation,'NIRS');
                nirs_batch_coreg(NIRS,newNIRSlocation);
                rend_file = fullfile(spm_dir,'TopoData.mat');
                NIRS.Dt.ana.rend = rend_file;                
            end
            NIRS.flags.GLMspec_OK = 1;
            save(newNIRSlocation,'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not specify GLM for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
    
end
out.NIRSmat = job.NIRSmat;