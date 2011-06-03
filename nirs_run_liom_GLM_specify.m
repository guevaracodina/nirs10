function out = nirs_run_liom_GLM_specify(job)

%TOO ANNOYING: was changed back to wls_bglm_specify -- just kept below as a
%reminder:
%Specify a GLM for each subject, each session and each chromophore
%Note: if a previously saved batch from a version of nirs10 prior to
%April 22 2011 fails to run, this maybe due to a renaming of the field
%wls_bglm_specify to liom_GLM_specify in job. If one needs to use such
%a batch, one could replace the wls_bglm_specify field by liom_GLM_specify
%using 
%job.('liom_GLM_specify')=job.('wls_bglm_specify')
%job = rmfield(job, 'wls_bglm_specify'); 
%ModulePosition=1;
%where job is the whole matlabbatch{ModulePosition}.spm.tools.nirs10.model_specify
%and saving the modified batch job

%MP=1;
%matlabbatch{MP}.spm.tools.nirs10.model_specify.('liom_GLM_specify')=matlabbatch{MP}.spm.tools.nirs10.model_specify.('wls_bglm_specify')
%matlabbatch{MP}.spm.tools.nirs10.model_specify = rmfield(matlabbatch{MP}.spm.tools.nirs10.model_specify, 'wls_bglm_specify'); 
%save('','matlabbatch')
%physiological confounds
try
    NIRSconfounds = job.NIRSchannelsConfound.NIRSconfounds;
    NIRSconfoundsOn = 1;
    NumChConfounds = NIRSconfounds.NumChConfounds;
    MinChDist = NIRSconfounds.MinChDist;
    MaxChDist = NIRSconfounds.MaxChDist;
catch
    NIRSconfoundsOn = 0;
end
%to generate display of design matrix
try
    flag_window = job.flag_window;
catch
    flag_window = 0;
end
%Option to skip generation of trRV. TrRV is required for statistics
try 
    generate_trRV = job.generate_trRV;
catch
    generate_trRV = 1;
end
%Option only used to search for a bug; correct option to use is 
%always filter_design_matrix = 0; -- This comment seems incorrect now --
%filter_design_matrix = 1 is required now to eliminate the bias when using
%the Butterworth high pass filter prior to the GLM
try 
    filter_design_matrix = job.filter_design_matrix;
catch
    filter_design_matrix = 0;
end

%Currently, only NIRS_SPM method works well
%Specify WLS, BGLM or NIRS_SPM parameters
meth0=job.wls_or_bglm;
try 
    meth0.NIRS_SPM;
    meth1 = 3;
catch
    try
        meth0.WLS;
        meth1 = 1;
    catch
        try
            meth0.BGLM;
            meth1 = 2;
        catch
            disp('Unrecognized method');
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
try 
    PCA = job.channel_pca;
catch
    PCA = 0; %no PCA
end

%HPF - Butterworth infinite impulse response filter
try
    hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPFbutter = 1;
catch
    HPFbutter = 0; %no high pass filter
    hpf_butter_freq = 0;
end

%LPF - Butterworth infinite impulse response filter
try    
    lpf_butter_freq = job.lpf_butter.lpf_butter_On.lpf_butter_freq; 
    LPFbutter = 1;
catch
    LPFbutter = 0; %no low pass filter
    lpf_butter_freq = 0;
end

%HPF - filter from NIRS_SPM - note that Butterworth HPF can be used with it
if meth1 ==3 
    try
        HPF = ['DCT, ' int2str(job.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_dct.hpf_dct_cutoff)];
    catch
        try job.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_wavelet;
            HPF = ['wavelet,' int2str(job.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_wavelet.hpf_wavelet_iter)];
            %wavelet_depth = job.nirs_hpf.hpf_wavelet.hpf_wavelet_depth;
        catch
            try job.wls_or_bglm.NIRS_SPM.nirs_hpf.hpf_none;
                HPF = 'none';
            catch
                disp('Unrecognized high pass filter');
            end
        end
    end
end

%LPF - filter from NIRS_SPM
if meth1 == 3
    try    
        FWHM = job.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_gauss.fwhm1;
        LPF = 'gaussian'; 
    catch
        try 
            job.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_hrf;
            LPF = 'hrf';
        catch
            try
                job.wls_or_bglm.NIRS_SPM.nirs_lpf.lpf_none;
                LPF = 'none';
            catch
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
        NIRS = [];
        load(job.NIRSmat{Idx,1});
           
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        nsess = size(rDtp,1);

        [dir1, fil1, ext1] = fileparts(rDtp{1,1});
        %create directory for stats for this subject
        spm_dir = fullfile(dir1,job.dir1);
        if ~exist(spm_dir,'dir'), mkdir(spm_dir); end

        %Find onsets
        try
            SPM.Sess = NIRS.Dt.fir.Sess;
        catch
            try
                %Ignore parametric modulations - cf spm_run_fmri_design.m
                P.name = 'none';
                P.h    = 0;
                for f=1:nsess
                    %load onset file
                    clear names onsets durations
                    load(job.subj(Idx,1).input_onsets{f});
                    for kk = 1:size(names, 2)
                        SPM.Sess(f).U(kk).name = names(kk);
                        SPM.Sess(f).U(kk).ons = onsets{kk};
                        SPM.Sess(f).U(kk).dur = durations{kk};
                        SPM.Sess(f).U(kk).P = P;
                    end
                end
            catch
                %Could not load onset
                disp('Could not find onsets');
            end
        end
        %Adding confound regressors 
        C = [];
        Cname = {};        
        for f=1:nsess
            try 
                if job.GLM_include_cardiac
                    %heart rate regressor
                    C = NIRS.Dt.fir.Sess(f).fR{1};
                    Cname = {'H'};
                end
                if job.GLM_include_Mayer
                    %Mayer wave regressor
                    C = [C NIRS.Dt.fir.Sess(f).mR{1}];
                    Cname = [Cname {'M'}];
                end
                if NIRSconfoundsOn
                    wl = NIRS.Cf.dev.wl;
                    HbO_like = [];
                    for i=1:length(wl)
                        if wl(i) > 750 %in nanometer
                            %found a wavelength that is HbO-like
                            HbO_like = [HbO_like i];
                        end
                    end
                    %HbO channels
                    chHbO = NIRS.Cf.H.C.wl== HbO_like;
                    fullHbO = NIRS.Cf.H.C.gp(chHbO); 
                    [B_HbO IX] = sort(fullHbO); %sort in ascending order
                    %Impose minimum and maximum bounds
                    IXmax = IX(B_HbO <= MaxChDist);
                    IXmin = IX(B_HbO >= MinChDist);
                    IX2 = intersect(IXmin,IXmax);
                    %need to sort again
                    tHbO = fullHbO(IX2);
                    [~, IX3] = sort(tHbO);
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
                    d = fopen_NIR(rDtp{f,1},NC);
                    d_conf = d(HbOIX4,:);
                    for j1=1:length(HbOIX3)
                        C = [C d_conf(j1,:)'];
                        Cname = [Cname {['C' int2str(j1)]}];
                    end 
                    NIRSconfounds.NumChConfoundsActual = length(HbOIX4);
                end
            catch
                try 
                    if job.GLM_include_cardiac
                        C = NIRS.Dt.fir.Sess(f).cR{1};
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
            try
                %only use of data for design specification. By storing 
                %size(d,2) in NIRS, we would avoid loading all the data!
                d = fopen_NIR(rDtp{f,1},NC);
            catch
                disp(['Aborting. Could not load data file for session ' int2str(f)]);
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
        
        % separate specifications for non-replicated sessions
        %--------------------------------------------------------------------------
        rep     = 0;
        %use field job.bases only to see if gamma was selected, otherwise
        %use old job.derivs field
        
        if strcmp(fieldnames(job.bases),'hrf')
            if all(job.derivs == [0 0])
                    SPM.xBF.name = 'hrf';
                elseif all(job.derivs == [1 0])
                    SPM.xBF.name = 'hrf (with time derivative)';
                elseif all(job.derivs == [1 1])
                    SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
                else
                    disp('Unrecognized hrf derivative choices.')
            end
        else
            nambase = fieldnames(job.bases);
            if ischar(nambase)
                nam=nambase;
            else
                nam=nambase{1};
            end
            switch nam,
                case 'fourier',
                    SPM.xBF.name = 'Fourier set';
                case 'fourier_han',
                    SPM.xBF.name = 'Fourier set (Hanning)';
                case 'gamma',
                    SPM.xBF.name = 'Gamma functions';
                case 'fir',
                    SPM.xBF.name = 'Finite Impulse Response';
                otherwise
                    error('Unrecognized hrf derivative choices.')
            end
            SPM.xBF.length = job.bases.(nam).length;
            SPM.xBF.order  = job.bases.(nam).order;
        end
        %the rest is very close to spm_fmri_design.m

        % get basis functions
        %--------------------------------------------------------------------------
        try
            bf      = SPM.xBF.bf;
        catch
            SPM.xBF = spm_get_bf(SPM.xBF);
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
                [X,Xn,Fc] = spm_Volterra(U,bf,V);

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
        %Butterworth low pass filter
        SPM.xX.LPFbutter = LPFbutter;
        SPM.xX.lpf_butter_freq = lpf_butter_freq;
        
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
                %SPM.xX.K.HParam.wavelet_depth = wavelet_depth; %PP
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
        %This is a longer calculation, that can potentially enlarge
        %considerably the SPM.mat structure, hence it is better left to the
        %estimate step
%         for s=1:nsess
%             K = struct( 'HParam', SPM.xX.K.HParam,...
%                 'row', SPM.Sess(s).row,...
%                 'RT', SPM.xY.RT,...
%                 'LParam', SPM.xX.K.LParam);
%             SPM.xX.K(s).K = spm_filter_HPF_LPF_WMDL(K); %???Indexing
%         end
% % %         %PP or instead:
% % %         allrow = [];
% % %         for s=1:nsess
% % %             allrow = [allrow SPM.Sess(s).row];
% % %         end
% % %         K = struct( 'HParam', SPM.xX.K.HParam,...
% % %                 'row', allrow,...
% % %                 'RT', SPM.xY.RT,...
% % %                 'LParam', SPM.xX.K.LParam);
% % %         SPM.xX.K = spm_filter_HPF_LPF_WMDL(K); %???Indexing
            
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

        %SPM.nirs.step = 'specification';
        %SPM.nirs.fname = NIRS.NIRS_SPM_Concfile{f};
        %SPM.nirs.Hb = hb;
        %SPM.nirs.level = 'individual';
        %SPM_nirs = SPM;
        
        %location of the HbO and HbR (combined) files - note that
        %for the purpose of GLM estimation, we do not care if a channel
        %is HbO or HbR, so we can loop over all the channels
        SPM.xY.P = NIRS.Dt.fir.pp(lst).p; 
        SPM.generate_trRV = generate_trRV;
        SPM.filter_design_matrix = filter_design_matrix;
        SPM.job = job;
        save(fullfile(spm_dir,'SPM.mat'),'SPM');
        %store path to SPM, after possible prior GLMs 
        try
            l1 = length(NIRS.SPM);
            NIRS.SPM{l1+1} = spm_dir;
        catch
            NIRS.SPM{1} = spm_dir;
        end
    catch
        disp(['Could not specify GLM for subject' int2str(Idx)]);
    end
    %NIRS is now modified - it includes a link to the GLM   
    newNIRSlocation = fullfile(spm_dir,'NIRS.mat');
    save(newNIRSlocation,'NIRS');
    job.NIRSmat{Idx,1} = newNIRSlocation;   
end
out.NIRSmat = job.NIRSmat;