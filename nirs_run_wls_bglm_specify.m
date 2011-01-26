function out = nirs_run_wls_bglm_specify(job)
%Specify a GLM for each subject, each session and each chromophore

%to generate display of design matrix
flag_window = 1;
NIRS.multi_reg = job.subj.multi_reg;
%group_sessions = job.group_sessions;
%NIRS.onset_files = job.subj.input_onsets;

%Begin by specifying filters, though we recommend filtering earlier,
%at the preprocessing stage, but maintained for convenience

%Specify WLS or BGLM parameters
meth0=job.wls_or_bglm;
try 
    meth0.NIRS_SPM;
    meth1 = 'NIRS_SPM';
catch
    try
        meth0.WLS;
        meth1 = 'WLS';
    catch
        try
            meth0.BGLM;
            meth1 = 'BGLM';
        catch
            disp('Unrecognized method');
        end
    end
end

switch meth1
    case 1
        Opt.meth = 'WLS';
        Opt.Design.L0=job.wls_or_bglm.WLS.WLS_L0; % maximum scale for signal decomposition=J-L0
        Opt.Design.J0=wls_or_bglm.WLS.WLS_J0;     % minimum scale to model physiology
        Opt.Design.threshold_drift=wls_or_bglm.WLS.WLS_threshold_drift; % Threshold for correlation analysis
    case 2
        Opt.meth = 'BGLM';
        Opt.Design.fmax=job.wls_or_bglm.BGLM.BGLM_fmax;   % maximum frequency for cosinusoidal drifts
        Opt.Design.degre=job.wls_or_bglm.BGLM.BGLM_degre; % maximum degree for polynomial drifts
        Opt.Design.threshold_drift=job.wls_or_bglm.BGLM.BGLM_threshold_drift; % Threshold for correlation analysis
    case 3
        Opt.meth = 'NIRS_SPM';       
    otherwise
end

%PCA
try 
    PCA = job.channel_pca;
catch
    PCA = 0; %no PCA
end

%HPF
try
    hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPFbutter = 1;
catch
    HPFbutter = 0; %no high pass filter
    hpf_butter_freq = 0;
end

%LPF
try    
    lpf_butter_freq = job.lpf_butter.lpf_butter_On.lpf_butter_freq; 
    LPFbutter = 1;
catch
    LPFbutter = 0; %no low pass filter
    lpf_butter_freq = 0;
end

%HPF
if strcmp(meth1,'NIRS_SPM')
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

%LPF
if strcmp(meth1,'NIRS_SPM')
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
        %copy NIRS.mat and datafile to spm_dir, for storage only
        if ~job.LiomDeleteLarge
            copyfile(job.NIRSmat{Idx,1},fullfile(spm_dir,'NIRS.mat'));
            for f=1:nsess
                copyfile(rDtp{f},fullfile(spm_dir,[fil1 ext1]));
            end
        end
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
                    load(job.subj{Idx,1}.input_onsets{f});
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
 
        %PP added from spm_run_fmri_design for multi confound regressors
%         C = [];
%         Cname = {};
%         try
%             % User specified regressors - currently not used
%             %-------------------------------------------------------------           
%             Cname = cell(1,numel(sess.regress));
%             for q = 1:numel(sess.regress),
%                 Cname{q} = sess.regress(q).name;
%                 C         = [C, sess.regress(q).val(:)];
%             end
%         catch
%         end

        % Augment the singly-specified regressors with the multiple regressors
        % specified in the regressors.mat file 
        %------------------------------------------------------------
        
        %Needs to be worked on, many loose ends
        %- generating movement regressors
        %- combining with heart rate regressor
        %loading from files vs available already in NIRS
        %interpolating on linspace to get the regressors
        %looping over all the files
% % %         try
% % %             
% % %             sess.multi_reg = {NIRS.multi_reg{f}};
% % %             if ~strcmp(sess.multi_reg,'')
% % %                 tmp=load(char(sess.multi_reg{:}));
% % %                 if isstruct(tmp) && isfield(tmp,'R')
% % %                     R = tmp.R;
% % %                 elseif isnumeric(tmp)
% % %                     % load from e.g. text file
% % %                     R = tmp;
% % %                 else
% % %                     warning('Can''t load user specified regressors in %s', ...
% % %                         char(sess.multi_reg{:}));
% % %                     R = [];
% % %                 end
% % % 
% % %                 C=[C, R];
% % %                 nr=size(R,2);
% % %                 nq=length(Cname);
% % %                 for inr=1:nr,
% % %                     Cname{inr+nq}=['R',int2str(inr)];
% % %                 end
% % %             end
% % %             SPM.Sess.C.C    = C;
% % %             SPM.Sess.C.name = Cname;
% % %         catch 
% % %             disp(['Could not load multi-regressors for subject ' int2str(Idx)]);
% % %         end

        %quick fix for now
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
        %For more precise specification of the design matrix,
        %ignored here
%                 if length(SPM.nscan) > 1 && ~any(diff(SPM.nscan)) && ~isfield(SPM,'Sess')
%                     rep = spm_input('are sessions replications','+1','yes|no',[1 0]);
%                 end

        if all(job.derivs == [0 0])
                SPM.xBF.name = 'hrf';
            elseif all(job.derivs == [1 0])
                SPM.xBF.name = 'hrf (with time derivative)';
            elseif all(job.derivs == [1 1])
                SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
            else
                disp('Unrecognized hrf derivative choices.')
        end

        %PP the rest is very close to spm_fmri_design.m, except for the
        %call of nirs_spm_get_ons_batch.m - now removed

        % get basis functions
        %--------------------------------------------------------------------------

        %SPM.xBF = spm_get_bf(SPM.xBF); %PP removed call to nirs_spm_get_bf

        try
            bf      = SPM.xBF.bf;
        catch
            SPM.xBF = spm_get_bf(SPM.xBF);
            bf      = SPM.xBF.bf;
        end

%                 switch hb
%                     case 'HbO'
%                         bf = SPM.xBF.bf;
%                     case 'HbR'
%                         bf = SPM.xBF.bf ; %PP removed flipped sign as too confusing * (-1);
%                     case 'HbT'
%                         bf = SPM.xBF.bf;
%                 end

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

        %PP Is the code treating each channel like a separate session?
        for s = 1:length(SPM.nscan)
            % number of scans for this session
            %----------------------------------------------------
            k   = SPM.nscan(s);

            if (s == 1) || ~rep %PP always true


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
        method_cor = job.wls_or_bglm.NIRS_SPM.nirs_noise;
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
    save(job.NIRSmat{Idx,1},'NIRS');    
end
out.NIRSmat = job.NIRSmat;