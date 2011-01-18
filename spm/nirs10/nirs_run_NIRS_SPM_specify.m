function out = nirs_run_NIRS_SPM_specify(job)
%Load NIRS.mat information
clear NIRS
load(job.NIRSmat{1,1});
%NIRS.derivs = job.derivs;
NIRS.onset_files = job.subj.input_onsets;
NIRS.multi_reg = job.subj.multi_reg;

%hb = 'hbo'; loop through all
%HPF_temp = job.nirs_hpf{1,1}; %'wavelet, 4';
%LPF_temp = job.nirs_lpf{1,1}; % 'hrf';

%HPF
try
    job.nirs_hpf.hpf_dct;
    HPF = ['DCT, ' int2str(job.nirs_hpf.hpf_dct.hpf_dct_cutoff)];
catch
    try job.nirs_hpf.hpf_wavelet;
        HPF = ['wavelet,' int2str(job.nirs_hpf.hpf_wavelet.hpf_wavelet_iter)];
        wavelet_depth = job.nirs_hpf.hpf_wavelet.hpf_wavelet_depth;
    catch
        try job.nirs_hpf.hpf_none;
            HPF = 'none';
        catch
            disp('Unrecognized high pass filter');
        end
    end
end

%LPF
try
    job.nirs_lpf.lpf_gauss;
    LPF = 'gaussian'; %, ' int2str(job.nirs_lpf.lpf_gauss.fwhm1)];
    FWHM = job.nirs_lpf.lpf_gauss.fwhm1;
catch
    try 
        job.nirs_lpf.lpf_hrf;
        LPF = ['hrf'];
    catch
        try
            job.nirs_lpf.lpf_none;
            LPF = 'none';
        catch
            disp('Unrecognized low pass filter');
        end
    end
end

method_cor = job.nirs_noise; %0;
flag_window = 1;
units = job.units; %1; %seconds
Volterra = job.volt;

%loop over HbO, HbR, HbT
for hbInd = 1:3
    switch hbInd
        case 1
            hb = 'HbO';
        case 2
            hb = 'HbR';
        case 3
            hb = 'HbT';
    end
    %loop over sessions
    for sessInd =1:size(NIRS.onset_files,1)
        clear SPM
        try
            load(NIRS.NIRS_SPM_Concfile{sessInd});
        catch
            disp(['Could not load data file for session' int2str(sessInd)]);
        end
        
        %PP added from spm_run_fmri_design for multi confound regressors
        C = [];
        Cname = {};
        try
            % User specified regressors
            %-------------------------------------------------------------           
            Cname = cell(1,numel(sess.regress));
            for q = 1:numel(sess.regress),
                Cname{q} = sess.regress(q).name;
                C         = [C, sess.regress(q).val(:)];
            end
        catch
        end

        % Augment the singly-specified regressors with the multiple regressors
        % specified in the regressors.mat file
        %------------------------------------------------------------
        try
            sess.multi_reg = {NIRS.multi_reg{sessInd}};
            if ~strcmp(sess.multi_reg,'')
                tmp=load(char(sess.multi_reg{:}));
                if isstruct(tmp) && isfield(tmp,'R')
                    R = tmp.R;
                elseif isnumeric(tmp)
                    % load from e.g. text file
                    R = tmp;
                else
                    warning('Can''t load user specified regressors in %s', ...
                        char(sess.multi_reg{:}));
                    R = [];
                end

                C=[C, R];
                nr=size(R,2);
                nq=length(Cname);
                for inr=1:nr,
                    Cname{inr+nq}=['R',int2str(inr)];
                end
            end
            SPM.Sess.C.C    = C;
            SPM.Sess.C.name = Cname;
        catch 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CODE from NIRS_SPM 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SPM.nscan = size(nirs_data.oxyData,1);
        SPM.xY.RT = 1/nirs_data.fs;
        SPM.xBF.T = job.time_res; %10; %???
        SPM.xBF.T0 = 1;
        SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
        if units == 0
                SPM.xBF.UNITS = 'scans';
            elseif units == 1
                SPM.xBF.UNITS = 'secs';
        end

        try
            %load onset file
            load(NIRS.onset_files{sessInd});
            for kk = 1:size(names, 2)
                SPM.Sess.U(kk).name = names(kk);
                SPM.Sess.U(kk).ons = onsets{kk};
                SPM.Sess.U(kk).dur = durations{kk};
                %Ignore parametric modulations - cf spm_run_fmri_design.m
                P.name = 'none';
                P.h    = 0;
                SPM.Sess.U(kk).P = P;
            end
                       
        catch
            %Could not load onset
            disp('Could not load onset file');
        end

        rep     = 0; %PP why?

        if all(job.derivs == [0 0])
                SPM.xBF.name = 'hrf';
            elseif all(job.derivs == [1 0])
                SPM.xBF.name = 'hrf (with time derivative)';
            elseif all(job.derivs == [1 1])
                SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
            else
                disp('Unrecognized hrf derivative choices.')
        end

        SPM.xBF = spm_get_bf(SPM.xBF); %PP removed call to nirs_spm_get_bf

        switch hb
            case 'HbO'
                bf = SPM.xBF.bf;
            case 'HbR'
                bf = SPM.xBF.bf ; %PP removed flipped sign as too confusing * (-1);
            case 'HbT'
                bf = SPM.xBF.bf;
        end

        V = Volterra;
        SPM.xBF.Volterra = V; % model interactions (Volterra) 

        %PP the rest is very close to spm_fmri_design.m, except for the
        %call of nirs_spm_get_ons_batch.m
        Xx    = [];
        Xb    = [];
        Xname = {};
        Bname = {};

        %PP Is the code treating each channel like a separate session?
        for s = 1:length(SPM.nscan)
            if (s == 1) || ~rep %PP always true
                k   = SPM.nscan(s);
                %U = nirs_spm_get_ons_batch(SPM, s, 1); %%PP not sure why need this special function
                U = spm_get_ons(SPM, s); %%PP not sure why need this special function

                [X,Xn,Fc] = spm_Volterra(U,bf,V);

                try
                    X = X([0:(k - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);
                end

                for i = 1:length(Fc)
                    X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
                end

                try
                    C     = SPM.Sess(s).C.C;
                    Cname = SPM.Sess(s).C.name;
                catch
        %             str   = sprintf('Session %d',s);
        %             spm_input('Other regressors',1,'d',str)
                    C     = [];
                    %             c     = spm_input('user specified','+1','w1',0);
                    c = 0;
                    while size(C,2) < c
                        str = sprintf('regressor %i',size(C,2) + 1);
                        C  = [C spm_input(str,2,'e',[],[k Inf])];
                    end

                    Cname = {};
                    for i = 1:size(C,2)
                        str      = sprintf('regressor %i',i);
                        Cname{i} = spm_input('name of','+0','s',str);
                    end
                end

                X      = [X spm_detrend(C)];
                Xn     = {Xn{:}   Cname{:}};

                B      = ones(k,1);
                Bn{1}  = sprintf('constant');

            end

            SPM.Sess(s).U      = U;
            SPM.Sess(s).C.C    = C;
            SPM.Sess(s).C.name = Cname;
            SPM.Sess(s).row    = size(Xx,1) + [1:k];
            SPM.Sess(s).col    = size(Xx,2) + [1:size(X,2)];
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
        SPM.xX.iC     = [1:size(Xx,2)];
        SPM.xX.iB     = [1:size(Xb,2)] + size(Xx,2);
        SPM.xX.iG     = [];
        SPM.xX.name   = {Xname{:} Bname{:}};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nscan = SPM.nscan;
        nsess = length(nscan);

        %%% updated for wavelet-MDL detrending 2009-03-19
        str = 'Detrending?';
        % try
        %     SPM.xX.K.HParam.type = cH.type;
        %     SPM.xX.K.HParam.M = cH.M;
        %     spm_input([str ' ' cH.type ',  ' num2str(cH.M)], '+1', 'd');
        % catch
        %     clear cH;
        %     cH = {'Wavelet-MDL', 'DCT'};
        %     cH = spm_input(str, '+1', 'b', cH);
        if isempty(strfind(HPF, 'wavelet')) == 0 % wavelet-MDL
            index_NT = find(HPF == ',');
            if isempty(index_NT) == 1
                NT = 4;
            else
                NT = str2num(HPF(index_NT+1:end));
            end
            SPM.xX.K.HParam.type = 'Wavelet-MDL';
            SPM.xX.K.HParam.M = NT;
            SPM.xX.K.HParam.wavelet_depth = wavelet_depth; %PP
        elseif isempty(strfind(HPF, 'DCT')) == 0 % DCT
            index_cutoff = find(HPF == ',');
            if isempty(index_cutoff) == 1
                cutoff = 128;
            else
                cutoff = str2num(HPF(index_cutoff+1:end));
            end
            SPM.xX.K.HParam.type = 'DCT';
            SPM.xX.K.HParam.M = cutoff;
        end

        if isempty(strfind(LPF, 'hrf')) == 0 % hrf smoothing
            SPM.xX.K.LParam.type = 'hrf';
        elseif isempty(strfind(LPF, 'gaussian')) == 0 % Gaussian smoothing
%             index_FWHM = find(LPF == ',');
%             if isempty(index_FWHM) == 1
%                 FWHM = 4;
%             else 
%                 FWHM = str2num(LPF(index_FWHM+1:end));
%             end
            SPM.xX.K.LParam.FWHM = FWHM;
            SPM.xX.K.LParam.type = 'Gaussian';
        else
            SPM.xX.K.LParam.type = 'none';
        end

        K = struct( 'HParam', SPM.xX.K.HParam,...
            'row', SPM.Sess.row,...
            'RT', SPM.xY.RT,...
            'LParam', SPM.xX.K.LParam);
        SPM.xX.K = spm_filter_HPF_LPF_WMDL(K);

        % related spm m-file : spm_fmri_spm_ui.m
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

        SPM.nirs.step = 'specification';
        SPM.nirs.fname = NIRS.NIRS_SPM_Concfile{sessInd};
        SPM.nirs.Hb = hb;
        SPM.nirs.level = 'individual';
        SPM_nirs = SPM;
        save(fullfile(job.dir{1,1},['SPM_indiv_' hb '_' int2str(sessInd) '.mat']), 'SPM_nirs');


    end % end for sessInd =1:size(NIRS.onset_files,1)
end %end for hbInd

save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
end