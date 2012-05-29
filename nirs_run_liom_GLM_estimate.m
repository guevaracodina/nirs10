function out = nirs_run_liom_GLM_estimate(job)
%NIRS_SPM GLM estimation - first level
which_GLM = job.NIRS_SPM_which_GLM;
%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'GLMest_OK') || job.force_redo)
            NC = NIRS.Cf.H.C.N;
            fs = NIRS.Cf.dev.fs;
            try
                switch which_GLM
                    case 1 %first
                        fGLM = NIRS.SPM(1);
                    case 2
                        %need to loop
                        fGLM = NIRS.SPM;
                    case 3
                        fGLM = NIRS.SPM(end);
                end
            catch
                try
                    fGLM = NIRS.SPM(1);
                catch
                    disp(['Could not find a GLM for subject ' int2str(Idx)]);
                end
            end
            
            %loop over GLMs to estimate for a given subject and set of sessions
            %(usually only one such GLM)
            for g=1:length(fGLM)
                SPM = [];
                iSPM = 0; %count total number of subsessions
                
                load(fGLM{g});
                [spm_dir dummy] = fileparts(fGLM{g});
                
                try
                    SPM.xY.P;
                catch
                    disp(['Could not find data file for subject ' int2str(Idx)]);
                end
                
                %Prepare SPM matrix -- if fields are present / removing field
                try
                    if isfield(SPM.xVi, 'V') == 1
                        SPM = rmfield(SPM, 'xVi');
                        precolor = 1;
                    elseif isfield(SPM_nirs.xVi, 'V') == 0
                        precolor = 0;
                    end
                catch
                    disp('No V field present - GLM already estimated. Estimating it again!');
                    precolor = 1;
                end
                
                
                %loop over sessions
                nsess = length(SPM.xY.P);
                for s=1:nsess
                    d = fopen_NIR(SPM.xY.P{s},NC);
                    %loop over subsessions, defined as period in between
                    %movement intervals -- to do later
                    
                    try
                        lst = length(NIRS.Dt.fir.pp);
                        %might not be referencing to the right data file
                        bpi = NIRS.Dt.fir.pp(lst).bpi{s,1}; %bad point indices
                        bpd = NIRS.Dt.fir.pp(lst).bpd{s,1}; %bad point durations
                        if isempty(bpi)
                            markers_available = 0;
                        else
                            markers_available = 1;
                            si = NIRS.Dt.fir.pp(lst).si{s,1};
                            ei = NIRS.Dt.fir.pp(lst).ei{s,1};
                        end
                    catch
                        markers_available = 0;
                    end
                    
                    if markers_available
                        nSubSess = length(si);
                    else
                        nSubSess = 1;
                    end
                    %nSubSess is the number of subsessions for session f
                    for iSubSess = 1:nSubSess
                        %Extract data
                        if markers_available
                            %Need to transpose
                            Y = d(:,si(iSubSess):ei(iSubSess))';
                        else
                            Y = d';
                        end
                        sigma_unf = std(Y,0,1);
                        % PCA
                        if SPM.xX.PCA
                            nComponents = 1;
                            %process HbO, HbR separately
                            which_channels = 1:(NC/2); %all HbOchannels
                            Y1 = makePca(Y,which_channels,nComponents);
                            Y2 = makePca(Y,which_channels+NC/2,nComponents);
                            Y = [Y1 Y2];
                        end
                        
                        %LPF
                        %                     if SPM.xX.LPFbutter
                        %                         cutoff=SPM.xX.lpf_butter_freq; %0.666; %Hz, or 1.5s
                        %                         FilterOrder=5; %Is this too weak?
                        %                         Y = ButterLPF(fs,cutoff,FilterOrder,Y);
                        %                     end
                        
                        %HPF
                        switch SPM.xX.HPFbutter
                            case 1
                                cutoff=SPM.xX.hpf_butter_freq; %in Hz,
                                FilterOrder=SPM.xX.hpf_butter_order; %Is this too weak?
                                Y = ButterHPF(fs,cutoff,FilterOrder,Y);
                                %need to filter the design matrix too,
                                %otherwise, the estimates will be significantly
                                %biased
                            case 2
                                nS = size(Y,1);
                                mX = linspace(0,round(nS/fs),nS);
                                mX = [mX' ones(nS,1)];
                                pmX = pinv(mX);    
                                Y = Y - mX * (pmX * Y);                                
                        end
                        
                        
                        %Last step before the GLM, adding HbO to HbR to get
                        %HbT - perhaps this should be done after the
                        %wavelets?
                        if SPM.GenerateHbT
                            tmp_ch = 1:(NC/2); %all HbOchannels
                            Y1 = Y(:,tmp_ch);
                            Y2 = Y(:,tmp_ch+NC/2);
                            Y = [Y1 Y2 Y1+Y2];
                        end
                        %carefully extract SPM info
                        tSPM = [];
                        tSPM.Sess = SPM.Sess(s);
                        tSPM.xX = SPM.xX;
                        try
                            %This only worked for constant number of regressors
                            %across sessions
                            %                         %find elements of X for session s
                            %                         nbeta = size(SPM.xX.X,2);
                            %                         nbetaS = (nbeta-nsess)/nsess;
                            %                         %last entry is the constant regressor
                            %                             beta = [(s-1)*nbetaS+1:s*nbetaS nbeta-nsess+s];
                            
                            beta = [SPM.Sess(s).col SPM.Sess(nsess).col(end)+s];
                            
                            if markers_available
                                svec = SPM.Sess(s).row(si(iSubSess):ei(iSubSess));
                                
                            else
                                svec = SPM.Sess(s).row;
                            end
                            K = struct( 'HParam', SPM.xX.K.HParam,...
                                'row', svec ,...
                                'RT', SPM.xY.RT,...
                                'LParam', SPM.xX.K.LParam);
                            tSPM.xX.K = spm_filter_HPF_LPF_WMDL(K);
                            
                            tSPM.xX.X = SPM.xX.X(svec,beta);
                            tSPM.xX.K.row = 1:length(svec);
                        catch
                            %Did not correctly extract this session's design
                            %matrix?
                        end
                        %Could lead to problems - got huge t-stats -
                        %perhaps this introduces correlations between data
                        %and protocol - perhaps because of boundary
                        %interpolation errors? - or more likely: because it
                        %introduces correlations into the design matrix
                        %that falsify the calculation of the number of
                        %degrees of freedom
                        try
                            if SPM.filter_design_matrix
                                switch SPM.xX.HPFbutter
                                    case 1
                                    %filter the design matrix
                                    cutoff=SPM.xX.hpf_butter_freq; %Hz, or 100s
                                    FilterOrder=SPM.xX.hpf_butter_order; %Is this too weak?
                                    %exclude the constant
                                    tX=ButterHPF(fs,cutoff,FilterOrder,tSPM.xX.X(:,1:end-1));
                                    %add back the constant
                                    tSPM.xX.X = [tX tSPM.xX.X(:,end)];
                                    case 2
                                        %do nothing -- no linear trend in
                                        %design matrix
                                    case 3
                                        %add a linear trend to the design
                                        %matrix
                                        nS = size(tSPM.xX.X,1);
                                        mX = linspace(0,round(nS/fs),nS)/(nS/fs);
                                        tSPM.xX.X = [tSPM.xX.X(:,1:end-1) mX' tSPM.xX.X(:,end)];
                                    case 4 %add cosines like in SPM
                                        % make high pass filter
                                        %------------------------------------------------------------------
                                        k       = size(tSPM.xX.X,1);
                                        n       = fix(2*(k/fs)*SPM.xX.hpf_butter_freq + 1);
                                        X0      = spm_dctmtx(k,n);
                                        X0 = X0(:,2:end)./repmat(max(X0(:,2:end),[],1),[k 1]);
                                        tSPM.xX.X = [tSPM.xX.X(:,1:end-1) X0 tSPM.xX.X(:,end)];                                        
                                end
                                %                             %LPF
                                %                             if SPM.xX.LPFbutter
                                %                                 cutoff=SPM.xX.lpf_butter_freq; %0.666; %Hz, or 1.5s
                                %                                 FilterOrder=5; %Is this too weak?
                                %                                 tX=ButterLPF(fs,cutoff,FilterOrder,Wn,tSPM.xX.X(:,1:end-1));            % buterworth filter
                                %                                 %add back the constant
                                %                                 tSPM.xX.X = [tX tSPM.xX.X(:,end)];
                                %                             end
                            end
                        end
                        try tSPM.generate_trRV = SPM.generate_trRV; end
                        switch SPM.xX.opt.meth
                            case 'BGLM'
                                nScan = size(Y,1);
                                time = (1:nScan)/fs;
                                fmax = tSPM.xX.opt.Design.fmax;
                                degre = tSPM.xX.opt.Design.degre;
                                threshold_corr = tSPM.xX.opt.Design.threshold_drift;
                                [D correlated_infos] = makeDCT(time',fmax,degre,...
                                    threshold_corr,tSPM.xX.X);
                                tSPM.Sess.C.C = [tSPM.Sess.C.C D];
                                
                                [Betas,Thetas,Sigma2,Modelisation] = ...
                                    glm(time,Y,tSPM.xX.X(:,1),[tSPM.xX.X(:,2:end-1) tSPM.Sess.C.C]);
                                %glm(time,Y,tSPM.xX.X(:,1:end-1),tSPM.Sess.C.C); %exclude constant in X
                            case 'WLS'
                                %remove constant regressor - assume it is the last entry
                                %if size(tSPM.xX.X,2) > 1
                                 %   tmpX = tSPM.xX.X(:,1:end-1);
                                %else
                                %    tmpX = tSPM.xX.X;
                                %end

                                %**********************************************************
                                %Temporarily added by Ke Peng
                                %**********************************************************                                
                          
                                %remove constant regressor - assume it is the last entry
                                if size(tSPM.xX.X,2) > 1
                                    tmpX = tSPM.xX.X(:,1:end-1);
                                else
                                    tmpX = tSPM.xX.X;
                                end
                                
                                nScan = size(tSPM.xX.X,1);
                                
                                switch tSPM.xX.K.HParam.type
                                    case 'Wavelet-MDL'
                                        tmp_K = tSPM.xX.K;
                                        tmp_K.HParam.type = '';
                                        tSPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(tmp_K, tSPM.xX.X)); % KX
                                        tSPM.xX.K.X = tSPM.xX.X;
                                        clear tmp_K;
                                    case 'DCT'
                                        tSPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(tSPM.xX.K, tSPM.xX.X)); % KX
                                    case 'none'
                                        tSPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(tSPM.xX.K, tSPM.xX.X)); % KX ?
                                end

                                tSPM.xX.xKXs.X = full(tSPM.xX.xKXs.X);
                                tSPM.xX.pKX = spm_sp('x-', tSPM.xX.xKXs);
                                
                                
                                switch tSPM.xX.K.LParam.type
                                case {'hrf', 'Gaussian'}
                                    S = tSPM.xX.K.KL;
                                case 'none'
                                    S = speye(nScan);
                                end

                                switch tSPM.xX.K.HParam.type
                                    case 'DCT'
                                        S = S - tSPM.xX.K.X0 * (tSPM.xX.K.X0' * S);
                                        %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
                                end
                                if ~isfield(tSPM,'TrRVRVexact')
                                    tSPM.TrRVRVexact = 0; %approximation
                                end

                                if tSPM.generate_trRV
                                    [trRV trRVRV ] = approx_trRV(tSPM.xX.xKXs.X,tSPM.xX.pKX,S,[],tSPM.TrRVRVexact);
                                else
                                    trRV = 0;
                                    trRVRV = 0;
                                end
                                tSPM.xX.trRV = trRV; % <R'*y'*y*R>
                                tSPM.xX.trRVRV = trRVRV; %- Satterthwaite
                                try 
                                    tSPM.xX.erdf = trRV^2/trRVRV; 
                                catch exception
                                    disp(exception.identifier);
                                    disp(exception.stack(1));
                                    disp('Problem calculating degrees of freedom');    
                                end
                            
                                tSPM.xX.Bcov = (tSPM.xX.pKX * S);
                                tSPM.xX.Bcov = tSPM.xX.Bcov * tSPM.xX.Bcov';
                                
                                [Betas,spectralExponents,Modelisation,Design,Ress,ResSS] = ...
                                    wls(fs,Y,tmpX,tSPM.xX.opt.Design,tSPM.xX.K);
                                
                                %[Betas,spectralExponents,Modelisation,Design] = ...
                                %    wls(fs,Y,tmpX,tSPM.xX.opt.Design);
                                
                                tSPM.KY = Y;
                                tSPM.xX.beta = Betas;
                                res = Ress;
                                tSPM.xX.ResSS = sum(res.^2);
                                tSPM.xX.ResSSch = ResSS; % Residual sum of squares
                                
                                
                                tSPM.xX.xKXs.X = tSPM.xX.xKXs.X(:,1:(end-1));

                                %**********************************************************                                
                                                                 
                            case 'NIRS_SPM'
                                if precolor
                                    [tSPM res] = precoloring_batch(tSPM,Y);
                                else
                                    %not done yet
                                    tSPM = prewhitening(tSPM, Y);
                                end
                        end
                        %fill in beta and var into tSPM.xX
                        switch SPM.xX.opt.meth
                            case {'BGLM', 'WLS'}
                                %fill tSPM
                                %tSPM.xX.beta is regressors times channels
                                tSPM.xX.beta = Betas.betas';
                                tSPM.xX.Bvar = Betas.betasVariances';
                                %tstat
                                try
                                    tSPM.xX.t=tSPM.xX.beta./sqrt(tSPM.xX.Bvar);
                                end
                            case 'NIRS_SPM'
                                for r=1:size(tSPM.xX.beta,1)
                                    try
                                        tSPM.xX.t(r,:) = tSPM.xX.beta(r,:)./...
                                            sqrt(tSPM.xX.ResSS(:)'*tSPM.xX.Bcov(r,r)/tSPM.xX.trRV);
                                    catch exception
                                        disp(exception.identifier)
                                        disp(exception.stack(1));
                                        disp('Problem calculating t stat');
                                    end
                                end
                                %b2 =SPM.xXn{1}.beta(1,:)./sqrt(SPM.xXn{1}.ResSS/SPM.xXn{1}.trRV);
                        end
                        %A posteriori estimation of the amplitude of the response
                        try
                            %filtered data less estimated beta times filtered design
                            %matrix = adjusted data
                            if length(beta) == 6 %HARD CODED - Careful! to catch case with 2 types of onsets for patient 1
                                V2r = 3; %position of 2nd Volterra regressor
                            else
                                V2r = 2;
                            end
                            

                            
                            %Need filtered Y (KY) and filtered X ()
                            %sigma = std(tSPM.KY- (tSPM.xX.X(:,1)*tSPM.xX.beta(1,:) + tSPM.xX.X(:,V2r)*tSPM.xX.beta(V2r,:)),0,1);
                            %sigma = std(tSPM.KY- (tSPM.xX.xKXs.X(:,1)*tSPM.xX.beta(1,:) + tSPM.xX.xKXs.X(:,V2r)*tSPM.xX.beta(V2r,:)),0,1);
                            

                            sigma = std(tSPM.KY- (tSPM.xX.xKXs.X*tSPM.xX.beta),0,1);
                            
                            %**********************************************************
                            %Temporarily added by Ke Peng
                            %**********************************************************

                            if strcmp(SPM.xX.opt.meth, 'WLS')
                                V2r = size(tSPM.xX.xKXs.X,2);
                                tSPM.xX.xKXs.X = [tSPM.xX.xKXs.X ones(length(tSPM.xX.xKXs.X), 1)];% to fill the blanc left for constant regressor
                                tSPM.xX.beta = [tSPM.xX.beta' zeros(length(tSPM.xX.beta), 1)]';
                            end % Very bad coding! Need to verify whether WLS removed the constant regressor. What if there are multiple kinds of spks?
                            %**********************************************************
                            
                            tSPM.xX.sigma = sigma;
                            tSPM.xX.sigma_unf = sigma_unf;
                            tSPM.xX.sigma_filt_unadj = std(tSPM.KY,0,1);
                            if strcmp(SPM.xBF.name,'Gamma Functions')
                                try
                                    if NIRS.Dt.fir.calculate_bf_norm
                                        tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma)/4.4631;
                                    else
                                        tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma);
                                    end
                                catch
                                    tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma);
                                end
                            else %For canonical HRF
                                try
                                    if NIRS.Dt.fir.calculate_bf_norm
                                        tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma)/4.7506;
                                    else
                                        tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma);
                                    end
                                catch
                                    tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma);
                                end %norm_bf = 4.7506 = 1/max(X(:,1)) for a protocol with only one spike
                            end
                            
                            tSPM.xX.b  = tSPM.xX.beta(V2r,:) ./ tSPM.xX.beta(1,:);
                            tSPM.xX.b2  = tSPM.xX.beta(V2r,:) ./ sigma;
                            %tSPM.xX.a2 and b2 give the sigma-normalized
                            %amplitudes, for each channel
                            %Careful, HARD-CODED, this might not be correct
                            %order of HbO and HbR channels
                            CHbO = 1:(NC/2);
                            CHbR = (1+(NC/2)):NC;
                            [tSPM.xX.S.OtaVmax, tSPM.xX.S.OtaCmax] = max(tSPM.xX.t(1,CHbO));
                            [tSPM.xX.S.OtaVmin, tSPM.xX.S.OtaCmin] = min(tSPM.xX.t(1,CHbO));
                            [tSPM.xX.S.OtbVmax, tSPM.xX.S.OtbCmax] = max(tSPM.xX.t(V2r,CHbO));
                            [tSPM.xX.S.OtbVmin, tSPM.xX.S.OtbCmin] = min(tSPM.xX.t(V2r,CHbO));
                            [tSPM.xX.S.RtaVmax, tSPM.xX.S.RtaCmax] = max(tSPM.xX.t(1,CHbR));
                            [tSPM.xX.S.RtaVmin, tSPM.xX.S.RtaCmin] = min(tSPM.xX.t(1,CHbR));
                            [tSPM.xX.S.RtbVmax, tSPM.xX.S.RtbCmax] = max(tSPM.xX.t(V2r,CHbR));
                            [tSPM.xX.S.RtbVmin, tSPM.xX.S.RtbCmin] = min(tSPM.xX.t(V2r,CHbR));
                            [tSPM.xX.S.OaVmax, tSPM.xX.S.OaCmax] = max(tSPM.xX.beta(1,CHbO)./ sigma(CHbO));
                            [tSPM.xX.S.OaVmin, tSPM.xX.S.OaCmin] = min(tSPM.xX.beta(1,CHbO)./ sigma(CHbO));
                            [tSPM.xX.S.ObVmax, tSPM.xX.S.ObCmax] = max(tSPM.xX.beta(V2r,CHbO)./ sigma(CHbO));
                            [tSPM.xX.S.ObVmin, tSPM.xX.S.ObCmin] = min(tSPM.xX.beta(V2r,CHbO)./ sigma(CHbO));
                            [tSPM.xX.S.RaVmax, tSPM.xX.S.RaCmax] = max(tSPM.xX.beta(1,CHbR)./ sigma(CHbR));
                            [tSPM.xX.S.RaVmin, tSPM.xX.S.RaCmin] = min(tSPM.xX.beta(1,CHbR)./ sigma(CHbR));
                            [tSPM.xX.S.RbVmax, tSPM.xX.S.RbCmax] = max(tSPM.xX.beta(V2r,CHbR)./ sigma(CHbR));
                            [tSPM.xX.S.RbVmin, tSPM.xX.S.RbCmin] = min(tSPM.xX.beta(V2r,CHbR)./ sigma(CHbR));

                            [tSPM.xX.S.OtaVmax, tSPM.xX.S.OtaCmax] = max(tSPM.xX.t(1,CHbO));
                            [tSPM.xX.S.OtaVmin, tSPM.xX.S.OtaCmin] = min(tSPM.xX.t(1,CHbO));
                           % [tSPM.xX.S.OtbVmax, tSPM.xX.S.OtbCmax] = max(tSPM.xX.t(V2r,CHbO));
                           % [tSPM.xX.S.OtbVmin, tSPM.xX.S.OtbCmin] = min(tSPM.xX.t(V2r,CHbO));
                            [tSPM.xX.S.RtaVmax, tSPM.xX.S.RtaCmax] = max(tSPM.xX.t(1,CHbR));
                            [tSPM.xX.S.RtaVmin, tSPM.xX.S.RtaCmin] = min(tSPM.xX.t(1,CHbR));
                           % [tSPM.xX.S.RtbVmax, tSPM.xX.S.RtbCmax] = max(tSPM.xX.t(V2r,CHbR));
                           % [tSPM.xX.S.RtbVmin, tSPM.xX.S.RtbCmin] = min(tSPM.xX.t(V2r,CHbR));

                            
                        catch exception
                            disp(exception.identifier)
                            disp(exception.stack(1));
                            disp('Problem generating OtaVmax, etc.');
                        end
                        %Add piece of SPM to the whole SPM
                        try
                            if markers_available
                                tSPM.xX.dur = (ei(iSubSess)-si(iSubSess))/fs;
                            end
                        end
                        iSPM = iSPM + 1;
                        SPM.xXn{iSPM} = tSPM.xX;
                        %save
                        save_dataON = 1; %we need the residuals for statistics calculations and the
                        %filtered data is often useful too
                        if save_dataON
                            try
                                if iSPM == 1
                                    nlst = lst;
                                end
                                temp  = tSPM.KY';
                                outfile = fullfile(spm_dir,['Sess' int2str(iSPM) '.nir']);
                                %Careful, this data may include HbT, therefore may
                                %have 50% more channels than the user expected...
                                fwrite_NIR(outfile,temp(:));
                                SPM.xY.Pf{iSPM,1} = outfile; %filtered
                                res_outfile = fullfile(spm_dir,['res_Sess' int2str(iSPM) '.nir']);
                                fwrite_NIR(res_outfile,res(:));
                                SPM.xXn{iSPM}.res = res_outfile;
                                SPM.xY.Cf = size(temp,1); %number of filtered channels stored
                                NIRS.Dt.fir.pp(nlst+1).p{iSPM,1} = outfile;
                            catch exception
                                disp(exception.identifier);
                                disp(exception.stack(1));
                                disp('Failed to save residuals or filtered data');
                            end
                        end
                    end %end for 1:nSubSess
                end
                %group weighted and unweighted amplitudes
                try
                    tmpa2 = zeros(size(SPM.xXn{1}.a2));
                    tmpb2 = zeros(size(SPM.xXn{1}.b2));
                    for n1=1:iSPM
                        %unweighted
                        tmpa2 = tmpa2 + SPM.xXn{n1}.a2;
                        tmpb2 = tmpb2 + SPM.xXn{n1}.b2;
                    end
                    tmpa2 = tmpa2/iSPM;
                    tmpb2 = tmpb2/iSPM;
                    %weighted
                    SPM.Gr.a2_uw = tmpa2;
                    SPM.Gr.b2_uw = tmpb2;
                    
                    %only do calculation if trRV, trRVRV and erdf are available
                    if ~isnan(SPM.xXn{1}.trRV) && ~(SPM.xXn{1}.trRV == 0)
                        nch = length(SPM.xXn{1}.ResSS);
                        cova = zeros(iSPM,nch);
                        covb = zeros(iSPM,nch);
                        betaa = zeros(1,nch);
                        betab = zeros(1,nch);
                        %fill data
                        for n1=1:iSPM
                            cova(n1,:) = SPM.xXn{n1}.ResSS./SPM.xXn{n1}.trRV*SPM.xXn{n1}.Bcov(1,1);
                            covb(n1,:) = SPM.xXn{n1}.ResSS./SPM.xXn{n1}.trRV*SPM.xXn{n1}.Bcov(V2r,V2r);
                            betaa(n1,:) = SPM.xXn{n1}.beta(1,:);
                            betab(n1,:) = SPM.xXn{n1}.beta(V2r,:);
                        end
                        %group - 1st Volterra
                        numera = betaa./cova;
                        SPM.Gr.betaa = sum(numera);
                        denuma = sqrt(sum(1./cova));
                        SPM.Gr.ta = SPM.Gr.betaa./denuma;
                        %2nd Volterra
                        numerb = betab./covb;
                        SPM.Gr.betab = sum(numerb);
                        denumb = sqrt(sum(1./covb));
                        SPM.Gr.tb = SPM.Gr.betab./denumb;
                    end
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1));
                    disp('Problem generating betaa, etc.');
                end
                try
                    K = SPM.xX.K;
                    K = rmfield(K, 'X');
                    K = rmfield(K, 'KL');
                    SPM.xX.K = K;
                end
                
                try
                    %add duration of subsessions, for convenience
                    if markers_available
                        SPM.dur = [];
                        for i1=1:length(SPM.xXn)
                            SPM.dur = [SPM.dur SPM.xXn{i1}.dur];
                        end
                    end
                end
                
                save(fGLM{g}, 'SPM');
                if save_dataON
                    save(job.NIRSmat{Idx,1},'NIRS');
                end
            end
            NIRS.flags.GLMest_OK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1));
        disp(['Could not estimate GLM for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;