function out = nirs_run_wls_bglm_estimate(job)
% NIRS_SPM GLM estimation - first level

% Which GLM to estimate {1='First', 2='All', 3='Last'} of all specified GLMs;
which_GLM = job.NIRS_SPM_which_GLM;

% Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        try
            switch which_GLM
                case 1 % first
                    fGLM = NIRS.SPM(1);
                case 2
                    % need to loop (all GLMs specified)
                    fGLM = NIRS.SPM;
                case 3 % last
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
            try
                load(fullfile(fGLM{g},'SPM.mat'));
            end

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

                    %Y = Y(:,[25 85]);
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
                    if SPM.xX.LPFbutter                
                        cutoff=SPM.xX.lpf_butter_freq; %0.666; %Hz, or 1.5s
                        FilterOrder=3;
                        Wn=cutoff*2/fs;                           % normalised cutoff frequency
                        [fb,fa]=butter(FilterOrder,Wn);            % buterworth filter
                        Y=filtfilt(fb,fa,Y);                            
                    end

                    %HPF
                    if SPM.xX.HPFbutter
                        cutoff=SPM.xX.hpf_butter_freq; %Hz, or 100s 
                        FilterOrder=3;
                        Wn=cutoff*2/fs;
                        [fb,fa]=butter(FilterOrder,Wn,'high');
                        Y=filtfilt(fb,fa,Y);
                        %need to filter the design matrix too,
                        %otherwise, the estimates will be significantly
                        %biased 

                        %This is done later
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
                        %find elements of X for session s
                        nbeta = size(SPM.xX.X,2);
                        nbetaS = (nbeta-nsess)/nsess;
                        %last entry is the constant regressor
                            beta = [(s-1)*nbetaS+1:s*nbetaS nbeta-nsess+s];

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
                            if SPM.xX.HPFbutter
                                %filter the design matrix
                                cutoff=SPM.xX.hpf_butter_freq; %Hz, or 100s 
                                FilterOrder=3;
                                Wn=cutoff*2/fs;
                                [fb,fa]=butter(FilterOrder,Wn,'high');
                                %exclude the constant
                                tX=filtfilt(fb,fa,tSPM.xX.X(:,1:end-1));
                                %add back the constant
                                tSPM.xX.X = [tX tSPM.xX.X(:,end)];
                            end
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
                            if size(tSPM.xX.X,2) > 1
                                tmpX = tSPM.xX.X(:,1:end-1);
                            else
                                tmpX = tSPM.xX.X;
                            end
                            [Betas,spectralExponents,Modelisation,Design] = ...
                                wls(fs,Y,tmpX,tSPM.xX.opt.Design);
                            tSPM.xX.beta = Betas;

                        case 'NIRS_SPM'
                             if precolor
                                tSPM = precoloring_batch(tSPM,Y);                                 
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
                            tSPM.xX.t=tSPM.xX.beta./sqrt(tSPM.xX.Bvar);
                        case 'NIRS_SPM'
                            for r=1:size(tSPM.xX.beta,1)
                                tSPM.xX.t(r,:) = tSPM.xX.beta(r,:)./...
                                    sqrt(tSPM.xX.ResSS(:)'*tSPM.xX.Bcov(r,r)/tSPM.xX.trRV);
                            end
                            %b2 =SPM.xXn{1}.beta(1,:)./sqrt(SPM.xXn{1}.ResSS/SPM.xXn{1}.trRV);
                    end
                    %A posteriori estimation of the amplitude of the response
                    try
                        %filtered data less estimated beta times filtered design
                        %matrix
                        sigma = std(Y- (tSPM.xX.X(:,1)*tSPM.xX.beta(1,:) + tSPM.xX.X(:,2)*tSPM.xX.beta(2,:)),0,1);
                        tSPM.xX.a2 = (tSPM.xX.beta(1,:) ./ sigma)/4.7506; %norm_bf = 4.7506 = 1/max(X(:,1)) for a protocol with only one spike

                        tSPM.xX.b  = tSPM.xX.beta(2,:) ./ tSPM.xX.beta(1,:);
                    catch
                    end
                    %Add piece of SPM to the whole SPM 
                    try 
                        if markers_available
                            tSPM.xX.dur = (ei(iSubSess)-si(iSubSess))/fs; 
                        end
                    end
                    iSPM = iSPM + 1;
                    SPM.xXn{iSPM} = tSPM.xX;
                end %end for 1:nSubSess
            end 
            try
                K = SPM.xX.K;
                K = rmfield(K, 'X');
                K = rmfield(K, 'KL');
                SPM.xX.K = K;
                %clear K;
            end
%             switch SPM.xX.opt.meth
%                 case {'BGLM', 'WLS'}
%                     SPM.beta      = 0; 
%                     SPM.xX.Bvar   = [];
%                     for n1=1:length(SPM.xXn)
%                         SPM.xX.Bvar   = blkdiag(SPM.xX.Bvar,sparse(SPM.xXn{iSPM}.Bvar));
%                         SPM.beta      = SPM.beta + SPM.xXn{iSPM}.beta; %incorrect
%                     end
%                     SPM.xX.t = SPM.beta ./sqrt(SPM.xX.Bvar);
%                 case 'NIRS_SPM'
%                     %Reconstruct SPM
%                     SPM.xX.trRV   = 0;
%                     SPM.xX.trRVRV = 0;
%                     SPM.xX.Bcov   = [];
%                     SPM.beta      = 0;
%                     SPM.ResSS     = 0;
%                     for n1=1:length(SPM.xXn)
%                         SPM.xX.trRV   = SPM.xX.trRV   + SPM.xXn{iSPM}.trRV;
%                         SPM.xX.trRVRV = SPM.xX.trRVRV + SPM.xXn{iSPM}.trRVRV;
%                         SPM.xX.Bcov   = blkdiag(SPM.xX.Bcov,sparse(SPM.xXn{iSPM}.Bcov));
%                         SPM.beta      = SPM.beta+ SPM.xXn{iSPM}.beta; %incorrect
%                         SPM.ResSS     = SPM.ResSS + SPM.xXn{iSPM}.ResSS;
%                     end
%                     SPM.xX.erdf = (SPM.xX.trRV)^2/SPM.xX.trRVRV;
%                     SPM.xX.t = zeros(size(SPM.beta));
%                     for r=1:size(SPM.beta,1)
%                         SPM.xX.t(r,:) = SPM.beta(r,:)./sqrt(SPM.ResSS(:)'*SPM.xX.Bcov(r,r)/SPM.xX.trRV);
%                     end
% 
%             end
            try 
                %add duration of subsessions, for convenience
                if markers_available
                    SPM.dur = [];
                    for i1=1:length(SPM.xXn)
                        SPM.dur = [SPM.dur SPM.xXn{i1}.dur]; 
                    end
                end
            end

            save(fullfile(fGLM{g},'SPM.mat'), 'SPM');

        end
    catch
        disp(['Could not estimate GLM for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;