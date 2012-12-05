function SPM = nirs_liom_average(NIRS,SPM)
try
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    iSPM = 0;
    try
        SPM.xY.P;
    catch
        disp(['Could not find data file for subject ' int2str(Idx)]);
    end
    
    %loop over sessions
    nsess = length(SPM.xY.P);
    first_pass = 1;
    %for the case of baseline_choice,'baseline_block_whole_session'
    %need to obtain and store the data from this session before proceeding
    if isfield(SPM.job.baseline_choice,'baseline_block_whole_session')
        base_choice = 2;
        baseline_session = SPM.job.baseline_choice.baseline_block_whole_session.baseline_session;
        srun = [baseline_session 1:nsess];
        storeY = 1;
    else
        if isfield(SPM.job.baseline_choice,'unique_baseline')
            base_choice = 3;
            baseline_session = SPM.job.baseline_choice.unique_baseline.baseline_session;
            srun = [baseline_session 1:nsess];
            storeY = 1;
        else
            base_choice = 1;
            srun = 1:nsess;
            storeY = 0;
        end
    end
    for s=srun
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
            tSPM.Idx = SPM.Idx;
            tSPM.Sess = SPM.Sess(s);
            tSPM.xX = SPM.xX;
            tSPM.fs = fs;
            tSPM.job = SPM.job;
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
            try
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
                tSPM.xX.K.row = 1:length(svec);
            end
            try
                try
                    beta = [SPM.Sess(s).col SPM.Sess(nsess).col(end)+s];
                catch
                    beta = 1;
                end
                tSPM.xX.X = SPM.xX.X(svec,beta);
            catch
                %Did not correctly extract this session's design
                %matrix?
            end
            
            %onsets
            tU = SPM.Sess(s).U;
            if markers_available
                for u0=1:length(tU)
                    %subtract time of subsession
                    tU(u0).ons = tU(u0).ons-si(iSubSess);
                    %remove negative onsets and onsets beyond end of subsession
                    id0 = tU(u0).ons > 0 & tU(u0).ons < (ei(iSubSess)-si(iSubSess));
                    tU(u0).ons = tU(u0).ons(id0);
                    tU(u0).dur = tU(u0).dur(id0);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Filtering and averaging
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Avg = [];
            if storeY
                if (s == baseline_session) && first_pass
                    baselineY = Y;
                else
                    tSPM.baselineY = baselineY;
                    [tSPM Avg] = averaging_core(tSPM,Y,tU);
                end
            else
                [tSPM Avg] = averaging_core(tSPM,Y,tU);
            end
            
            %Add piece of SPM to the whole SPM
            try
                if markers_available
                    tSPM.xX.dur = (ei(iSubSess)-si(iSubSess))/fs;
                end
            end
            
            if first_pass
                first_pass = 0;
                if storeY
                    %storeY = 0;
                else
                    iSPM = iSPM + 1;
                    if (storeY == 0 || (storeY == 1 && ((~(s == baseline_session) && base_choice == 2) || base_choice == 3)))
                        SPM.xXn{iSPM} = tSPM.xX;
                    end
                end
            else
                iSPM = iSPM + 1;
                if (storeY == 0 || (storeY == 1 && ((~(s == baseline_session) && base_choice == 2) || base_choice == 3)))
                    SPM.xXn{iSPM} = tSPM.xX;
                end
            end
            if iSPM == 1
                xX0 = tSPM.xX;
            end
            
            %save
            save_dataON = 1; %we need the residuals for statistics calculations and the
            %filtered data is often useful too
            if save_dataON && ~isempty(Avg) && (storeY == 0 || (storeY == 1 && ((~(s == baseline_session) && base_choice == 2) || base_choice == 3)))
                try
                    spm_dir = NIRS.spm_dir;
                    %if iSPM == 1
                    nlst = lst;
                    %end
                    temp  = Avg.KY';
                    Avg = rmfield(Avg,'KY');
                    %if storeY
                    %    outfile = fullfile(spm_dir,['Sess' int2str(iSPM-1) '.nir']);
                    %else
                    outfile = fullfile(spm_dir,['Sess' int2str(iSPM) '.nir']);
                    %end
                    %Careful, this data may include HbT, therefore may
                    %have 50% more channels than the user expected...
                    fwrite_NIR(outfile,temp(:));
                    SPM.xY.Pf{iSPM,1} = outfile; %filtered
                    %                     res_outfile = fullfile(spm_dir,['res_Sess' int2str(iSPM) '.nir']);
                    %                     fwrite_NIR(res_outfile,res(:));
                    %                     SPM.xXn{iSPM}.res = res_outfile;
                    %if storeY
                    %    Avg_outfile = fullfile(spm_dir,['Avg' int2str(iSPM-1) '.mat']);
                    %else
                    Avg_outfile = fullfile(spm_dir,['Avg' int2str(iSPM) '.mat']);
                    %end
                    save(Avg_outfile,'Avg');
                    SPM.xY.Cf = size(temp,1); %number of filtered channels stored
                    NIRS.Dt.fir.pp(nlst+1).p{iSPM,1} = outfile;
                    %******************************************************
                    %To save the pseudo-residuals to file
                    %Ke Peng, 2012-07-17
                    %******************************************************
                    res_outfile = fullfile(spm_dir,['res_Sess' int2str(iSPM) '.nir']);
                    fwrite_NIR(res_outfile,Avg.res(:));
                    SPM.xXn{iSPM}.res = res_outfile;
                    SPM.xXn{iSPM}.AvgFile = Avg_outfile;
                    %******************************************************                    
                catch exception
                    disp(exception.identifier);
                    disp(exception.stack(1));
                    disp('Failed to save filtered data');
                end
            end
        end %end for 1:nSubSess
    end
    % try
    %     K = SPM.xX.K;
    %     K = rmfield(K, 'X');
    %     K = rmfield(K, 'KL');
    %     SPM.xX.K = K;
    % end
    
    SPM.xX = xX0;
    try
        %add duration of subsessions, for convenience
        if markers_available
            SPM.dur = [];
            for i1=1:length(SPM.xXn)
                SPM.dur = [SPM.dur SPM.xXn{i1}.dur];
            end
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1));
    disp(['Could not do average']);
end
