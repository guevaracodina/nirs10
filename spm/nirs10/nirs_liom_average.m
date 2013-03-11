function [SPM NIRS] = nirs_liom_average(NIRS,SPM)
try
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    iSPM = 0;
    try
        SPM.xY.P;
    catch
        disp(['Could not find data file for subject ' int2str(Idx)]);
    end
    if SPM.GenerateHbT
        NCt = 3/2*NC;
    else
        NCt = NC;
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
            if isfield(SPM.job.baseline_choice,'no_baseline')
                base_choice = 4;
                srun = 1:nsess;
                storeY = 0;               
                a = zeros(length(srun),NCt);
            else
                if isfield(SPM.job.baseline_choice,'baseline_average_sessions')
                    base_choice = 5;
                    baseline_session = SPM.job.baseline_choice.baseline_average_sessions.sessions_to_average;
                    srun = [baseline_session 1:nsess];
                    nBaseline = length(baseline_session);
                    first_pass = length(baseline_session);
                    baselineB = zeros(1,NCt);
                    B_start = SPM.job.baseline_choice.baseline_average_sessions.SL_start;
                    B_duration = SPM.job.baseline_choice.baseline_average_sessions.SL_duration;
                    B_end = SPM.job.baseline_choice.baseline_average_sessions.SL_end;   
                    B_start = nirs_rep_array(B_start,nBaseline); 
                    B_duration = nirs_rep_array(B_duration,nBaseline);
                    B_end = nirs_rep_array(B_end,nBaseline);
                    storeY = 1;
                else
                    base_choice = 1;
                    srun = 1:nsess;
                    storeY = 0;
                end
            end
        end
    end
    sc = 0; %session counter for baseline -- only for base_choice = 5
    sc2 = 0; %session counter for sessions to average -- only for option combine_sessionsEXT
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
            
            if isfield(SPM.job.averaging_choice,'average_all_data')
                need_onsets = 0;
                need_res = 0;
                tU = [];
            else
                need_res = 1;
                need_onsets = 1;
            end
            if need_onsets
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
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Filtering and averaging
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Avg = [];
            if storeY
                if any(s == baseline_session) && first_pass
                    sc = sc + 1;
                    if base_choice == 2 || base_choice == 3
                        baselineY = Y;
                    else
                        if base_choice == 5
                            %store a single value per channel                            
                            if isempty(B_end)
                                Yb = Y;
                            else
                                Yb = Y(1:end-round(fs*B_end(sc)),:);
                            end
                            if ~isempty(B_start)
                                Yb = Yb(round(fs*B_start(sc)):end,:);
                            end
                            if ~isempty(B_duration)
                                if size(Yb,1) > round(fs*B_duration(sc))                                                                 
                                    Yb = Yb(1:round(fs*B_duration(sc)),:);
                                end
                            end
                            baselineB = baselineB + mean(Yb,1)/nBaseline;
                        end
                    end
                else
                    sc2 = sc2 + 1;
                    try tSPM.baselineY = baselineY; end
                    try tSPM.baselineB = baselineB; end
                    Ya = Y;
                    if isfield(SPM.job.averaging_choice,'average_all_data')
                        if isfield(SPM.job.averaging_choice.average_all_data,'combine_sessionsEXT')
                            SL = SPM.job.averaging_choice.average_all_data.combine_sessionsEXT.session_list2EXT;
                            %loop through SL to find start/duration/end for this session
                            for s0=1:length(SL)
                                sl0 = SL(s0).SL_List;
                                sl = find(sc2 == sl0);
                                A_start = '';
                                A_duration = '';
                                A_end = '';
                                if ~isempty(sl)
                                    A_start = SL(s0).SL_start;
                                    A_duration = SL(s0).SL_duration;
                                    A_end = SL(s0).SL_end;
                                    nSL = length(sl0);
                                    A_start = nirs_rep_array(A_start,nSL);
                                    A_duration = nirs_rep_array(A_duration,nSL);
                                    A_end = nirs_rep_array(A_end,nSL);
                                end
                                if ~isempty(A_end)
                                    Ya = Ya(1:end-round(fs*A_end(sl)),:);
                                end
                                if ~isempty(A_start)
                                    Ya = Ya(round(fs*A_start(sl)):end,:);
                                end
                                if ~isempty(A_duration)
                                    if size(Ya,1) > round(fs*A_duration(sl))
                                        Ya = Ya(1:round(fs*A_duration(sl)),:);
                                    end
                                end
                            end
                        end
                    end                    
                    [tSPM Avg] = averaging_core(tSPM,Ya,tU);
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
                first_pass = first_pass-1;
                if storeY
                    %storeY = 0;
                else
                    iSPM = iSPM + 1;
                    if (storeY == 0 || (storeY == 1 && ((~any(s == baseline_session) && (base_choice == 2 )) || base_choice == 3)))
                        SPM.xXn{iSPM} = tSPM.xX;
                    end
                end
            else
                iSPM = iSPM + 1;
                if (storeY == 0 || (storeY == 1 && ((~any(s == baseline_session) && (base_choice == 2 )) || base_choice == 3 || base_choice == 5)))
                    SPM.xXn{iSPM} = tSPM.xX;
                end
            end
            if iSPM == 1
                xX0 = tSPM.xX;
            end
            save_dataON = 1; %we need the residuals for statistics calculations and the
            %filtered data is often useful too
            if save_dataON && ~isempty(Avg) && (storeY == 0 || (storeY == 1 && ((~any(s == baseline_session) && (base_choice == 2 )) || base_choice == 3 || base_choice ==5)))
                try
                    spm_dir = NIRS.spm_dir;
                    nlst = lst;
                    temp  = Avg.KY';
                    Avg = rmfield(Avg,'KY');
                    outfile = fullfile(spm_dir,['Sess' int2str(iSPM) '.nir']);
                    fwrite_NIR(outfile,temp(:));
                    SPM.xY.Pf{iSPM,1} = outfile; %filtered
                    if base_choice == 4 || base_choice == 5
                        a(iSPM,:) = Avg.a;
                        NIRS.Avg.a = a;
                    else
                        Avg_outfile = fullfile(spm_dir,['Avg' int2str(iSPM) '.mat']);
                        %end
                        save(Avg_outfile,'Avg');
                        SPM.xXn{iSPM}.AvgFile = Avg_outfile;
                    end
                    SPM.xY.Cf = size(temp,1); %number of filtered channels stored
                    NIRS.Dt.fir.pp(nlst+1).p{iSPM,1} = outfile;
                    if need_res
                        res_outfile = fullfile(spm_dir,['res_Sess' int2str(iSPM) '.nir']);
                        fwrite_NIR(res_outfile,Avg.res(:));
                        SPM.xXn{iSPM}.res = res_outfile;
                    end                    
                catch exception
                    disp(exception.identifier);
                    disp(exception.stack(1));
                    disp('Failed to save filtered data');
                end
            end
        end %end for 1:nSubSess
    end

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
    %
    averaging_choice = SPM.job.averaging_choice;
    do_SLavg = 0;
    if isfield(averaging_choice,'average_all_data')
        if isfield(averaging_choice.average_all_data,'combine_sessions')
            SL = averaging_choice.average_all_data.combine_sessions.session_list2;
            do_SLavg = 1;
        else
            if isfield(averaging_choice.average_all_data,'combine_sessionsEXT')
                SL = averaging_choice.average_all_data.combine_sessionsEXT.session_list2EXT;
                do_SLavg = 1;
            end
        end
        if do_SLavg
            %combine the sessions 
            oldSPM = SPM;
            SPM = [];
            SPM.xY.Cf = oldSPM.xY.Cf;
            iSPM = 0;           
            for s0=1:length(SL)
                try 
                    iSPM = iSPM+1;
                    temp = oldSPM.xXn{SL(s0).SL_List(1)}.beta;
                    nSess = length(SL(s0).SL_List);
                    for k0=2:nSess
                        temp = temp + oldSPM.xXn{SL(s0).SL_List(k0)}.beta;
                    end
                    temp = temp/nSess;
                    SPM.xXn{iSPM}.beta = temp;
                    SPM.xXn{iSPM}.Sname = SL(s0).SL_Label;
                    SPM.xXn{iSPM}.X = oldSPM.xXn{1}.X; 
                catch
                    iSPM = iSPM-1;
                    disp(['Could not process session list ' int2str(s0)])
                end
            end
        end
    end    
catch exception
    disp(exception.identifier)
    disp(exception.stack(1));
    disp(['Could not do average']);
end
