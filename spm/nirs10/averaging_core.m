function [SPM Avg] = averaging_core(SPM,Y,U)
try
    [nScan nreg] = size(SPM.xX.X);
    fs = SPM.fs;
    
    if nreg > 1
        % Set up of the filters / Design space and projector matrix [pseudoinverse] for WLS
        %==================================================================
        switch SPM.xX.K.HParam.type
            case 'Wavelet-MDL'
                tmp_K = SPM.xX.K;
                tmp_K.HParam.type = '';
                SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(tmp_K, SPM.xX.X)); % KX
                SPM.xX.K.X = SPM.xX.X;
                clear tmp_K;
            case 'DCT'
                SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.xX.X)); % KX
            case 'none'
                SPM.xX.xKXs = spm_sp('Set', spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.xX.X)); % KX ?
        end
        SPM.xX.xKXs.X = full(SPM.xX.xKXs.X);
        SPM.xX.pKX = spm_sp('x-', SPM.xX.xKXs); % projector
    end
        
    %baseline
    baseline_choice = SPM.job.baseline_choice;
    if isfield(baseline_choice,'baseline_block_averaging')
        base_choice = 1;
        baseline_offset = baseline_choice.baseline_block_averaging.baseline_offset;
        baseline_duration = baseline_choice.baseline_block_averaging.baseline_duration;
        base_offset = round(baseline_offset*fs);
        base_duration = round(baseline_duration*fs);
    else
        if isfield(baseline_choice,'baseline_block_whole_session')
            base_choice = 2;
            %baseline_session = baseline_choice.baseline_block_whole_session.baseline_session;
            baseline_offset = baseline_choice.baseline_block_whole_session.baseline_offset;
            base_offset = round(baseline_offset*fs);
            baseline_duration = baseline_choice.baseline_block_whole_session.baseline_duration;
            base_duration = round(baseline_duration*fs);
        else
            if isfield(baseline_choice,'unique_baseline')
                base_choice = 3;
                baseline_start = baseline_choice.unique_baseline.baseline_start;
                base_start = round(baseline_start*fs);
                baseline_duration = baseline_choice.unique_baseline.baseline_duration;
                base_duration = round(baseline_duration*fs);
            else
                if isfield(baseline_choice,'no_baseline')
                    base_choice = 4;
                else
                    if isfield(baseline_choice,'baseline_average_sessions')
                        base_choice = 5;
                    end
                end
            end
        end
    end
    
    %filtering of the data
    if base_choice ~= 5
        KY = spm_filter_HPF_LPF_WMDL(SPM.xX.K, Y);
        if isfield(SPM,'baselineY')
            baselineY = spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.baselineY);
        end
        if nreg > length(U)+1
            %Further filtering of the confounds
            KY = KY - SPM.xX.xKXs.X * (SPM.xX.pKX * KY);
        end
    else
        KY = Y;
    end
    SPM.KY = KY; %Store
    
    
    %averaging
    averaging_choice = SPM.job.averaging_choice;
    if isfield(averaging_choice,'block_averaging')
        avg_choice = 1;
        onset_delay = averaging_choice.block_averaging.onset_delay;
        onset_duration = averaging_choice.block_averaging.onset_duration;
        %Quick fix for Sarah's study -- if NIRS files are too short, set
        %duration from 120 s to 60 s (due to protocol duration changing
        %according to condition (COG: 60 s, MOT, TDT: 120 s).
        if size(KY,1)/fs < 110 %will affect only very short NIRS files of about 1 minute or less
            onset_duration = min(55,size(KY,1)/fs-20); %subtract 1.5 s for protocol to start, plus 5 seconds of cross on screen, plus 5 seconds for beginning of hemodynamic response
        end
        ons_delay = round(onset_delay*fs);
        ons_duration = round(onset_duration*fs);
    else
        if isfield(averaging_choice,'average_all_data')
            avg_choice = 2;
        else
            if isfield(averaging_choice,'seizure_evolution')
               avg_choice = 3; 
               onset_delays = averaging_choice.seizure_evolution.onset_delays;
               onset_duration = averaging_choice.seizure_evolution.onset_duration;
               ons_delays = round(onset_delays*fs);
               ons_duration = round(onset_duration*fs);
            end
        end
    end
    switch avg_choice
        case 3 %seizure_evolution
            LOD = length(onset_delays);
            nst = length(U.ons)*LOD;
            nch = size(KY,2);
            SPM.xX.beta = zeros(nst,nch);%stimulus types by channels
            SPM.xX.covbeta = zeros(nst,nch);
            for u0=1:length(U.ons)
                %for covbeta
                switch base_choice
                    case 1                       
                        bo = round(fs*(U.ons(u0)-baseline_offset));
                        bidx = (bo-base_duration):bo;
                    case 3
                        bidx = base_start:(base_start+base_duration);
                end
                B0 = std(KY(bidx,:),0,1).^2;
                M0 = mean(KY(bidx,:),1);
                for t0=1:LOD
                    idx1 = round(fs*(U.ons(u0)+onset_delays(t0)));
                    idx = idx1:(idx1+ons_duration-1);
                    SPM.xX.beta(t0+(u0-1)*LOD,:) = mean(KY(idx,:),1)-M0;
                    SPM.xX.covbeta(t0+(u0-1)*LOD,:) = B0;                    
                end
            end
            SPM.xX.t = SPM.xX.beta./SPM.xX.covbeta.^(1/2);
            Avg.a = SPM.xX.beta;
            Avg.b = B0; %not quite right
            Avg.KY = KY;
            %pseudo-residuals -- ?
        case 2 %basic, extended            
            a = mean(KY,1);
            Avg.a = a;
            Avg.b = SPM.baselineB;
            Avg.KY = KY;
            SPM.xX.beta = a-Avg.b;
        case 1
            nst = length(U);
            nch = size(KY,2);
            SPM.xX.beta = zeros(nst,nch);%stimulus types by channels
            SPM.xX.covbeta = zeros(nst,nch);
            
            tb{nst} = []; ta{nst} = [];
            db{nst} = []; da{nst} = []; d{nst} = [];
            a = zeros(nst,nch); %average
            s = zeros(nst,nch); %standard deviation
            for u0=1:nst
                nons = length(U(u0).ons);
                %data: array after
                ta{u0} = zeros(ons_duration,nons,nch);
                %subtracted data
                ts{u0} = zeros(ons_duration,nons,nch);
                switch base_choice
                    case 1
                        %baseline: array before
                        tb{u0} = zeros(base_duration,nons,nch);
                        %baseline before
                        db{u0} = zeros(nons,nch);
                    case 2
                        tb = zeros(base_duration,nch);
                        db = zeros(1,nch);
                    case 3
                        if u0 == 1
                            tb = baselineY((base_start+1):(base_start+base_duration),:);
                            db = mean(tb,1);
                        end
                end
                %data for onset after
                da{u0} = zeros(nons,nch);
                %subtracted data:
                ds{u0} = zeros(nons,nch);
                good_onsets = [];
                bad_onsets = [];
                %******************************************************************
                %To create an ideal averaged response
                %Ke Peng. 2012-07-20
                %******************************************************************
                V{u0} = zeros(size(KY,1),nch);
                response_good_onsets = [];
                %******************************************************************
                
                for k0=1:nons
                    co =round(U(u0).ons(k0)*fs); %current onset
                    %co
                    k0_OK = 1;
                    try
                        switch base_choice
                            case 1
                                tb{u0}(:,k0,:) = KY((co-base_offset-base_duration+1):(co-base_offset),:);
                            case 2
                                if k0 == 1 && u0 == 1
                                    end_data = min(size(baselineY,1),1-base_offset+base_duration);
                                    if end_data < 1-base_offset+base_duration
                                        disp(['Working with truncated data for baseline' ...
                                            ' for stim type ' int2str(u0) ' for subject Idx=' int2str(SPM.Idx)]);
                                    end
                                    tb = baselineY((1-base_offset):end_data,:);
                                end
                            case 3
                        end
                    catch
                        disp(['Could not include baseline ' int2str(k0) ' for stim type ' int2str(u0) ' for subject Idx=' int2str(SPM.Idx)]);
                        k0_OK = 0;
                        bad_onsets = [bad_onsets k0];
                    end
                    try
                        %                 end_data = min(co+ons_delay+ons_duration,size(KY,1));
                        %                 if end_data < co+ons_delay+ons_duration
                        %                     disp(['Working with truncated data for onset ' int2str(k0) ...
                        %                         ' for stim type ' int2str(u0) ' for subject Idx=' int2str(SPM.Idx)]);
                        %                 end
                        ta{u0}(:,k0,:) = KY((co+ons_delay+1):co+ons_delay+ons_duration,:);
                    catch
                        disp(['Could not include onset ' int2str(k0) ' for stim type ' int2str(u0) ' for subject Idx=' int2str(SPM.Idx)]);
                        k0_OK = 0;
                        bad_onsets = [bad_onsets k0];
                    end
                    if k0_OK
                        good_onsets = [good_onsets k0];
                        %**********************************************************
                        %To restore the good onsets response for pseudo-residuals
                        %Ke Peng, 2012-07-16
                        %**********************************************************
                        %response_good_onsets = [response_good_onsets co+ons_delay+1];
                        %**********************************************************
                        
                        %temporal averaging for each onset
                        da{u0}(k0,:) = mean(ta{u0}(:,k0,:),1);
                        switch base_choice
                            case 1
                                tmp_baseline_array = tb{u0}(:,k0,:);
                                db{u0}(k0,:) = squeeze(mean(tmp_baseline_array,1));
                                ds{u0}(k0,:) = da{u0}(k0,:) - db{u0}(k0,:);
                                ts{u0}(:,k0,:) = ta{u0}(:,k0,:) - permute(repmat(db{u0}(k0,:),[ons_duration 1]),[1 3 2]);
                            case 2
                                db = mean(tb,1);
                                ds{u0}(k0,:) = da{u0}(k0,:) - db(1,:);
                                ts{u0}(:,k0,:) = ta{u0}(:,k0,:) - permute(repmat(db(1,:),[ons_duration 1]),[1 3 2]);
                            case 3
                                ds{u0}(k0,:) = da{u0}(k0,:) - db(1,:);
                                ts{u0}(:,k0,:) = ta{u0}(:,k0,:) - permute(repmat(db(1,:),[ons_duration 1]),[1 3 2]);
                        end
                    end
                    %**********************************************************
                    %To calculate the ideal averaged response
                    %Ke Peng, 2012-07-17
                    %**********************************************************
                    % for l0 = 1 : nch
                    %for d0 = 1 : (ons_duration)
                    % if V{u0}((co+ons_delay+d0),l0)  == 0
                    offset = co+ons_delay+1;
                    range = offset:min(offset+ons_duration-1,size(KY,1));
                    if offset+ons_duration -1 > size(KY,1)
                        adj_duration = size(KY,1) - offset + 1;
                    else
                        adj_duration = ons_duration;
                    end
                    V{u0}(range,:) = repmat(ds{u0}(k0,:),adj_duration,1);
                    %                     else
                    %                         V{u0}((co+ons_delay+d0),l0) = (V{u0}((co+ons_delay+d0),l0) + da{u0}(k0,l0))/2;
                    %      end
                    %    end
                    
                    %**********************************************************
                end
                %averaging over onsets
                
                %**********************************************************
                %To calculate the pseudo-residuals
                %Ke Peng, 2012-07-17
                %**********************************************************
                
                res = KY - V{u0};
                ResSS = sum(res.*2);
                ResSSch = res' * res; %V{u0}'*V{u0};
                
                SPM.xX.ResSS = ResSS;
                SPM.xX.ResSSch = ResSSch;
                Avg.res = res';
                
                %for l0 = 1 : nch
                %    V{u0}((co+ons_delay+1):(co+ons_delay+ons_duration),l0) = da{u0}(k0,:);
                %end
                %res = KY(response_good_onsets,:) - da{u0};
                %ResSS = sum(res.^2);
                %ResSSch = KY(response_good_onsets,:)' * KY(response_good_onsets,:) - da{u0}'*da{u0};
                %SPM.xX.ResSS = ResSS;
                %SPM.xX.ResSSch = ResSSch;
                %Avg.res = res';
                %**********************************************************
                tmean = mean(ds{u0}(good_onsets,:),1);
                tstd = std(ds{u0}(good_onsets,:),[],1);
                nfactor = length(good_onsets)^0.5;
                tstd(tstd<(abs(tmean)/nfactor)) = tmean(tstd<abs(tmean)/nfactor)/nfactor^4; %enforce std to be at least comparable to mean, depending on number of onsets
                a(u0,:) = mean(ds{u0}(good_onsets,:),1);
                s(u0,:) = tstd; %std(ds{u0}(good_onsets,:),[],1);
                
                Avg.go{u0} = good_onsets;
                Avg.bo{u0} = unique(bad_onsets);
            end
            Avg.KY = KY;
            SPM.xX.beta = a;
            SPM.xX.covbeta = s.^2;
            SPM.xX.t = a./(s/length(Avg.go{1})^(1/2)); %incorrect, one needs 
            %to divide s by the number of stims^0.5 -- this assumes a constant number of stims and that 
            %Avg.go{1} is this number
            %Fill Avg
            Avg.a = a;
            %     if sum(isnan(a(:)))
            %         disp(['Problem: NaN values for subject with Idx=' int2str(SPM.Idx)]);
            %     end
            Avg.s = s;
            Avg.ds = ds;
            Avg.ts = ts;
            Avg.da = da;
            Avg.ta = ta;
            Avg.db = db;
            Avg.tb = tb;
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem with averaging');
    Avg = [];
end