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
    %filtering of the data
    KY = spm_filter_HPF_LPF_WMDL(SPM.xX.K, Y);
    if nreg > 1
        %Further filtering of the confounds
        KY = KY - SPM.xX.xKXs.X * (SPM.xX.pKX * KY);
    end
    SPM.KY = KY; %Store
    
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
        end
    end
    %averaging
    averaging_choice = SPM.job.averaging_choice;
    if isfield(averaging_choice,'block_averaging')
        avg_choice = 1;
        onset_delay = averaging_choice.block_averaging.onset_delay;
        onset_duration = averaging_choice.block_averaging.onset_duration;
        ons_delay = round(onset_delay*fs);
        ons_duration = round(onset_duration*fs);
    end
    
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
        end
        %data for onset after
        da{u0} = zeros(nons,nch);
        %subtracted data:
        ds{u0} = zeros(nons,nch);
        good_onsets = [];
        bad_onsets = [];
        for k0=1:nons
            co =round(U(u0).ons(k0)*fs); %current onset
            k0_OK = 1;
            try
                switch base_choice
                    case 1
                        tb{u0}(:,k0,:) = KY((co-base_offset-base_duration+1):(co-base_offset),:);
                    case 2 
                        if k0 == 1 && u0 == 1
                        end_data = min(size(SPM.baselineY,1),1-base_offset+base_duration);
                        if end_data < 1-base_offset+base_duration
                            disp(['Working with truncated data for baseline' ...
                                ' for stim type ' int2str(u0) ' for subject Idx=' int2str(SPM.Idx)]);
                        end
                        tb = SPM.baselineY((1-base_offset):end_data,:);
                        end
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
                %temporal averaging for each onset
                da{u0}(k0,:) = mean(ta{u0}(:,k0,:),1);
                switch base_choice
                    case 1
                        db{u0}(k0,:) = mean(tb{u0}(:,k0,:),1);
                        ds{u0}(k0,:) = da{u0}(k0,:) - db{u0}(k0,:);
                        ts{u0}(:,k0,:) = ta{u0}(:,k0,:) - permute(repmat(db{u0}(k0,:),[ons_duration 1]),[1 3 2]);
                    case 2
                        db = mean(tb,1);
                        ds{u0}(k0,:) = da{u0}(k0,:) - db(1,:);
                        ts{u0}(:,k0,:) = ta{u0}(:,k0,:) - permute(repmat(db(1,:),[ons_duration 1]),[1 3 2]);
                end
             end
        end
        %averaging over onsets
        a(u0,:) = mean(ds{u0}(good_onsets,:),1);
        s(u0,:) = std(ds{u0}(good_onsets,:),1);
        Avg.go{u0} = good_onsets;
        Avg.bo{u0} = unique(bad_onsets);
    end
    SPM.xX.beta = a;
    SPM.xX.covbeta = s;
    SPM.xX.t = a./s;
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
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem with averaging');
    Avg = [];
end