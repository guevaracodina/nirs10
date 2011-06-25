function out = nirs_run_addTestStimuli(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

outNIRSmat = {};
%Various options for the simulations
%The final options chosen are:
%GammaOn = 0: use mostly canonical HRF
%filter_X = 1: need to filter design matrix
%calculate_bf_norm = 0: no additional normalization
%std_or_power = 1; normalize with respect to std of filtered baseline
%normalize_impulse_area_to_unity = 0; should not be used, not giving what
%it was intended for
%save_all_data = 0; only for generating data for pedagogical figures
%set_baseline_to_zero = 0: only for testing, do not use, removes all the
%baseline data
try
    switch job.testGamma
        case 1
            GammaOn = 1;
        case 2
            GammaOn = 0;
    end
catch
    GammaOn = 0; %Use canonical HRF rather than gamma HRF
end
try
    filter_data = job.testFilterData;
catch
    filter_data = 1;
end
try
    filter_X = job.testFilterX;
catch
    filter_X = 1;
end
try
    calculate_bf_norm = job.testBfNorm;
catch
    calculate_bf_norm = 0;
end
try
    std_or_power = job.testStdvsPower;
catch
    std_or_power = 1;
end
try 
    use_Butter_HPF = job.testHPFButterOn;
catch
    use_Butter_HPF = 1;
end
try
    HPFButt = job.testHPFbutterCutoff;
catch
    HPFButt = 0.004;
end
try
    HPFButt_order = job.testHPFbutterOrder;
catch
    HPFButt_order = 3;
end
try 
    use_gaussian_LPF = job.testLPFGaussianOn;
catch
    use_gaussian_LPF = 1;
end
try 
    LPFGaussianFWHM = job.testLPFGaussianFWHM;
catch
    LPFGaussianFWHM = 1.5;
end
try
    use_wavelet_filter = job.testWaveletMDLOn;
catch
    use_wavelet_filter = 0;
end
normalize_impulse_area_to_unity = 0; %no longer used
remove_X_mean = 0; %no longer used
use_OLD_wrong_Butter_LPF = 0; %no longer used
LFP_butter = 0.67; %0.067 was a typo %no longer used
%To save data for figures: spikes simulated, response with and without noise
save_all_data = 0; 
%Ensure reproducibility of results by resetting the random seed
set_baseline_to_zero = 0;
for Idx=1:size(job.NIRSmat,1)
    
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});       
        %Parameters
        tc = job.testChannels;
        try
            job.testPType.testEP;
            %Paradigm type 0: event
            tp = 0;
            try
                job.testPType.testEP.NoFrequentSpikes;
                tss = job.testPType.testEP.NoFrequentSpikes.testExpSlowSpike1;
                tsc = job.testPType.testEP.NoFrequentSpikes.testRescaleOn1;
                rseed = job.testPType.testEP.NoFrequentSpikes.testSeed1;
                %Event paradigm type 0
                tp2 = 0;
            catch
                try
                    job.testPType.testEP.FrequentSpikes;
                    tfs = job.testPType.testEP.FrequentSpikes.testExpFastSpike;
                    tss = job.testPType.testEP.FrequentSpikes.testExpSlowSpike2;
                    tsc = job.testPType.testEP.FrequentSpikes.testRescaleOn2;
                    tsfg = job.testPType.testEP.FrequentSpikes.testAvgNumFastSpikes_perGroup;
                    tssg = job.testPType.testEP.FrequentSpikes.testAvgNumSlowSpikes_perGroup;
                    rseed = job.testPType.testEP.FrequentSpikes.testSeed2;
                    %Event paradigm type 1
                    tp2 = 1;
                catch
                    disp('unrecognized event paradigm, aborting');
                    return
                end
            end            
        catch           
            try
                job.testPType.testBP;
                %Paradigm type 1: block
                tp = 1;
            catch
                disp('unrecognized paradigm, aborting');
                return
            end
        end
        tn = job.testStimuliNumber;
        ts = job.testSessionNumber;
        tl = job.testWavelength;
        try
            ta = job.testAmplitudeTarget.testAmplitude/100;
            AmpTargetMethod = 1;
        catch
            tSNR = job.testAmplitudeTarget.testSNR;
            AmpTargetMethod = 0;
        end
        
        try tb = job.testAmplitude2/100; catch; end
        tname = job.testStimulusName;
        volt = job.voltAddStim;
        try
            job.keepAllChannels.AllChannels;
            %keeping all channels
            AllChannels = 1;
        catch
            AllChannels = 0;
            tk = job.keepAllChannels.keepChannels;
        end
        
        %NIRS total number of channels
        NC = NIRS.Cf.H.C.N;
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        try
            f=1;
            bpi = NIRS.Dt.fir.pp(lst).bpi{f,1}; %bad point indices
            bpd = NIRS.Dt.fir.pp(lst).bpd{f,1}; %bad point durations
            si = NIRS.Dt.fir.pp(lst).si{f,1};
            ei = NIRS.Dt.fir.pp(lst).ei{f,1};
            if ~isempty(bpi)
                markers_available = 1;
            else
                markers_available = 0;
            end
        catch
            bpi = [];
            bpd = [];
            markers_available = 0;
        end
        %Here could add loop over sessions - to be done?
        [dir1,fil1,ext1] = fileparts(NIRS.Dt.fir.pp(lst).p{ts,:});
        d = fopen_NIR(fullfile(dir1,[fil1,ext1]),NC,ext1);
        if job.testDupChannels
            %store a copy of the data
            d_copy = d;
        end
        %channels of interest
        chn = [];
        for i=1:length(tl)
            chn = [chn tc+(tl(i)-1)*NC/length(NIRS.Cf.dev.wl)];
        end
        
        dc = d(chn,:); %we keep d, as we will write dc over d, and save
        ns = size(d,2);
        %frequency
        fs = NIRS.Cf.dev.fs;
        
        %now repeat with more time precision to calculate bf
        xBF.T = 10;
        xBF.dt = 1/(fs*xBF.T); % - time bin length {seconds}
        if GammaOn
            xBF.name = 'Gamma functions';
            xBF.length = 16;
            xBF.order  = 1;
        else
            xBF.name = 'hrf'; %description of basis functions specified
        end
        xBF = spm_get_bf(xBF);
        if normalize_impulse_area_to_unity
            %unit impulse
            xBF.bf = xBF.bf/(xBF.dt*sum(xBF.bf));
        else
            %unit area
            xBF.bf = xBF.bf/sum(xBF.bf); %normalize
        end
        bf  = xBF.bf;
        
        %factor to normalize bf to 1
        if calculate_bf_norm
            %fill SPM's xBF structure from spm_get_bf %fill SPM's U structure from spm_get_ons
            SPM = [];
            %SPM's session
            s = 1;
            SPM.nscan(s) = ns;
            SPM.xBF = xBF;
            SPM.xBF.T = xBF.T; %get more precision on onsets position
            SPM.xBF.T0 = 1; %shouldn't need it - no offset
            SPM.xBF.UNITS = 'secs';
            SPM.xBF.Volterra = volt;
            SPM.xY.RT = 1/fs; %xBF.dt;
            
            U.name = {'one'};
            U.dt = xBF.dt;
            U.ons = 1;
            U.dur = 0;
            U.P.name = 'none';
            U.P.h    = 0;
            SPM.Sess(s).U = U;
            
            %copied code from spm_get_ons(SPM,s):
            % create stimulus functions (32 bin offset)
            %======================================================================
            U = spm_get_ons(SPM,s);
            V = SPM.xBF.Volterra;
            %convolve stimuli U with basis functions
            [X,Xn,Fc] = spm_Volterra(U,bf,V);
            try
                X = X((0:(ns - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
            end
            
            %             % and orthogonalise (within trial type)
            %             %--------------------------------------
            %             for i = 1:length(Fc)
            %                 X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
            %             end
            bf_norm = 1/max(X(:,1)); %For more precision, could additionally
            %filter X before calculating bf_norm, but will be a small
            %effect
        end
        
        if tp
            %Block paradigm type
            %interval length
            il = floor(ns/(tn+1));
            %stimuli array (in seconds)
            as = (il*(1:tn)+1)*xBF.dt;
            count = tn;
        else
            try
                mtstream = RandStream('mt19937ar','Seed',rseed);
            catch
                mtstream = RandStream('mt19937ar','Seed',0);
            end
            RandStream.setDefaultStream(mtstream);
            
            %Event paradigm - exponential process with possibly 2 phases with
            %different decay parameters
            count = 0; %total spikes
            tlen = 0; %cumulative time in arbitrary units
            ds = []; %spike positions in arbitrary units
            flen = ns/fs; %length of file in seconds
            if tp2
                %Two types of spikes
                fcount = 0; %frequent spikes
                %while not all desired spikes have been generated
                while count < tn %subtract one from tn because will add one
                    %infrequent spike at the end
                    %alternate between infrequent and frequent spikes
                    p1 = poissrnd(tssg);
                    for i1=1:p1
                        e1 = exprnd(tss);
                        tlen = tlen + e1;
                        if count > tn-1, break; end
                        %check if exceed file length only if not rescaling
                        if ~tsc && tlen > flen, break; end
                        ds = [ds tlen];
                        count = count + 1;
                    end
                    %to escape while count < tn-1
                    if count > tn-1, break; end
                    if ~tsc && tlen > flen, break; end
                    p2 = poissrnd(tsfg);
                    for i1=1:p2
                        e2 = unifrnd(tfs(1),tfs(2));
                        tlen = tlen + e2;
                        if count > tn-1, break; end
                        if ~tsc && tlen > flen, break; end
                        ds = [ds tlen];
                        count = count + 1;
                        fcount = fcount + 1;
                    end
                    %to escape while count < tn-1
                    if count > tn-1, break; end
                    if ~tsc && tlen > flen, break; end
                end %end while
            else
                %One type of spikes
                while count < tn
                    e1 = exprnd(tss);
                    tlen = tlen + e1;
                    ds = [ds tlen];
                    count = count + 1;
                    if ~tsc && tlen > flen, break; end
                end
            end
            %add a final interval after the last spike from the slow
            %distribution
            e1 = exprnd(tss);
            tlen = tlen + e1;
            if tsc
                %rescale to fit size of file and convert to seconds
                as = ds*size(d,2)/(tlen*fs);
            else
                %ds is already in seconds
                as = ds;
            end
        end
        if tsc && (count > tn || count < tn)
            disp('problem with number of spikes');
        end
        
        %fill SPM's xBF structure from spm_get_bf %fill SPM's U structure from spm_get_ons
        SPM = [];
        %SPM's session
        s = 1;
        SPM.nscan(s) = ns;
        SPM.xBF = xBF;
        SPM.xBF.T = xBF.T; %get more precision on onsets position
        SPM.xBF.T0 = 1; %shouldn't need it - no offset
        SPM.xBF.UNITS = 'secs';
        SPM.xBF.Volterra = volt;
        SPM.xY.RT = 1/fs; %xBF.dt;
        
        U.name = {tname};
        U.dt = xBF.dt;
        U.ons = as;
        U.dur = 0;
        U.P.name = 'none';
        U.P.h    = 0;
        SPM.Sess(s).U = U;
        
        %copied code from spm_get_ons(SPM,s):
        % create stimulus functions (32 bin offset)
        %======================================================================
        U = spm_get_ons(SPM,s);
        V = SPM.xBF.Volterra;
        %convolve stimuli U with basis functions
        [X,Xn,Fc] = spm_Volterra(U,bf,V);
        try
            X = X((0:(ns - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
        end
        
                % and orthogonalise (within trial type)
                %--------------------------------------
                for i = 1:length(Fc) %This does not orthogonalize the Volterras
                    X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
                end
        
        %X = X(33:end,:);
        if remove_X_mean
            X = X-repmat(mean(X,1),[size(X,1),1]);
        end
        NIRS.Dt.fir.filter_X = filter_X;
        if filter_X
            if use_gaussian_LPF
                switch use_wavelet_filter
                    case 1
                        HParam.type = 'Wavelet-MDL';
                        HParam.M = 4;
                    case 0
                        HParam.type = 'none'; %'Wavelet-MDL';
                end                
                LParam.FWHM = LPFGaussianFWHM; %in seconds 0.667;
                LParam.type = 'Gaussian';
                K = struct('HParam', HParam,'row', 1:ns ,'RT', 1/fs,'LParam', LParam);
                K = spm_filter_HPF_LPF_WMDL(K);
                switch HParam.type
                    case 'Wavelet-MDL'
                        K.X = X;
                end
                switch volt
                    case 1
                        fX1 = spm_filter_HPF_LPF_WMDL(K, X(:,1));
                    case 2
                        fX1 = spm_filter_HPF_LPF_WMDL(K, X(:,1)+tb*X(:,2));
                        FX1 = spm_filter_HPF_LPF_WMDL(K, X(:,1));
                        FX2 = spm_filter_HPF_LPF_WMDL(K, X(:,2));
                end
                if use_Butter_HPF
                    %HPF
                    cutoff=HPFButt;
                    FilterOrder=HPFButt_order;
                    Wn=cutoff*2/fs;
                    [fb,fa]=butter(FilterOrder,Wn,'high');
                    switch volt
                    case 1
                        fX1 =  filtfilt(fb,fa,fX1);
                    case 2
                        fX1 =  filtfilt(fb,fa,fX1);
                        FX1 =  filtfilt(fb,fa,FX1);
                        FX2 =  filtfilt(fb,fa,FX2);
                    end
                end
                if use_OLD_wrong_Butter_LPF
                    %LPF
                    cutoff=LFP_butter;
                    FilterOrder=5;
                    Wn=cutoff*2/fs;
                    [fb,fa]=butter(FilterOrder,Wn);
                    switch volt
                    case 1
                        fX1 =  filtfilt(fb,fa,fX1);
                    case 2
                        fX1 =  filtfilt(fb,fa,fX1);
                        FX1 =  filtfilt(fb,fa,FX1);
                        FX2 =  filtfilt(fb,fa,FX2);
                    end
                end
            else
                %LPF
                cutoff=0.667;
                FilterOrder=5;
                Wn=cutoff*2/fs;
                [fb,fa]=butter(FilterOrder,Wn);
                
                switch volt
                    case 1
                        fX1=filtfilt(fb,fa,X(:,1));
                    case 2
                        fX1 = filtfilt(fb,fa,X(:,1)+tb*X(:,2));
                        FX1 = filtfilt(fb,fa,X(:,1));
                        FX2 = filtfilt(fb,fa,X(:,2));
                end
                %HPF
                cutoff=HPFButt;
                FilterOrder=5;
                Wn=cutoff*2/fs;
                [fb,fa]=butter(FilterOrder,Wn,'high');
                fX1 = filtfilt(fb,fa,fX1);
            end            
        else
            switch volt
                case 1
                    fX1 = X(:,1);
                case 2
                    fX1 = X(:,1)+tb*X(:,2);
            end
        end
        
        %calculate power of protocol
        %this will not work if we have subsessions
        %Assumes that first column is the canonical HRF
        %m1 = mean(fX1);
        NIRS.Dt.fir.Ep = std(fX1)^2; %sum((fX1-m1).^2)/ns; %length(fX1);
        %Various options - should all be chosen except calculate_bf_norm
        %and normalize_impulse_area_to_unity
        NIRS.Dt.fir.filter_data = filter_data;
        NIRS.Dt.fir.calculate_bf_norm = calculate_bf_norm;
        NIRS.Dt.fir.GammaOn = GammaOn;
        NIRS.Dt.fir.filter_X = filter_X;
        NIRS.Dt.fir.std_or_power = std_or_power;
        NIRS.Dt.fir.remove_X_mean = remove_X_mean;
        NIRS.Dt.fir.use_gaussian_LPF = use_gaussian_LPF;
        NIRS.Dt.fir.use_Butter_HPF = use_Butter_HPF;
        NIRS.Dt.fir.HPFButt_order = HPFButt_order;
        NIRS.Dt.fir.normalize_impulse_area_to_unity = normalize_impulse_area_to_unity;
        %rescale X by channel:
        for Cidx=1:length(chn)
            %take absolute value of median, or sign of response is
            %unpredictable
            %m = abs(median(dc(Cidx,:)));
            %Take standard deviation instead
            
            %calculate power of baseline
            
            %filter raw data to calculate power between two frequencies
            if filter_data
                %LPF
                if use_gaussian_LPF
                    tdc = spm_filter_HPF_LPF_WMDL(K, dc(Cidx,:)')';    
                    if use_Butter_HPF
                        %HPF
                        cutoff=HPFButt;
                        FilterOrder=HPFButt_order;
                        Wn=cutoff*2/fs;
                        [fb,fa]=butter(FilterOrder,Wn,'high');
                        tdc = filtfilt(fb,fa,tdc);
                    end
                    if use_OLD_wrong_Butter_LPF
                    %LPF
                        cutoff=LFP_butter;
                        FilterOrder=5;
                        Wn=cutoff*2/fs;
                        [fb,fa]=butter(FilterOrder,Wn);
                        tdc = filtfilt(fb,fa,tdc);
                    end
                else
                    %Butterworth - careful: number of degrees of freedom
                    %will be incorrectly calculated!
                    cutoff=0.667; %SPM.xX.lpf_butter_freq; %0.666; %Hz, or 1.5s
                    FilterOrder=5;
                    Wn=cutoff*2/fs;                           % normalised cutoff frequency
                    [fb,fa]=butter(FilterOrder,Wn);            % buterworth filter
                    tdc=filtfilt(fb,fa,dc(Cidx,:));
                    %HPF
                    cutoff=HPFButt; %SPM.xX.hpf_butter_freq; %Hz, or 100s (250s may be optimal)
                    FilterOrder=5;
                    Wn=cutoff*2/fs;
                    [fb,fa]=butter(FilterOrder,Wn,'high');
                    tdc=filtfilt(fb,fa,tdc);
                  
                end
                
            else
                tdc = dc(Cidx,:);
            end
            if Cidx==1
                stdc = tdc;
            end
            NIRS.Dt.fir.Eb(Cidx) = std(tdc)^2; %sum((tdc-m1).^2)/ns;

            if AmpTargetMethod
                %target amplitude
                if length(ta) > 1
                    a2 = ta(Cidx);
                else
                    a2 = ta;
                end
                %Boolean - for testing
                if std_or_power
                    %if filter_data
                        m = NIRS.Dt.fir.Eb(Cidx)^0.5; %std(tdc);
                    %else
                    %    m = std(dc(Cidx,:));
                    %end
                    %a = a2*m/std(fX1);
                    if calculate_bf_norm
                        a = a2*m*bf_norm;
                    else
                        a = a2*m;
                    end
                    NIRS.Dt.fir.SNR(Cidx) = 10*log10(a^2*NIRS.Dt.fir.Ep/NIRS.Dt.fir.Eb(Cidx));
                    %same as 10*log10(a2^2*NIRS.Dt.fir.Ep);
                else
                    m = NIRS.Dt.fir.Eb(Cidx)/NIRS.Dt.fir.Ep;
                    if calculate_bf_norm
                        a = a2*m^(0.5)*bf_norm;
                    else
                        a = a2*m^(0.5);
                    end
                    NIRS.Dt.fir.SNR(Cidx) = 10*log10(a2);
                end
                NIRS.Dt.fir.a(Cidx) = a;
                NIRS.Dt.fir.a2(Cidx) = a2; %/m;
                %NIRS.Dt.fir.a3(Cidx) = a/median(tdc);
                %NIRS.Dt.fir.a4(Cidx) = a2*m/std(fX1);
                %NIRS.Dt.fir.a(Cidx) = a;
                NIRS.Dt.fir.std_or_power = std_or_power;
            else
                %target SNR
                a = (10^(tSNR/10)*NIRS.Dt.fir.Eb(Cidx)/NIRS.Dt.fir.Ep)^(0.5);
                NIRS.Dt.fir.SNR(Cidx) = tSNR;
                NIRS.Dt.fir.a(Cidx) = a;
                m = std(dc(Cidx,:));
                NIRS.Dt.fir.a2(Cidx) = a/m;
                NIRS.Dt.fir.a3(Cidx) = a/median(dc(Cidx,:));
            end
            try NIRS.Dt.fir.AmpTargetMethod = AmpTargetMethod; end
            
            %is only a few percent point-by-point on the stimuli
            %a = 100; for a test...
            %dc(Cidx,:) = zeros(size(dc(Cidx,:)));
            switch volt
                case 1
                    if ~set_baseline_to_zero
                        d(chn(Cidx),:) = dc(Cidx,:)+X'*a;
                    else
                        d(chn(Cidx),:) = X'*a;
                        d_copy(chn(Cidx),:) = zeros(size(d_copy(chn(Cidx),:)));
                    end
                case 2
                    if length(tb)>1
                        b = tb(Cidx); %relative weight of 2nd Volterra regressor to the first
                    else
                        b = tb;
                    end
                    NIRS.Dt.fir.b = b;
                    a = [a b*a];
                    if ~set_baseline_to_zero
                        d(chn(Cidx),:) = dc(Cidx,:)+a*(X');
                    else
                        d(chn(Cidx),:) = a*(X');
                        d_copy(chn(Cidx),:) = zeros(size(d_copy(chn(Cidx),:)));
                    end
            end
        end %end for Cidx
        %[dir1 fil1 ext1] = fileparts(NIRS.Dt.fir.raw.p{ts,:});
        if ~AllChannels
            %channels to keep
            chnk = [tk tk+(1:length(NIRS.Cf.dev.wl)-1)*NC/length(NIRS.Cf.dev.wl)];
            d = d(chnk,:);
        end
        if job.testDupChannels
            if ~AllChannels
                %add a copy of the channels kept, but without adding stimuli
                d = [d; d_copy(chnk,:)];
            else
                d = [d; d_copy];
            end
        end
        if save_all_data %careful, may not work if first index not chosen
            Z = [];
            Z.d = d(1,:); %data with response
            Z.n = d_copy(1,:); %data without response = noise
            Z.U = U; %Onsets
            Z.X = X; %Design matrix
            Z.f = stdc; %filtered data
            Z.fX1 = fX1; %filtered combined regressors
            try
                Z.FX1 = FX1; %filtered first Volterra
                Z.FX2 = FX2; %filtered second Volterra
                save('All_data_for_figure','Z');
            end
        end
        %save
        testp = fullfile(dir1,tname);
        if ~exist(testp,'dir'), mkdir(testp); end
        testpn = fullfile(testp,[fil1 ext1]);
        fwrite_NIR(testpn,d);
        NIRS.Dt.fir.pp(lst+1).p{1} = testpn;
        NIRS.Dt.fir.pp(lst+1).pre = 'addTestStimuli';
        NIRS.Dt.fir.pp(lst+1).job = job;
        try
            NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
            NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
            NIRS.Dt.fir.pp(lst+1).si{f,1} = si;
            NIRS.Dt.fir.pp(lst+1).ei{f,1} = ei;
        catch
        end
        %Generate onset files - to do - actually better to simply store
        %in NIRS.mat the names, onsets, durations structure
        NIRS.Dt.fir.Sess(f).U = U; %only one session for test, and only one stimulus
        %number of frequent spikes and total number of spikes:
        try if tp2
                NIRS.Dt.fir.NfrSpk = fcount;
            end
        catch
            
        end
        NIRS.Dt.fir.NtotSpk = count;
        if job.testDupChannels
            if ~AllChannels
                NIRS.Cf.H.C.N = 2*length(chnk);
                try NIRS.Cf.H.C.n = [NIRS.Cf.H.C.n(chnk) NIRS.Cf.H.C.n(chnk)]; end
                try NIRS.Cf.H.C.id = [NIRS.Cf.H.C.id(:,chnk) NIRS.Cf.H.C.id(:,chnk)]; end
                try NIRS.Cf.H.C.wl = [NIRS.Cf.H.C.wl(chnk) NIRS.Cf.H.C.wl(chnk)]; end
                try NIRS.Cf.H.C.gp = [NIRS.Cf.H.C.gp(chnk) NIRS.Cf.H.C.gp(chnk)]; end
                try NIRS.Cf.H.C.ok = [NIRS.Cf.H.C.ok(chnk) NIRS.Cf.H.C.ok(chnk)]; end
            else
                NIRS.Cf.H.C.N = 2*NIRS.Cf.H.C.N;
                try NIRS.Cf.H.C.n = [NIRS.Cf.H.C.n NIRS.Cf.H.C.n]; end
                try NIRS.Cf.H.C.id = [NIRS.Cf.H.C.id NIRS.Cf.H.C.id]; end
                try NIRS.Cf.H.C.wl = [NIRS.Cf.H.C.wl NIRS.Cf.H.C.wl]; end
                try NIRS.Cf.H.C.gp = [NIRS.Cf.H.C.gp NIRS.Cf.H.C.gp]; end
                try NIRS.Cf.H.C.ok = [NIRS.Cf.H.C.ok NIRS.Cf.H.C.ok]; end
            end
        else
            if ~AllChannels
                NIRS.Cf.H.C.N = length(chnk);
                try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(chnk); end
                try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,chnk); end
                try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(chnk); end
                try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(chnk); end
                try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(chnk); end
            end
        end
        NIRSmat = fullfile(testp,'NIRS.mat');
        save(NIRSmat,'NIRS');
        outNIRSmat = [outNIRSmat; NIRSmat];
        %disp('*************************************************************');
        disp('Recall that NIRS.mat for testing is now in the testing folder');
        %disp('*************************************************************');
    catch exception
        disp(exception);
        disp(['Adding stimuli for testing failed for subject' int2str(Idx)]);
    end
end
out.NIRSmat = outNIRSmat;