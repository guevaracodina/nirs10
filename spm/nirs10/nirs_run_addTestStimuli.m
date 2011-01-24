function out = nirs_run_addTestStimuli(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Here NIRS.mat in the pipeline will be replaced by a new one in a new
%location
outNIRSmat = {};

%Ensure reproducibility of results by resetting the random seed
  
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
        %channels of interest
        chn = tc+(tl-1)*NC/length(NIRS.Cf.dev.wl);
        
        dc = d(chn,:); %we keep d, as we will write dc over d, and save
        ns = size(d,2);
        %frequency
        fs = NIRS.Cf.dev.fs; 
        %fill SPM's xBF structure from spm_get_bf
        xBF.dt = 1/fs; % - time bin length {seconds}
        xBF.name = 'hrf'; %description of basis functions specified            
        xBF = spm_get_bf(xBF);
        bf  = xBF.bf;
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
            flen = ns*xBF.dt; %length of file in seconds
            if tp2
                %Two types of spikes               
                fcount = 0; %frequent spikes
                %while not all desired spikes have been generated
                while count < tn-1 %subtract one from tn because will add one
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
                        if count > tn, break; end
                        if tlen > flen, break; end
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
                while count < tn-1
                    e1 = exprnd(tss);
                    tlen = tlen + e1;
                    ds = [ds tlen];                 
                    count = count + 1;
                    if tlen > flen, break; end
                end               
            end
            %add a final interval after the last spike from the slow
            %distribution
            e1 = exprnd(tss);
            tlen = tlen + e1;
            if tsc 
                %rescale to fit size of file and convert to seconds
                as = ds*size(d,2)*xBF.dt/tlen;
            else
                %ds is already in seconds
                as = ds;
            end
        end
        
        
        %fill SPM's U structure from spm_get_ons
        SPM = [];
        %SPM's session
        s = 1;
        SPM.nscan(s) = ns;
        SPM.xBF = xBF;
        SPM.xBF.T = 1;
        %SPM.xBF.T0 = 1; %shouldn't need it
        SPM.xBF.UNITS = 'secs';
        SPM.xBF.Volterra = volt;
        SPM.xY.RT = xBF.dt;

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
        X = X(33:end,:);
        
        %calculate power of protocole
        switch volt
            case 1 
                %this will not work if we have subsessions
                %Assumes that first column is the canonical HRF
                m1 = mean(X(:,1));
                NIRS.Dt.fir.Ep = sum((X(:,1)-m1).^2)/size(X,1);
            case 2
                %Assumes first column is canonical HRF and second column
                %is second Volterra kernel
                m1 = mean(X(:,1)+tb*X(:,2));
                %this assumes that tb is the same for each channel
                NIRS.Dt.fir.Ep = sum((X(:,1)+tb*X(:,2)-m1).^2)/size(X,1);
        end
        %sum contributions from each regressor with a weight of 1
        %X = sum(X,2);
        %X should be our regressors
        %rescale X by channel:
        for Cidx=1:length(chn)
            %take absolute value of median, or sign of response is
            %unpredictable
            %m = abs(median(dc(Cidx,:)));
            %Take standard deviation instead
            
            %calculate power of baseline
            m1 = mean(dc(Cidx,:));
            NIRS.Dt.fir.Eb(Cidx) = sum((dc(Cidx,:)-m1).^2)/ns;
            
            
            if AmpTargetMethod
                %target amplitude
                if length(ta) > 1
                    a2 = ta(Cidx);
                else
                    a2 = ta;
                end
                
                %Boolean - for testing
                std_or_power = 0;
                if std_or_power
                    m = std(dc(Cidx,:)); 
                    a = a2*m;
                    NIRS.Dt.fir.SNR(Cidx) = 10*log10((a2*m)^2*NIRS.Dt.fir.Ep/NIRS.Dt.fir.Eb(Cidx));
                else
                    m = NIRS.Dt.fir.Eb(Cidx)/NIRS.Dt.fir.Ep;
                    a = a2*m^(0.5);
                    NIRS.Dt.fir.SNR(Cidx) = 10*log10(a2); 
                end
                                
                                                           
                NIRS.Dt.fir.a(Cidx) = a;
                NIRS.Dt.fir.std_or_power = std_or_power;
            else
                %target SNR
                a = (10^(tSNR/10)*NIRS.Dt.fir.Eb(Cidx)/NIRS.Dt.fir.Ep)^(0.5);
                NIRS.Dt.fir.SNR(Cidx) = tSNR;
                NIRS.Dt.fir.a(Cidx) = a;
                
            end
                               
            %is only a few percent point-by-point on the stimuli
            switch volt
                case 1                    
                    d(chn(Cidx),:) = dc(Cidx,:)+X'*a; 
                case 2                    
                    if length(tb)>1
                        b = tb(Cidx); %relative weight of 2nd Volterra regressor to the first
                    else
                        b = tb;
                    end
                    a = [a b*a]; 
                    d(chn(Cidx),:) = dc(Cidx,:)+a*(X'); 
            end               
        end %end for Cidx
        %[dir1 fil1 ext1] = fileparts(NIRS.Dt.fir.raw.p{ts,:});
        if ~AllChannels
            %channels to keep 
            chnk = [tk tk+(1:length(NIRS.Cf.dev.wl)-1)*NC/length(NIRS.Cf.dev.wl)];
            d = d(chnk,:);
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
        NIRS.Dt.fir.Sess.U = U; %only one session for test, and only one stimulus
        %number of frequent spikes and total number of spikes:
        try if tp2
            NIRS.Dt.fir.NfrSpk = fcount;
            end
        catch
            
        end
        NIRS.Dt.fir.NtotSpk = count;
        if ~AllChannels
            NIRS.Cf.H.C.N = length(chnk);
            try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(chnk); end
            try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,chnk); end
            try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(chnk); end 
            try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(chnk); end
            try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(chnk); end 
        end
        NIRSmat = fullfile(testp,'NIRS.mat');
        save(NIRSmat,'NIRS');
        outNIRSmat = [outNIRSmat; NIRSmat];
        %disp('*************************************************************');
        disp('Recall that NIRS.mat for testing is now in the testing folder');
        %disp('*************************************************************');
    catch
        disp(['Adding stimuli for testing failed for subject' int2str(Idx)]);
    end
end
out.NIRSmat = outNIRSmat;