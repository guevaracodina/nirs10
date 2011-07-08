function out = nirs_run_ROCtest(job)
dirROC = pwd;
%Hard-coded
try 
    chan_len = job.ROCnumCh;
catch
    chan_len = 40;
end
all_channels = 1;
half_chan_len = chan_len/2;
quarter_chan_len = chan_len/4;
generate_ROC_curves = 1; %Boolean
%specify some of the directory structure
try 
    dir_dataSPM = job.dir_dataSPM;
catch
    dir_dataSPM = 'dataSPM';
end
%dir_stat = 'StatV';
%alpha_unc = 0.05; %uncorrected threshold
exp_th_list = [0 10 19 22 24 28 31 35 39 42 47 66 72 80 90 100 200];
%exp_th_list = [0 10 19 22 24 28 31 35 39 42 47 52 59 66 69 72 74 76 78 80 85 90 95 100 105 110 120 130 150 180 200];
%options for the plots
try 
    byIter = job.byIter;
catch
    byIter = 0; %Boolean
end
try
    compute_OR = job.compute_OR;
catch
    compute_OR = 0; %Boolean
end
try
    compute_LU = job.compute_LU;
catch
    compute_LU = 0; %Boolean
end
try
    runFtest = job.runFtest;
catch
    runFtest = 0; %Boolean
end
%Here specifify whether run Jidx for Subject Idx will be later looked at
%with a positive (1) or negative (0) t-test for the 2nd Volterra
%this value will be in time assigned to Volt2_positive_ttest

IN = job.ROCiternum;
try 
    RunGLMorFigures = job.RunGLMorFigures;
catch
    RunGLMorFigures = 3;
end
switch RunGLMorFigures
    case 1
        run_GLM = 1;
        run_ROC = 0;
    case 2
        run_GLM = 0;
        run_ROC = 1;
    case 3
        run_GLM = 1;
        run_ROC = 1;
end
nSubj = size(job.NIRSmat,1);
nJob = size(job.ROCLoopJob,1);
Volt2_positive_ttest = false;

try
    v2a = job.Volt2;
    if size(v2a,1)*size(v2a,2) == 1
        for i1=1:size(v2a,1)
            for j1=1:size(v2a,2)
                Volt2{i1,j1} = v2a;
            end
        end
    else
        for i1=1:size(v2a,1)
            for j1=1:size(v2a,2)
                %careful - transposed here... usually only one subject
                Volt2{j1,i1} = v2a(i1,j1);
            end
        end 
    end       
catch
    Volt2{1,1} = 0;
    Volt2{2,1} = 0;
    Volt2{3,1} = 0;
    Volt2{4,1} = 0;
    Volt2{5,1} = 0;
    Volt2{6,1} = 0;
    Volt2{7,1} = 0;
    Volt2{8,1} = 0;
    Volt2{9,1} = 0;
    Volt2{10,1} = 0;
end

%Loop over subjects
for Idx=1:nSubj
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir0,dummy,dummy2] = fileparts(job.NIRSmat{Idx,1});
        dir1 = [dir0 filesep dir_dataSPM];
        
        %loop over jobs
        for Jidx=1:nJob
            %loop over the specified job
            for ROCiter=1:IN
                clear LoopJob
                %Load/read job
                LoopJob = load(job.ROCLoopJob{Jidx});
                %replace whatever is specified in the LoopJob NIRSmat by
                %subject Idx
                LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets.addTestStimuli.NIRSmat = {job.NIRSmat{Idx,1}};
                %set the seed
                try %Frequent spikes
                    LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets.addTestStimuli. ...
                        testPType.testEP.FrequentSpikes.testSeed2; %check if testSeed2 exists
                    LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets.addTestStimuli. ...
                        testPType.testEP.FrequentSpikes.testSeed2 = ROCiter;
                catch %No Frequent Spikes
                    LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets.addTestStimuli. ...
                        testPType.testEP.NoFrequentSpikes.testSeed1 = ROCiter;               
                end
                %build folder name
                testName = LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets. ...
                        addTestStimuli.testStimulusName;
                %N_IN = ceil(log10(IN));
                %if ROCiter < 10                    
                testFullName = [testName int2str(ROCiter)];
                LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets. ...
                    addTestStimuli.testStimulusName = testFullName;

                if run_GLM
                    %run the job
                    spm_jobman('run',LoopJob.matlabbatch);
                    %delete files
                    if job.ROCDeleteLarge == 1
                        testDir = fullfile(dir1,testFullName);
                        largefiles = spm_select('FPList',testDir,'.nir');
                        try 
                            for i1=1:size(largefiles,1)
                                delete(deblank(largefiles(i1,:)));
                            end
                        end
                    end
                    disp(['Job: ' int2str(Jidx) ' -- Iter: ' int2str(ROCiter)]);
                end
            end

            if run_ROC
                %preloading for array-size assignments
                clear LoopJob
                %Load/read job
                LoopJob = load(job.ROCLoopJob{Jidx});
                %build folder name
                testName = LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets. ...
                        addTestStimuli.testStimulusName;
                testFullName = [testName int2str(1)];
                dir_stat = LoopJob.matlabbatch{2}.spm.tools.nirs10.model_specify. ...
                    wls_bglm_specify.dir1;
                dir_spm = [dir1 filesep testFullName filesep dir_stat];
                SPM = [];
                load(fullfile(dir_spm,'SPM.mat'));                
                %For F-test
                if runFtest
                try
                    jobf = [];
                    jobf.consess{1}.fcon.name = 'F'; 
                    jobf.consess{1}.fcon.convec = {eye(2)}; 
                    jobf.consess{1}.fcon.sessrep = 'none'; 
                    SPM.xCon = [];
                    SPM.xX = SPM.xXn{1};
                    SPM = nirs_spm_run_con(jobf,SPM);
                    xCon = SPM.xCon; %F-contrast
                end
                end
                try 
                    SPM.xX.Bvar;
                    %ROC
                    sz_beta = size(SPM.beta);
                    sz_Bvar = size(SPM.xX.Bvar);
                    try
                        if SPM.GenerateHbT
                            sz_beta(2) = chan_len; %quick fix for bug due to
                            %introducing HbT
                            sz_ResSS(2) = chan_len; %quick fix for bug due to introducing HbT
                        end
                    end
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bvar   = zeros(IN,sz_Bvar(1),sz_Bvar(2));
                    tF   = zeros(IN,sz_beta(2)); %F-test
                    RE = zeros(IN,3,sz_beta(2)/2);
                    SNR = zeros(IN,sz_beta(2)/2);
                    CORR = zeros(IN,1); %correlation between 1st and 2nd Volterra
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        
                        try
                            if SPM.GenerateHbT
                                beta(ROCiter,:,:) = SPM.beta(:,1:chan_len); %quick fix for bug due to introducing HbT
                                Bvar(ROCiter,:,:) = full(SPM.xX.Bvar(:,1:chan_len)); %quick fix for bug due to introducing HbT
                            else
                                beta(ROCiter,:,:) = SPM.beta; %(:,1:chan_len); %quick fix for bug due to introducing HbT
                                Bvar(ROCiter,:,:) = full(SPM.xX.Bvar) ; %(:,1:chan_len)); %quick fix for bug due to introducing HbT
                            end
                        catch
                            beta(ROCiter,:,:) = SPM.beta; %(:,1:chan_len); %quick fix for bug due to introducing HbT
                            Bvar(ROCiter,:,:) = full(SPM.xX.Bvar) ; %(:,1:chan_len)); %quick fix for bug due to introducing HbT
                        end
                        NIRS = [];
                        load(fullfile(dir_spm,'NIRS.mat'));
                        %1st Volterra kernel
                        RE(ROCiter,1,:) = SPM.beta(1,1:(sz_beta(2)/2))./ NIRS.Dt.fir.a -1;
                        %2nd Volterra kernel
                        try
                            RE(ROCiter,2,:) = SPM.beta(2,1:(sz_beta(2)/2))./ (NIRS.Dt.fir.b*NIRS.Dt.fir.a) -1;
                        end
                        SNR(ROCiter,:) = NIRS.Dt.fir.SNR;
                        CORR(ROCiter,1) = corr(SPM.xX.X(:,1),SPM.xX.X(:,2));
                    end

                    %get t stat for each iteration
                    t = zeros(IN,2,sz_beta(2));
                    %loop over 1st and 2nd Volterra regressors
                    for r=1:2
                        for i=1:IN
                            t(i,r,:) = squeeze(beta(i,r,:)) ./sqrt(squeeze(Bvar(i,r,:)));
                        end
                    end
                    %very approximate value for erdf - assume enough
                    %degrees of freedom to approximate t-stat by Z-score
                    erdf   = 1000*ones(IN,1);                    
                catch
                    %ROC
                    sz_beta = size(SPM.xXn{1}.beta);
                    
                    sz_Bcov = size(SPM.xXn{1}.Bcov);
                    sz_ResSS = size(SPM.xXn{1}.ResSS);
                    try
                        if SPM.GenerateHbT
                            sz_beta(2) = chan_len; %quick fix for bug due to
                            %introducing HbT
                            sz_ResSS(2) = chan_len; %quick fix for bug due to introducing HbT
                        end
                    end
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bcov   = zeros(IN,sz_Bcov(1),sz_Bcov(1));
                    tF     = zeros(IN,sz_beta(2)); %F-test
                    ResSS  = zeros(IN,sz_ResSS(2));
                    trRV   = zeros(IN,1);
                    trRVRV = zeros(IN,1);
                    erdf   = zeros(IN,1);
                    RE = zeros(IN,3,sz_beta(2)/2);
                    SNR = zeros(IN,sz_beta(2)/2);
                    CORR = zeros(IN,1); %correlation between 1st and 2nd Volterra
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        %spm_DesRep('DesRepUI',SPM);
                        Bcov(ROCiter,:,:)   = SPM.xXn{1}.Bcov;
                        try
                            if SPM.GenerateHbT
                                beta(ROCiter,:,:)   = SPM.beta(:,1:chan_len); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.ResSS(1,1:chan_len); %quick fix for bug due to introducing HbT
                            else
                                beta(ROCiter,:,:)   = SPM.xXn{1}.beta; %(:,1:chan_len); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.xXn{1}.ResSS; %(1,1:chan_len); %quick fix for bug due to introducing HbT
                            end
                        catch
                            if SPM.GenerateHbT
                                beta(ROCiter,:,:)   = SPM.xXn{1}.beta(:,1:chan_len); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.xXn{1}.ResSS(:,1:chan_len);
                            else
                                beta(ROCiter,:,:)   = SPM.xXn{1}.beta; %(:,1:chan_len); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.xXn{1}.ResSS; %(1,1:chan_len); %quick fix for bug due to introducing HbT
                            end
                        end
                        CORR(ROCiter,1) = corr(SPM.xXn{1}.X(:,1),SPM.xXn{1}.X(:,2));
                        trRV(ROCiter)   = SPM.xXn{1}.trRV;
                        trRVRV(ROCiter) = SPM.xXn{1}.trRVRV;
                        erdf(ROCiter)   = trRV(ROCiter)^2/trRVRV(ROCiter);
                        if isnan(erdf(ROCiter))
                            erdf(ROCiter) = 491;
                            trRV(ROCiter) = 343;
                            trRVRV(ROCiter) = 239;
                        end
                        NIRS = [];
                        load(fullfile(dir_spm,'NIRS.mat'));
                        %1st Volterra kernel
                        RE(ROCiter,1,:) = SPM.xXn{1}.beta(1,1:(sz_beta(2)/2))./ NIRS.Dt.fir.a; % -1;
                        %2nd Volterra kernel
                        try
                            RE(ROCiter,2,:) = SPM.xXn{1}.beta(2,1:(sz_beta(2)/2))./ (NIRS.Dt.fir.b*NIRS.Dt.fir.a);% -1;
                            %for 3rd component, calculate the ratio of 2nd
                            %to 1st Volterra estimated amplitudes, less simulated 
                            RE(ROCiter,3,:) = SPM.xXn{1}.beta(2,1:(sz_beta(2)/2))./ ...
                                SPM.xXn{1}.beta(1,1:(sz_beta(2)/2));%  - NIRS.Dt.fir.b;
                        end  
                        SNR(ROCiter,:) = NIRS.Dt.fir.SNR;
                        
                        %For F-test
                        if runFtest
                        try 
                            %hF = spm_FcUtil('Hsqr',xCon(1), SPM.xX.X);
                            hF = spm_FcUtil('Hsqr',xCon(1), SPM.xXn{1}.xKXs.X);
                            switch SPM.xXn{1}.K.LParam.type
                                case {'hrf', 'Gaussian'}
                                    S = SPM.xXn{1}.K.KL;
                                case 'none'
                                    S = speye(nScan);
                            end
                            switch SPM.xXn{1}.K.HParam.type 
                                case 'DCT'
                                     S = S - SPM.xXn{1}.K.X0 * (SPM.xXn{1}.K.X0' * S);                                    
                                    %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
                            end
                            %Calculate modified nubmer of degrees of freedom due to filtering
                            %for the sum of squares difference between the full and reduced models
                            trRV2 = approx_trRV(SPM.xXn{1}.xKXs.X,SPM.xXn{1}.pKX,S,xCon(1).c);
                            %xCon(1).eidf_mod = trRV2^2/trRVRV2; %This is
                            %identical to eidf!
                            xCon(1).eidf_mod = xCon(1).eidf;
                            if SPM.GenerateHbT
                                %tF1 = zeros(1,chan_len);
                                for chn=1:chan_len
                                    %tF1 = (squeeze(beta(ROCiter,:,chn))*hF')*(hF*squeeze(beta(ROCiter,:,chn))');%-ResSS(ROCiter,chn); 
                                    %tF(ROCiter,chn) = tF1/(ResSS(ROCiter,chn)*xCon(1).eidf_mod/erdf(ROCiter));
                                    tF1 = (squeeze(beta(ROCiter,:,chn))*hF')*(hF*squeeze(beta(ROCiter,:,chn))')/trRV2;%-ResSS(ROCiter,chn); 
                                    tF(ROCiter,chn) = tF1/(ResSS(ROCiter,chn)/trRV(ROCiter)); %xCon(1).eidf=2                                   
                                end
                            else
                                %tF1 = zeros(1,sz_beta(2));
                                for chn=1:sz_beta(2)
                                    tF1 = (squeeze(beta(ROCiter,:,chn))*hF')*(hF*squeeze(beta(ROCiter,:,chn))')/trRV2;%-ResSS(ROCiter,chn); 
                                    tF(ROCiter,chn) = tF1/(ResSS(ROCiter,chn)/trRV(ROCiter)); %xCon(1).eidf=2
                                end
                            end
                        end    
                        end
                    end

                    %get t stat for each iteration
                    t = zeros(IN,2,sz_beta(2));
                    %loop over 1st and 2nd Volterra regressors
                    for r=1:2
                        for i=1:IN
                            t(i,r,:) = squeeze(beta(i,r,:)) ./sqrt(squeeze(ResSS(i,:)'*Bcov(i,r,r)/trRV(i)));                        
                        end
                    end
                end
                %Store
                T{Jidx,Idx}.t = t;
                T{Jidx,Idx}.b = beta;
                T{Jidx,Idx}.RE = RE;
                T{Jidx,Idx}.SNR = SNR;
                T{Jidx,Idx}.CORR = CORR;
                T{Jidx,Idx}.a2 = NIRS.Dt.fir.a2;
                t1 = squeeze(T{Jidx,Idx}.t(:,1,:));
                t2 = squeeze(T{Jidx,Idx}.t(:,2,:));
                %Add bonferroni and choice of t-stat value.                             
                if generate_ROC_curves
                    %Enormous list of all the variables of interest, by run
                    %and subject
                    
                    %need to loop over thresholds
                    %1st Volterra
                    TPb1{Jidx,Idx} = [];
                    TPu1{Jidx,Idx} = [];
                    FPb1{Jidx,Idx} = [];
                    FPu1{Jidx,Idx} = [];
                    %2nd Volterra
                    TPb2{Jidx,Idx} = [];
                    TPu2{Jidx,Idx} = [];
                    FPb2{Jidx,Idx} = [];
                    FPu2{Jidx,Idx} = [];
                    %for 25th and 75th percentiles (lower and upper) 
                    if compute_LU
                        TPu1u{Jidx,Idx} = [];
                        FPu1u{Jidx,Idx} = [];
                        TPu1l{Jidx,Idx} = [];
                        FPu1l{Jidx,Idx} = [];
                        TPu2u{Jidx,Idx} = [];
                        FPu2u{Jidx,Idx} = [];
                        TPu2l{Jidx,Idx} = [];
                        FPu2l{Jidx,Idx} = [];
                    end
                    %for HbO and HbR separately
                    if compute_OR
                        TPu1O{Jidx,Idx} = [];
                        TPu1R{Jidx,Idx} = [];
                        TPu2O{Jidx,Idx} = [];
                        TPu2R{Jidx,Idx} = [];
                        FPu1O{Jidx,Idx} = [];
                        FPu1R{Jidx,Idx} = [];
                        FPu2O{Jidx,Idx} = [];
                        FPu2R{Jidx,Idx} = [];
                        if compute_LU
                            TPu1Ou{Jidx,Idx} = [];
                            TPu1Ol{Jidx,Idx} = [];
                            TPu2Ou{Jidx,Idx} = [];
                            TPu2Ol{Jidx,Idx} = [];
                            TPu1Ru{Jidx,Idx} = [];
                            TPu1Rl{Jidx,Idx} = [];
                            TPu2Ru{Jidx,Idx} = [];
                            TPu2Rl{Jidx,Idx} = [];
                            FPu1Ou{Jidx,Idx} = [];
                            FPu1Ol{Jidx,Idx} = [];
                            FPu2Ou{Jidx,Idx} = [];
                            FPu2Ol{Jidx,Idx} = [];
                            FPu1Ru{Jidx,Idx} = [];
                            FPu1Rl{Jidx,Idx} = [];
                            FPu2Ru{Jidx,Idx} = [];
                            FPu2Rl{Jidx,Idx} = [];
                        end
                    end
                    if all_channels
                        %1st Volterra
                        TPu1A{Jidx,Idx} = [];
                        FPu1A{Jidx,Idx} = [];
                        %2nd Volterra
                        TPu2A{Jidx,Idx} = [];
                        FPu2A{Jidx,Idx} = [];
                    end
                    %
                    try
                        Volt2_positive_ttest = Volt2{Jidx,Idx};
                    end
                    %Volt2_positive_ttest = false;
                    
                    for exp_alpha_unc = exp_th_list % [2:10 12:2:20 25]
                        alpha_unc = 10.^(-exp_alpha_unc/15);
                        %1st Volterra
                        TPn = half_chan_len; %NcTP*IN; %Number of data points for true positives and false negatives
                        alpha_bonf_TPn = alpha_unc/TPn;
                        [u v] = count_TP_FP(IN,[1:half_chan_len],...
                            t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                        TPb1{Jidx,Idx} = [TPb1{Jidx,Idx} median(u)];
                        TPu1{Jidx,Idx} = [TPu1{Jidx,Idx} median(v)];
                        if compute_LU
                            TPu1u{Jidx,Idx} = [TPu1u{Jidx,Idx} prctile(v,75)-median(v)];
                            TPu1l{Jidx,Idx} = [TPu1l{Jidx,Idx} median(v)-prctile(v,25)];
                        end
                        if all_channels
                            %1st Volterra
                            TPu1A{Jidx,Idx} = [TPu1A{Jidx,Idx} v];
                        end
                        %2nd Volterra
                        [u v] = count_TP_FP(IN,[1:half_chan_len],...
                            t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                        TPb2{Jidx,Idx} = [TPb2{Jidx,Idx} median(u)];
                        TPu2{Jidx,Idx} = [TPu2{Jidx,Idx} median(v)];
                        if compute_LU
                            TPu2u{Jidx,Idx} = [TPu2u{Jidx,Idx} prctile(v,75)-median(v)];
                            TPu2l{Jidx,Idx} = [TPu2l{Jidx,Idx} median(v)-prctile(v,25)];
                        end
                        if all_channels
                            %2nd Volterra
                            TPu2A{Jidx,Idx} = [TPu2A{Jidx,Idx} v];
                        end
                        %False Positive
                        %1st Volterra
                        FPn = half_chan_len;
                        alpha_bonf_FPn = alpha_unc/FPn;
                        [u v ] = count_TP_FP(IN,[(half_chan_len+1):chan_len],...
                            t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPb1{Jidx,Idx} = [FPb1{Jidx,Idx} median(u)];
                        FPu1{Jidx,Idx} = [FPu1{Jidx,Idx} median(v)];
                        if compute_LU
                            FPu1u{Jidx,Idx} = [FPu1u{Jidx,Idx} prctile(v,75)-median(v)];
                            FPu1l{Jidx,Idx} = [FPu1l{Jidx,Idx} median(v)-prctile(v,25)];
                        end
                        if all_channels
                            %1st Volterra
                            FPu1A{Jidx,Idx} = [FPu1A{Jidx,Idx} v];
                        end
                        %2nd Volterra
                        [u v ] = count_TP_FP(IN,[(half_chan_len+1):chan_len],...
                            t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPb2{Jidx,Idx} = [FPb2{Jidx,Idx} median(u)];
                        FPu2{Jidx,Idx} = [FPu2{Jidx,Idx} median(v)];
                        if compute_LU
                            FPu2u{Jidx,Idx} = [FPu2u{Jidx,Idx} prctile(v,75)-median(v)];
                            FPu2l{Jidx,Idx} = [FPu2l{Jidx,Idx} median(v)-prctile(v,25)];
                        end
                        if all_channels
                            %1st Volterra
                            FPu2A{Jidx,Idx} = [FPu2A{Jidx,Idx} v];
                        end
                        if compute_OR
                            %For HbO and HbR separately
                            %1st Volterra
                            TPn = quarter_chan_len; %NcTP*IN; %Number of data points for true positives and false negatives
                            alpha_bonf_TPn = alpha_unc/TPn;
                            %HbO
                            [u v] = count_TP_FP(IN,[1:quarter_chan_len],...
                                t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                            TPu1O{Jidx,Idx} = [TPu1O{Jidx,Idx} median(v)];     
                            if compute_LU
                                TPu1Ou{Jidx,Idx} = [TPu1Ou{Jidx,Idx} prctile(v,75)-median(v)];
                                TPu1Ol{Jidx,Idx} = [TPu1Ol{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %2nd Volterra
                            [u v] = count_TP_FP(IN,[1:quarter_chan_len],...
                                t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                            TPu2O{Jidx,Idx} = [TPu2O{Jidx,Idx} median(v)];                      
                            if compute_LU
                                TPu2Ou{Jidx,Idx} = [TPu2Ou{Jidx,Idx} prctile(v,75)-median(v)];
                                TPu2Ol{Jidx,Idx} = [TPu2Ol{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %False Positive
                            %1st Volterra
                            FPn = quarter_chan_len;
                            alpha_bonf_FPn = alpha_unc/FPn;
                            [u v] = count_TP_FP(IN,[(half_chan_len+1):(half_chan_len+quarter_chan_len)],...
                                t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                            FPu1O{Jidx,Idx} = [FPu1O{Jidx,Idx} median(v)];
                            if compute_LU
                                FPu1Ou{Jidx,Idx} = [FPu1Ou{Jidx,Idx} prctile(v,75)-median(v)];
                                FPu1Ol{Jidx,Idx} = [FPu1Ol{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %2nd Volterra
                            [u v ] = count_TP_FP(IN,[(half_chan_len+1):(half_chan_len+quarter_chan_len)],...
                                t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                            FPu2O{Jidx,Idx} = [FPu2O{Jidx,Idx} median(v)];
                            if compute_LU
                                FPu2Ou{Jidx,Idx} = [FPu2Ou{Jidx,Idx} prctile(v,75)-median(v)];
                                FPu2Ol{Jidx,Idx} = [FPu2Ol{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %HbR
                            [u v] = count_TP_FP(IN,[(quarter_chan_len+1):half_chan_len],...
                                t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                            TPu1R{Jidx,Idx} = [TPu1R{Jidx,Idx} median(v)];                        
                            if compute_LU
                                TPu1Ru{Jidx,Idx} = [TPu1Ru{Jidx,Idx} prctile(v,75)-median(v)];
                                TPu1Rl{Jidx,Idx} = [TPu1Rl{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %2nd Volterra
                            [u v] = count_TP_FP(IN,[(quarter_chan_len+1):half_chan_len],...
                                t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                            TPu2R{Jidx,Idx} = [TPu2R{Jidx,Idx} median(v)];
                            if compute_LU
                                TPu2Ru{Jidx,Idx} = [TPu2Ru{Jidx,Idx} prctile(v,75)-median(v)];
                                TPu2Rl{Jidx,Idx} = [TPu2Rl{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %False Positive
                            %1st Volterra
                            FPn = quarter_chan_len;
                            alpha_bonf_FPn = alpha_unc/FPn;
                            [u v] = count_TP_FP(IN,[(half_chan_len+quarter_chan_len+1):chan_len],...
                                t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                            FPu1R{Jidx,Idx} = [FPu1R{Jidx,Idx} median(v)];
                            if compute_LU
                                FPu1Ru{Jidx,Idx} = [FPu1Ru{Jidx,Idx} prctile(v,75)-median(v)];
                                FPu1Rl{Jidx,Idx} = [FPu1Rl{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                            %2nd Volterra
                            [u v ] = count_TP_FP(IN,[(half_chan_len+quarter_chan_len+1):chan_len],...
                                t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                            FPu2R{Jidx,Idx} = [FPu2R{Jidx,Idx} median(v)];
                            if compute_LU
                                FPu2Ru{Jidx,Idx} = [FPu2Ru{Jidx,Idx} prctile(v,75)-median(v)];
                                FPu2Rl{Jidx,Idx} = [FPu2Rl{Jidx,Idx} median(v)-prctile(v,25)];
                            end
                        end
                    end
%                 else
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %Sensitivity = true positives
%                    
%                     TPn = half_chan_len; %NcTP*IN; %Number of data points for true positives and false negatives
%                     alpha_bonf_TPn = alpha_unc/TPn;
% 
%                     %[TPb{Jidx,Idx} TPu{Jidx,Idx}] = count_TP_FP(IN,[1:2],t,alpha_bonf_TPn,alpha_unc,erdf,true,true); 
%                     %1st Volterra
%                     [TPb1{Jidx,Idx} TPu1{Jidx,Idx}] = count_TP_FP(IN,[1:half_chan_len],...
%                         t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
%                     %2nd Volterra
%                     [TPb2{Jidx,Idx} TPu2{Jidx,Idx}] = count_TP_FP(IN,[1:half_chan_len],...
%                         t2,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
%                     
%                     %Specificity = 1 - false positives
%                     FPn = half_chan_len; %(sz_beta(2)-NcTP)*IN; %Number of data points for false positives 
%                     %and true negatives
%                     alpha_bonf_FPn = alpha_unc/FPn;
%                     %1st Volterra
%                     [FPb1{Jidx,Idx} FPu1{Jidx,Idx}] = count_TP_FP(IN,[(half_chan_len+1):chan_len],...
%                         t1,alpha_bonf_FPn,alpha_unc,erdf,false,true,byIter); 
%                     %2nd Volterra
%                     [FPb2{Jidx,Idx} FPu2{Jidx,Idx}] = count_TP_FP(IN,[(half_chan_len+1):chan_len],...
%                         t2,alpha_bonf_FPn,alpha_unc,erdf,false,true,byIter); 
% 
%                     %tFPb = FPb; tFPu = FPu; tTPb = TPb; tTPu = TPu; 
%                     %FPb = FPb{Jidx,:}; tFPu = FPu{Jidx,:}; tTPb = TPb{Jidx,:}; tTPu = TPu{Jidx,:};                  
%                     %save(fullfile(dirROC,['ROC' testName int2str(Jidx) '.mat']),...
%                     %    'TPb','FPb','TPu','FPu');
%                     %FPb = tFPb; FPu = tFPu; TPb = tTPb; TPu = tTPu;
                end               
                %to display a design matrix:
                %spm_DesRep('DesRepUI',SPM)                              
            end %end if run_ROC           
        end %end for Jidx
        %if run_ROC
        %    save(fullfile(dirROC,['ROC' testName '.mat']),'TPb','FPb','TPu','FPu');
        %end                        
    catch exception
        disp(['Could not perform ROC iterations subject ' int2str(Idx)]);
        disp(exception.identifier);
    end   
end %end for Idx

if run_ROC
if compute_OR
    if compute_LU
        save(fullfile(dirROC,'ROC.mat'),'T','TPb1','FPb1','TPu1','FPu1','TPb2','FPb2','TPu2','FPu2',...
            'TPu1O','FPu1O','TPu2O','FPu2O','TPu1R','FPu1R','TPu2R','FPu2R',...
            'TPu1Ou','FPu1Ou','TPu2Ou','FPu2Ou','TPu1Ol','FPu1Ol','TPu2Ol','FPu2Ol',...
            'TPu1Ru','FPu1Ru','TPu2Ru','FPu2Ru','TPu1Rl','FPu1Rl','TPu2Rl','FPu2Rl');
    else
        save(fullfile(dirROC,'ROC.mat'),'T','TPb1','FPb1','TPu1','FPu1','TPb2','FPb2','TPu2','FPu2',...
            'TPu1O','FPu1O','TPu2O','FPu2O','TPu1R','FPu1R','TPu2R','FPu2R');
    end
else
    if compute_LU
        save(fullfile(dirROC,'ROC.mat'),'T','TPb1','FPb1','TPu1','FPu1','TPb2','FPb2','TPu2','FPu2',...
            'TPu1u','FPu1u','TPu2u','FPu2u','TPu1l','FPu1l','TPu2l','FPu2l');
    else
        if ~all_channels
            save(fullfile(dirROC,'ROC.mat'),'T','TPb1','FPb1','TPu1','FPu1','TPb2','FPb2','TPu2','FPu2');
        else
            save(fullfile(dirROC,'ROC.mat'),'T','TPb1','FPb1','TPu1','FPu1','TPb2','FPb2','TPu2','FPu2',...
                'TPu1A','FPu1A','TPu2A','FPu2A');
        end
    end
end
end 
out.NIRSmat = job.NIRSmat;
end

function [nB nu]  = count_TP_FP(IN,ch,t,alpha_bonf,alpha_unc,erdf,TPorFP,positive_ttest,byIter)

testBonferroni = 0;
%Inputs:
%IN: number of iterations
%ch: which channels 
%t: t-stats by iterations times channels
%alpha_bonf: Bonferroni corrected p-value
%alpha-unc: uncorrected p-value
%erdf: degrees of freedom, a vector of size IN
%TPorFP: Boolean, whether to compute true positives or false positives
%positive_ttest: Boolean, whether to perform a positive ttest or a negative ttest 
%by_Iter: Boolean, whether to test and to output for each protocole (iteration) 
%
%Outputs:
%nB: for each channel, percentage of TP or FP, Bonferroni corrected
%nu: same as nB, but not corrected for multiple counting
if byIter
    nB = zeros(length(ch),IN); nu = zeros(length(ch),IN);
else
    nB = zeros(length(ch),1); nu = zeros(length(ch),1);
end
for i=1:IN
    if testBonferroni
    th_z = spm_invTcdf(1-alpha_bonf, erdf(i));
    end
    th_zu = spm_invTcdf(1-alpha_unc, erdf(i));
    k = 0; %channel counter
    for j=ch
        k = k+1;
        if TPorFP
            if length(positive_ttest) > 1
                p_ttest = positive_ttest(k);
            else
                p_ttest = positive_ttest;
            end
            if p_ttest
                %one-sided test for true positives
                %careful, might need to change > for < depending on expected
                %sign of response
                if testBonferroni
                if t(i,j) > th_z
                    if byIter
                        nB(k,i) = 1;
                    else
                        nB(k) = nB(k)+1;
                    end
                end
                end
                if t(i,j) > th_zu
                    if byIter
                        nu(k,i) = 1;
                    else
                        nu(k) = nu(k)+1;
                    end
                end              
            else
                if testBonferroni
                if t(i,j) < -th_z
                    if byIter
                        nB(k,i) = 1;
                    else
                        nB(k) = nB(k)+1;
                    end
                end
                end
                if t(i,j) < -th_zu
                    if byIter
                        nu(k,i) = 1;
                    else
                        nu(k) = nu(k)+1;
                    end
                end
            end
        else
            %two-sided test - for false positives
            if testBonferroni
            if abs(t(i,j)) > th_z
                if byIter
                    nB(k,i) = 1;
                else
                    nB(k) = nB(k)+1;
                end
            end
            end
            if abs(t(i,j)) > th_zu
                if byIter
                    nu(k,i) = 1;
                else
                    nu(k) = nu(k)+1;
                end
            end
        end
    end
end
if ~byIter                   
    nB = nB/IN;
    nu = nu/IN;
end
end

