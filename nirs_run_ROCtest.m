function out = nirs_run_ROCtest(job)
%specify some of the directory structure
dir_dataSPM = 'dataSPM';
%dir_stat = 'StatV';
IN = job.ROCiternum;
run_GLM = 0;
run_ROC = 1;
nSubj = size(job.NIRSmat,1);
nJob = size(job.ROCLoopJob,1);
Volt2_positive_ttest = false;

%Loop over subjects
for Idx=1:nSubj
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir0,~,~] = fileparts(job.NIRSmat{Idx,1});
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
                        catch
                        end
                    end
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
                
                try 
                    SPM.xX.Bvar;
                    %ROC
                    sz_beta = size(SPM.beta);
                    sz_Bvar = size(SPM.xX.Bvar);
                    try
                        if SPM.GenerateHbT
                            sz_beta(2) = 40; %quick fix for bug due to
                            %introducing HbT
                            sz_ResSS(2) = 40; %quick fix for bug due to introducing HbT
                        end
                    end
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bvar   = zeros(IN,sz_Bvar(1),sz_Bvar(2));
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
                                beta(ROCiter,:,:)   = SPM.beta(:,1:40); %quick fix for bug due to introducing HbT
                                Bvar(ROCiter,:,:)   = full(SPM.xX.Bvar(:,1:40)); %quick fix for bug due to introducing HbT
                            else
                                beta(ROCiter,:,:)   = SPM.beta; %(:,1:40); %quick fix for bug due to introducing HbT
                                Bvar(ROCiter,:,:)   = full(SPM.xX.Bvar) ; %(:,1:40)); %quick fix for bug due to introducing HbT
                            end
                        catch
                            beta(ROCiter,:,:)   = SPM.beta; %(:,1:40); %quick fix for bug due to introducing HbT
                            Bvar(ROCiter,:,:)   = full(SPM.xX.Bvar) ; %(:,1:40)); %quick fix for bug due to introducing HbT
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
                    sz_beta = size(SPM.beta);
                    
                    sz_Bcov = size(SPM.xX.Bcov);
                    sz_ResSS = size(SPM.ResSS);
                    try
                        if SPM.GenerateHbT
                            sz_beta(2) = 40; %quick fix for bug due to
                            %introducing HbT
                            sz_ResSS(2) = 40; %quick fix for bug due to introducing HbT
                        end
                    end
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bcov   = zeros(IN,sz_Bcov(1),sz_Bcov(1));
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
                        Bcov(ROCiter,:,:)   = SPM.xX.Bcov;
                        try
                            if SPM.GenerateHbT
                                beta(ROCiter,:,:)   = SPM.beta(:,1:40); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.ResSS(1,1:40); %quick fix for bug due to introducing HbT
                            else
                                beta(ROCiter,:,:)   = SPM.beta; %(:,1:40); %quick fix for bug due to introducing HbT
                                ResSS(ROCiter,:)  = SPM.ResSS; %(1,1:40); %quick fix for bug due to introducing HbT
                            end
                        catch
                            beta(ROCiter,:,:)   = SPM.beta; %(:,1:40); %quick fix for bug due to introducing HbT
                            ResSS(ROCiter,:)  = SPM.ResSS; %(1,1:40); %quick fix for bug due to introducing HbT
                        end
                        CORR(ROCiter,1) = corr(SPM.xX.X(:,1),SPM.xX.X(:,2));
                        trRV(ROCiter)   = SPM.xX.trRV;
                        trRVRV(ROCiter) = SPM.xX.trRVRV;
                        erdf(ROCiter)   = trRV(ROCiter)^2/trRVRV(ROCiter);
                        if isnan(erdf(ROCiter))
                            erdf(ROCiter) = 491;
                            trRV(ROCiter) = 343;
                            trRVRV(ROCiter) = 239;
                        end
                        NIRS = [];
                        load(fullfile(dir_spm,'NIRS.mat'));
                        %1st Volterra kernel
                        RE(ROCiter,1,:) = SPM.beta(1,1:(sz_beta(2)/2))./ NIRS.Dt.fir.a; % -1;
                        %2nd Volterra kernel
                        try
                            RE(ROCiter,2,:) = SPM.beta(2,1:(sz_beta(2)/2))./ (NIRS.Dt.fir.b*NIRS.Dt.fir.a);% -1;
                            %for 3rd component, calculate the ratio of 2nd
                            %to 1st Volterra estimated amplitudes, less simulated 
                            RE(ROCiter,3,:) = SPM.beta(2,1:(sz_beta(2)/2))./ ...
                                SPM.beta(1,1:(sz_beta(2)/2));%  - NIRS.Dt.fir.b;
                        end  
                        SNR(ROCiter,:) = NIRS.Dt.fir.SNR;
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
                alpha_unc = 0.05; %uncorrected threshold
                
                %options for the plots
                byIter = false;
                
                generate_ROC_curves = 1;
                if generate_ROC_curves
                    %need to loop over thresholds
                    TPb1{Jidx,Idx} = [];
                    TPu1{Jidx,Idx} = [];
                    FPb1{Jidx,Idx} = [];
                    FPu1{Jidx,Idx} = [];
                    TPb2{Jidx,Idx} = [];
                    TPu2{Jidx,Idx} = [];
                    FPb2{Jidx,Idx} = [];
                    FPu2{Jidx,Idx} = [];
                    %for 25th and 75th percentiles
                    TPu1u{Jidx,Idx} = [];
                    FPu1u{Jidx,Idx} = [];
                    TPu1l{Jidx,Idx} = [];
                    FPu1l{Jidx,Idx} = [];
                    TPu2u{Jidx,Idx} = [];
                    FPu2u{Jidx,Idx} = [];
                    TPu2l{Jidx,Idx} = [];
                    FPu2l{Jidx,Idx} = [];
                    %for HbO and HbR separately
                    TPu1O{Jidx,Idx} = [];
                    TPu1R{Jidx,Idx} = [];
                    TPu2O{Jidx,Idx} = [];
                    TPu2R{Jidx,Idx} = [];
                    FPu1O{Jidx,Idx} = [];
                    FPu1R{Jidx,Idx} = [];
                    FPu2O{Jidx,Idx} = [];
                    FPu2R{Jidx,Idx} = [];
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
                    
                    %
                    try
                        Volt2{1,1} = 0;
                        Volt2{2,1} = 0;
                        Volt2{3,1} = 0;
                        Volt2{4,1} = 0;
                        Volt2{5,1} = 0;
                        Volt2{6,1} = 0;
                        Volt2{7,1} = 1;
                        Volt2{8,1} = 1;
                        Volt2_positive_ttest = Volt2{Jidx,Idx};
                    end
                    %Volt2_positive_ttest = false;
                    
                    for exp_alpha_unc = [19 22 24 28 31 35 39 42 47 66 72 80] % [2:10 12:2:20 25]
                        alpha_unc = 10.^(-exp_alpha_unc/15);
                        %1st Volterra
                        TPn = 20; %NcTP*IN; %Number of data points for true positives and false negatives
                        alpha_bonf_TPn = alpha_unc/TPn;
                        [u v] = count_TP_FP(IN,[1:20],...
                            t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                        TPb1{Jidx,Idx} = [TPb1{Jidx,Idx} median(u)];
                        TPu1{Jidx,Idx} = [TPu1{Jidx,Idx} median(v)];
                        
                        TPu1u{Jidx,Idx} = [TPu1u{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu1l{Jidx,Idx} = [TPu1l{Jidx,Idx} median(v)-prctile(v,25)];
                        %2nd Volterra
                        [u v] = count_TP_FP(IN,[1:20],...
                            t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                        TPb2{Jidx,Idx} = [TPb2{Jidx,Idx} median(u)];
                        TPu2{Jidx,Idx} = [TPu2{Jidx,Idx} median(v)];
                        
                        TPu2u{Jidx,Idx} = [TPu2u{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu2l{Jidx,Idx} = [TPu2l{Jidx,Idx} median(v)-prctile(v,25)];
                        
                        %False Positive
                        %1st Volterra
                        FPn = 20;
                        alpha_bonf_FPn = alpha_unc/FPn;
                        [u v ] = count_TP_FP(IN,[21:40],...
                            t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPb1{Jidx,Idx} = [FPb1{Jidx,Idx} median(u)];
                        FPu1{Jidx,Idx} = [FPu1{Jidx,Idx} median(v)];
                        
                        FPu1u{Jidx,Idx} = [FPu1u{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu1l{Jidx,Idx} = [FPu1l{Jidx,Idx} median(v)-prctile(v,25)];
                        
                        %2nd Volterra
                        [u v ] = count_TP_FP(IN,[21:40],...
                            t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPb2{Jidx,Idx} = [FPb2{Jidx,Idx} median(u)];
                        FPu2{Jidx,Idx} = [FPu2{Jidx,Idx} median(v)];
                        
                        FPu2u{Jidx,Idx} = [FPu2u{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu2l{Jidx,Idx} = [FPu2l{Jidx,Idx} median(v)-prctile(v,25)];
                        
                        %For HbO and HbR separately
                        %1st Volterra
                        TPn = 10; %NcTP*IN; %Number of data points for true positives and false negatives
                        alpha_bonf_TPn = alpha_unc/TPn;
                        %HbO
                        [u v] = count_TP_FP(IN,[1:10],...
                            t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                        TPu1O{Jidx,Idx} = [TPu1O{Jidx,Idx} median(v)];                        
                        TPu1Ou{Jidx,Idx} = [TPu1Ou{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu1Ol{Jidx,Idx} = [TPu1Ol{Jidx,Idx} median(v)-prctile(v,25)];
                        %2nd Volterra
                        [u v] = count_TP_FP(IN,[1:10],...
                            t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                        TPu2O{Jidx,Idx} = [TPu2O{Jidx,Idx} median(v)];                      
                        TPu2Ou{Jidx,Idx} = [TPu2Ou{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu2Ol{Jidx,Idx} = [TPu2Ol{Jidx,Idx} median(v)-prctile(v,25)];
                        %False Positive
                        %1st Volterra
                        FPn = 10;
                        alpha_bonf_FPn = alpha_unc/FPn;
                        [u v] = count_TP_FP(IN,[21:30],...
                            t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPu1O{Jidx,Idx} = [FPu1O{Jidx,Idx} median(v)];
                        FPu1Ou{Jidx,Idx} = [FPu1Ou{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu1Ol{Jidx,Idx} = [FPu1Ol{Jidx,Idx} median(v)-prctile(v,25)];
                        %2nd Volterra
                        [u v ] = count_TP_FP(IN,[21:30],...
                            t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPu2O{Jidx,Idx} = [FPu2O{Jidx,Idx} median(v)];
                        FPu2Ou{Jidx,Idx} = [FPu2Ou{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu2Ol{Jidx,Idx} = [FPu2Ol{Jidx,Idx} median(v)-prctile(v,25)];
                        %HbR
                        [u v] = count_TP_FP(IN,[11:20],...
                            t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                        TPu1R{Jidx,Idx} = [TPu1R{Jidx,Idx} median(v)];                        
                        TPu1Ru{Jidx,Idx} = [TPu1Ru{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu1Rl{Jidx,Idx} = [TPu1Rl{Jidx,Idx} median(v)-prctile(v,25)];
                        %2nd Volterra
                        [u v] = count_TP_FP(IN,[11:20],...
                            t2,alpha_bonf_TPn,alpha_unc,erdf,true,Volt2_positive_ttest,byIter); 
                        TPu2R{Jidx,Idx} = [TPu2R{Jidx,Idx} median(v)];
                        
                        TPu2Ru{Jidx,Idx} = [TPu2Ru{Jidx,Idx} prctile(v,75)-median(v)];
                        TPu2Rl{Jidx,Idx} = [TPu2Rl{Jidx,Idx} median(v)-prctile(v,25)];
                        %False Positive
                        %1st Volterra
                        FPn = 10;
                        alpha_bonf_FPn = alpha_unc/FPn;
                        [u v] = count_TP_FP(IN,[31:40],...
                            t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPu1R{Jidx,Idx} = [FPu1R{Jidx,Idx} median(v)];
                        FPu1Ru{Jidx,Idx} = [FPu1Ru{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu1Rl{Jidx,Idx} = [FPu1Rl{Jidx,Idx} median(v)-prctile(v,25)];

                        %2nd Volterra
                        [u v ] = count_TP_FP(IN,[31:40],...
                            t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,byIter); 
                        FPu2R{Jidx,Idx} = [FPu2R{Jidx,Idx} median(v)];
                        FPu2Ru{Jidx,Idx} = [FPu2Ru{Jidx,Idx} prctile(v,75)-median(v)];
                        FPu2Rl{Jidx,Idx} = [FPu2Rl{Jidx,Idx} median(v)-prctile(v,25)];
                    end
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Sensitivity = true positives
                   
                    TPn = 20; %NcTP*IN; %Number of data points for true positives and false negatives
                    alpha_bonf_TPn = alpha_unc/TPn;

                    %[TPb{Jidx,Idx} TPu{Jidx,Idx}] = count_TP_FP(IN,[1:2],t,alpha_bonf_TPn,alpha_unc,erdf,true,true); 
                    %1st Volterra
                    [TPb1{Jidx,Idx} TPu1{Jidx,Idx}] = count_TP_FP(IN,[1:20],...
                        t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                    %2nd Volterra
                    [TPb2{Jidx,Idx} TPu2{Jidx,Idx}] = count_TP_FP(IN,[1:20],...
                        t2,alpha_bonf_TPn,alpha_unc,erdf,true,true,byIter); 
                    
                    %Specificity = 1 - false positives
                    FPn = 20; %(sz_beta(2)-NcTP)*IN; %Number of data points for false positives 
                    %and true negatives
                    alpha_bonf_FPn = alpha_unc/FPn;
                    %1st Volterra
                    [FPb1{Jidx,Idx} FPu1{Jidx,Idx}] = count_TP_FP(IN,[21:40],...
                        t1,alpha_bonf_FPn,alpha_unc,erdf,false,true,byIter); 
                    %2nd Volterra
                    [FPb2{Jidx,Idx} FPu2{Jidx,Idx}] = count_TP_FP(IN,[21:40],...
                        t2,alpha_bonf_FPn,alpha_unc,erdf,false,true,byIter); 

                    %tFPb = FPb; tFPu = FPu; tTPb = TPb; tTPu = TPu; 
                    %FPb = FPb{Jidx,:}; tFPu = FPu{Jidx,:}; tTPb = TPb{Jidx,:}; tTPu = TPu{Jidx,:};                  
                    %save(fullfile(dirROC,['ROC' testName int2str(Jidx) '.mat']),...
                    %    'TPb','FPb','TPu','FPu');
                    %FPb = tFPb; FPu = tFPu; TPb = tTPb; TPu = tTPu;
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
%figure; plot(FP,TP);
% if run_ROC
%     save(fullfile(dirROC,'ROC.mat'),'TPb','FPb','TPu','FPu');
% end

%ROC plots

% linespec{1,1} = '.';
% linespec{2,1} = '+';
% linespec{3%,1} = 'o';
% linespec{1,2} = '*';
% linespec{2,2} = 'x';
% linespec{3,2} = 's';
% linespec{1,3} = 'd';
% linespec{2,3} = '^';
% linespec{3,3} = 'v';
% linespec{1,4} = '<';
% linespec{2,4} = '>';
% linespec{3,4} = 'p';
% linespec{1,5} = 'h';


linespec{1,1} = '.';
linespec{2,1} = '+';
linespec{3,1} = 'o';
linespec{4,1} = '^';
linespec{5,1} = 's';
linespec{6,1} = 'p';
linespec{7,1} = 'd';
linespec{8,1} = '*';
linespec{9,1} = 'x';
linespec{10,1} = 'h';

linespec{1,1} = '-b';
linespec{2,1} = '--r';
linespec{3,1} = ':k';
linespec{4,1} = '-.g';
linespec{5,1} = '-k';
linespec{6,1} = '--g';
linespec{7,1} = ':k';
linespec{8,1} = '-.b';

linespec{1,1} = '-';
linespec{2,1} = ':';
linespec{3,1} = '--';
linespec{4,1} = '-.';
linespec{5,1} = '-';
linespec{6,1} = '--';
linespec{7,1} = ':';
linespec{8,1} = '-.';

scatspec{1,1} = 'ob';
scatspec{2,1} = '.k';
scatspec{3,1} = '+r';
scatspec{4,1} = 'sg';

% % set(gcf,'PaperUnits','centimeters')
% % %This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
% % xSize = 9; ySize = 7;
% % %These are my size variables, width of 8 and a height of 12, will be used a lot later.
% % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % %Additional coordinates to center the figure on A4-paper
% % set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% % %This command sets the position and size of the figure on the paper to the desired values.
% % set(gcf,'Position',[X Y xSize*50 ySize*50])

if generate_ROC_curves
    %Correlation
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.RE(:,3,:),3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,max(T{Jidx,Idx}.RE(:,3,:),[],3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.t(:,1,:),3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.t(:,2,:),3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,min(T{Jidx,Idx}.t(:,2,:),[],3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,median(T{Jidx,Idx}.t(:,2,:),3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,prctile(T{Jidx,Idx}.t(:,2,:),25,3),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.SNR,2),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    %[b,bint,r,rint,stats] = regress(mean(T{Jidx,Idx}.SNR,2),T{Jidx,Idx}.CORR)
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            scatter(T{Jidx,Idx}.CORR,10.^(mean(T{Jidx,Idx}.SNR,2)/10),scatspec{Jidx,Idx}); hold on
        end   
    end
    hold off
    
    %1st Volterra
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(FPb1{Jidx,Idx},TPb1{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %2nd Volterra
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(FPb2{Jidx,Idx},TPb2{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %uncorrected only
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},linespec{Jidx,Idx}); hold on
        end   
    end
    hold off
    
    %1st and 2nd Volterra grouped with subplot - with error bars
    figure;
    subplot(2,1,1);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
             %add error bars
            errorbar(FPu1{Jidx,Idx},TPu1{Jidx,Idx},TPu1l{Jidx,Idx},TPu1u{Jidx,Idx},linespec{Jidx,Idx}); hold on
            herrorbar(FPu1{Jidx,Idx},TPu1{Jidx,Idx},FPu1l{Jidx,Idx},FPu1u{Jidx,Idx},linespec{Jidx,Idx}); hold on
        end   
    end
    hold off
    subplot(2,1,2);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
             errorbar(FPu2{Jidx,Idx},TPu2{Jidx,Idx},TPu2l{Jidx,Idx},TPu2u{Jidx,Idx},linespec{Jidx,Idx}); hold on
             herrorbar(FPu2{Jidx,Idx},TPu2{Jidx,Idx},FPu2l{Jidx,Idx},FPu2u{Jidx,Idx},linespec{Jidx,Idx}); hold on
        end   
    end
    hold off
    
    %1st and 2nd Volterra grouped with subplot
    figure;
    subplot(2,1,1);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(FPb1{Jidx,Idx},TPb1{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    subplot(2,1,2);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(FPb2{Jidx,Idx},TPb2{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    
    %figure for SNR
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            m = median(T{Jidx,Idx}.SNR(:));
            errorbar(median(T{Jidx,Idx}.a2),m,m-prctile(T{Jidx,Idx}.SNR(:),25),prctile(T{Jidx,Idx}.SNR(:),75)-m); hold on
        end   
    end
    hold off
    %better as a boxplot?
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp T{Jidx,Idx}.SNR(:)];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
    %COMPARE HbO and HbR
    %1st and 2nd Volterra grouped with subplot - with error bars
    figure;
    subplot(2,1,1);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
             %add error bars
            errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
            herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    subplot(2,1,2);
    for Idx=1:nSubj    
        for Jidx= 1:nJob
             errorbar(FPu2O{Jidx,Idx},TPu2O{Jidx,Idx},TPu2Ol{Jidx,Idx},TPu2Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
             herrorbar(FPu2O{Jidx,Idx},TPu2O{Jidx,Idx},FPu2Ol{Jidx,Idx},FPu2Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
             errorbar(FPu2R{Jidx,Idx},TPu2R{Jidx,Idx},TPu2Rl{Jidx,Idx},TPu2Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
             herrorbar(FPu2R{Jidx,Idx},TPu2R{Jidx,Idx},FPu2Rl{Jidx,Idx},FPu2Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %Custom:
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
             %add error bars
            errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
            herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    
    %Custom:
    figure;
    for Idx=1:nSubj    
        for Jidx=[2 4]; % 1:nJob
             %add error bars
            errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
            errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx+2,Idx} 'b']); hold on
            herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx+2,Idx} 'b']); hold on
        end   
    end
    hold off
    
    %SNR vs t-stat
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            tmp= T{Jidx,Idx}.t(:,1,1:20);
            scatter(T{Jidx,Idx}.SNR(:),tmp(:));
            scatter(mean(T{Jidx,Idx}.SNR,2),squeeze(mean(tmp,3)));
        end
    end
    addpath('J:\NIRS_nonlinear');
    %careful: order of series is the opposite!!!
    XMatrix1 = [FPu1R{2,1};FPu1O{2,1};FPu1R{1,1};FPu1O{1,1};]';
    YMatrix1 = [TPu1R{2,1};TPu1O{2,1};TPu1R{1,1};TPu1O{1,1}]';
    LMatrix1 = [TPu1Rl{2,1};TPu1Ol{2,1};TPu1Rl{1,1};TPu1Ol{1,1}]';
    UMatrix1 = [TPu1Ru{2,1};TPu1Ou{2,1};TPu1Ru{1,1};TPu1Ou{1,1}]';        
    Figure3_createfigure(XMatrix1, YMatrix1, LMatrix1, UMatrix1);
    %Figure3_createfigure(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},
    %TPu1Ou{Jidx,Idx}
else
    %1st Volterra
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(mean(FPb1{Jidx,Idx}),mean(TPb1{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
                 mean(FPu1{Jidx,Idx}),mean(TPu1{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %2nd Volterra
    figure;
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(mean(FPb2{Jidx,Idx}),mean(TPb2{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
                 mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %1st and 2nd Volterra grouped as subplots
    figure;
    subplot(1,2,1);
    %1st Volterra
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(mean(FPb1{Jidx,Idx}),mean(TPb1{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
                 mean(FPu1{Jidx,Idx}),mean(TPu1{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    subplot(1,2,2);
    %2nd Volterra
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(mean(FPb2{Jidx,Idx}),mean(TPb2{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
                 mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
    %Uncorrected only, no colors
    for Idx=1:nSubj    
        for Jidx= 1:nJob
            plot(mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
        end   
    end
    hold off
      
    %Relative error
    %1st Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
        end   
    end
    boxplot(tmp); 
    %2nd Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
    %Relative error grouped as subplots
    %1st Volterra
    figure;
    subplot(1,2,1);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %2nd Volterra
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
    %Relative error for 1st Volterra and ratio of 2nd to 1st Volterra for
    %second plot
    figure;
    subplot(1,2,1);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %ratio of 2nd to 1st estimated Volterra less simulated
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,3,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %absolute value
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.RE(:,3,:),3)))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
    %t-stat
    %1st Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %2nd Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %absolute value
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.t(:,2,:),3)))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
    %1st and 2nd Volterra grouped as subplots
    %1st Volterra
    figure;
    subplot(1,2,1);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    %2nd Volterra
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
end
%max(reshape(cell2mat(TPu),10,[]))
%br = squeeze(beta(:,2,11:20))./squeeze(beta(:,1,11:20));
%
%4 plots combined
 %1st Volterra
    figure;
    subplot(2,2,1)
    %for Idx=1:nSubj    
        Jidx= 1;
            plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-b',...
                FPu1{Jidx,Idx},TPu1{Jidx,Idx},'ob'); hold on
               
          Jidx = 2;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'--r',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
             
             Jidx = 3;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.g',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+g'); hold on
             
             Jidx = 4;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'--k',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+k'); hold on
%               Jidx = 5;
%            plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.r',...
%                  FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
%              
%              Jidx = 6;
%            plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},':b',...
%                  FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+b'); hold on
%    %   end   
    %end
    hold off   
    subplot(2,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    subplot(2,2,3);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    subplot(2,2,4);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp T{Jidx,Idx}.SNR(:)];          
        end   
    end
    boxplot(tmp,'notch','on'); 
%figures for box-whisker plot
%3 plots combined
 %1st Volterra
    figure;
    subplot(1,3,1)
    %for Idx=1:nSubj    
        Jidx= 1;
            plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-b',...
                FPu1{Jidx,Idx},TPu1{Jidx,Idx},'ob'); hold on
               
          Jidx = 2;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'--r',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
                Jidx = 3;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.g',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+g'); hold on
             
             Jidx = 4;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},':k',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+k'); hold on
                 Jidx = 5;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.r',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
             
             Jidx = 6;
           plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},':b',...
                 FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+b'); hold on
    %   end   
    %end
    hold off   
    subplot(1,3,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    subplot(1,3,3);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 

    %4 plots combined
 %2nd Volterra
    figure;
    subplot(2,2,1)
    %for Idx=1:nSubj    
        Jidx= 1;
            plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-b',...
                FPu2{Jidx,Idx},TPu2{Jidx,Idx},'ob'); hold on
               
          Jidx = 2;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--r',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
             
             Jidx = 3;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.g',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+g'); hold on
             
             Jidx = 4;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--k',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+k'); hold on
              Jidx = 5;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.r',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
             
             Jidx = 6;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':b',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+b'); hold on
             
            
    %   end   
    %end
    hold off   
    subplot(2,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    subplot(2,2,3);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.t(:,2,:),3)))];          
        end   
    end
     hold off   
     boxplot(tmp,'notch','on'); 
    subplot(2,2,4);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,3,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
%figures for box-whisker plot
%3 plots combined
 %2nd Volterra
    figure;
    subplot(1,3,1)
    %for Idx=1:nSubj    
        Jidx= 1;
            plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-b',...
                FPu2{Jidx,Idx},TPu2{Jidx,Idx},'ob'); hold on
               
          Jidx = 2;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--r',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
                Jidx = 3;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.g',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+g'); hold on
             
             Jidx = 4;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':k',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+k'); hold on
                Jidx = 5;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.r',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
             
             Jidx = 6;
           plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':b',...
                 FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+b'); hold on
    %   end   
    %end
    hold off   
    subplot(1,3,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    subplot(1,3,3);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
        end   
    end
    boxplot(tmp,'notch','on'); 
    
out.NIRSmat = job.NIRSmat;
end

function [nB nu]  = count_TP_FP(IN,ch,t,alpha_bonf,alpha_unc,erdf,TPorFP,positive_ttest,byIter)
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
    th_z = spm_invTcdf(1-alpha_bonf, erdf(i));
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
                if t(i,j) > th_z
                    if byIter
                        nB(k,i) = 1;
                    else
                        nB(k) = nB(k)+1;
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
                if t(i,j) < -th_z
                    if byIter
                        nB(k,i) = 1;
                    else
                        nB(k) = nB(k)+1;
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
            if abs(t(i,j)) > th_z
                if byIter
                    nB(k,i) = 1;
                else
                    nB(k) = nB(k)+1;
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

