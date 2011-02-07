function out = nirs_run_ROCtest(job)
%specify some of the directory structure
dir_dataSPM = 'dataSPM';
%dir_stat = 'StatV';
IN = job.ROCiternum;
run_GLM = 0;
run_ROC = 1;
nSubj = size(job.NIRSmat,1);
nJob = size(job.ROCLoopJob,1);

%Loop over subjects
for Idx=1:nSubj
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir1,~,~] = fileparts(job.NIRSmat{Idx,1});
        %dirROC = dir1;
        x = findstr(filesep,dir1);
        dirROC = dir1(1:x(end));
        
        dir1 = [dir1 filesep dir_dataSPM];
        
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
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bvar   = zeros(IN,sz_Bvar(1),sz_Bvar(2));
                    RE = zeros(IN,sz_beta(1),sz_beta(2)/2);
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        beta(ROCiter,:,:)   = SPM.beta;
                        Bvar(ROCiter,:,:)   = full(SPM.xX.Bvar);
                        NIRS = [];
                        load(fullfile(dir_spm,'NIRS.mat'));
                        %1st Volterra kernel
                        RE(ROCiter,1,:) = SPM.beta(1,1:(sz_beta(2)/2))./ NIRS.Dt.fir.a -1;
                        %2nd Volterra kernel
                        try
                            RE(ROCiter,2,:) = SPM.beta(2,1:(sz_beta(2)/2))./ (NIRS.Dt.fir.b*NIRS.Dt.fir.a) -1;
                        end
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
                    beta   = zeros(IN,sz_beta(1),sz_beta(2));
                    Bcov   = zeros(IN,sz_Bcov(1),sz_Bcov(1));
                    ResSS  = zeros(IN,sz_ResSS(2));
                    trRV   = zeros(IN,1);
                    trRVRV = zeros(IN,1);
                    erdf   = zeros(IN,1);
                    RE = zeros(IN,sz_beta(1),sz_beta(2)/2);
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        %spm_DesRep('DesRepUI',SPM);
                        beta(ROCiter,:,:)   = SPM.beta;
                        Bcov(ROCiter,:,:)   = SPM.xX.Bcov;
                        ResSS(ROCiter,:)  = SPM.ResSS;
                        trRV(ROCiter)   = SPM.xX.trRV;
                        trRVRV(ROCiter) = SPM.xX.trRVRV;
                        erdf(ROCiter)   = trRV(ROCiter)^2/trRVRV(ROCiter);
                        NIRS = [];
                        load(fullfile(dir_spm,'NIRS.mat'));
                        %1st Volterra kernel
                        RE(ROCiter,1,:) = SPM.beta(1,1:(sz_beta(2)/2))./ NIRS.Dt.fir.a -1;
                        %2nd Volterra kernel
                        try
                            RE(ROCiter,2,:) = SPM.beta(2,1:(sz_beta(2)/2))./ (NIRS.Dt.fir.b*NIRS.Dt.fir.a) -1;
                            %for 3rd component, calculate the ratio of 2nd
                            %to 1st Volterra estimated amplitudes, less simulated 
                            RE(ROCiter,3,:) = SPM.beta(2,1:(sz_beta(2)/2))./ ...
                                SPM.beta(1,1:(sz_beta(2)/2)) - NIRS.Dt.fir.b;
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
                t1 = squeeze(T{Jidx,Idx}.t(:,1,:));
                t2 = squeeze(T{Jidx,Idx}.t(:,2,:));
                %Add bonferroni and choice of t-stat value.
                alpha_unc = 0.05; %uncorrected threshold
                
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
                    for exp_alpha_unc = 1:25
                        alpha_unc = 10.^(-exp_alpha_unc/5);
                        %1st Volterra
                        TPn = 20; %NcTP*IN; %Number of data points for true positives and false negatives
                        alpha_bonf_TPn = alpha_unc/TPn;
                        [u v] = count_TP_FP(IN,[1:20],...
                            t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,false); 
                        TPb1{Jidx,Idx} = [TPb1{Jidx,Idx} mean(u)];
                        TPu1{Jidx,Idx} = [TPu1{Jidx,Idx} mean(v)];
                        %2nd Volterra
                        [u v] = count_TP_FP(IN,[1:20],...
                            t2,alpha_bonf_TPn,alpha_unc,erdf,true,false,false); 
                        TPb2{Jidx,Idx} = [TPb2{Jidx,Idx} mean(u)];
                        TPu2{Jidx,Idx} = [TPu2{Jidx,Idx} mean(v)];
                        %False Positive
                        %1st Volterra
                        FPn = 20;
                        alpha_bonf_FPn = alpha_unc/FPn;
                        [u v ] = count_TP_FP(IN,[21:40],...
                            t1,alpha_bonf_FPn,alpha_unc,erdf,false,false,false); 
                        FPb1{Jidx,Idx} = [FPb1{Jidx,Idx} mean(u)];
                        FPu1{Jidx,Idx} = [FPu1{Jidx,Idx} mean(v)];
                        %2nd Volterra
                        [u v ] = count_TP_FP(IN,[21:40],...
                            t2,alpha_bonf_FPn,alpha_unc,erdf,false,false,false); 
                        FPb2{Jidx,Idx} = [FPb2{Jidx,Idx} mean(u)];
                        FPu2{Jidx,Idx} = [FPu2{Jidx,Idx} mean(v)];
                    end
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Sensitivity = true positives
                   
                    TPn = 20; %NcTP*IN; %Number of data points for true positives and false negatives
                    alpha_bonf_TPn = alpha_unc/TPn;

                    %[TPb{Jidx,Idx} TPu{Jidx,Idx}] = count_TP_FP(IN,[1:2],t,alpha_bonf_TPn,alpha_unc,erdf,true,true); 
                    %1st Volterra
                    [TPb1{Jidx,Idx} TPu1{Jidx,Idx}] = count_TP_FP(IN,[1:20],...
                        t1,alpha_bonf_TPn,alpha_unc,erdf,true,true,false); 
                    %2nd Volterra
                    [TPb2{Jidx,Idx} TPu2{Jidx,Idx}] = count_TP_FP(IN,[1:20],...
                        t2,alpha_bonf_TPn,alpha_unc,erdf,true,true,false); 
                    
                    %Specificity = 1 - false positives
                    FPn = 20; %(sz_beta(2)-NcTP)*IN; %Number of data points for false positives 
                    %and true negatives
                    alpha_bonf_FPn = alpha_unc/FPn;
                    %1st Volterra
                    [FPb1{Jidx,Idx} FPu1{Jidx,Idx}] = count_TP_FP(IN,[21:40],...
                        t1,alpha_bonf_FPn,alpha_unc,erdf,false,true,false); 
                    %2nd Volterra
                    [FPb2{Jidx,Idx} FPu2{Jidx,Idx}] = count_TP_FP(IN,[21:40],...
                        t2,alpha_bonf_FPn,alpha_unc,erdf,false,true,false); 

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
if generate_ROC_curves
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
    boxplot(tmp); 
    
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
    boxplot(tmp); 
    %2nd Volterra
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
        end   
    end
    boxplot(tmp); 
    
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
    boxplot(tmp); 
    %ratio of 2nd to 1st estimated Volterra less simulated
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,3,:),3))];          
        end   
    end
    boxplot(tmp); 
    
    
    %t-stat
    %1st Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
        end   
    end
    boxplot(tmp); 
    %2nd Volterra
    figure;
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
        end   
    end
    boxplot(tmp); 
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
    boxplot(tmp); 
    %2nd Volterra
    subplot(1,2,2);
    tmp = [];    
    for Idx=1:nSubj    
        for Jidx= 1:nJob
           tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
        end   
    end
    boxplot(tmp); 
end
%max(reshape(cell2mat(TPu),10,[]))
%br = squeeze(beta(:,2,11:20))./squeeze(beta(:,1,11:20));
%

%figures for box-whisker plot

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
            if positive_ttest
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

