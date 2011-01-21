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
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        beta(ROCiter,:,:)   = SPM.beta;
                        Bvar(ROCiter,:,:)   = full(SPM.xX.Bvar);
                    end

                    %get t stat for each iteration
                    r = 2; t = zeros(IN,sz_beta(2));
                    for i=1:IN
                        t(i,:) = squeeze(beta(i,r,:)) ./sqrt(squeeze(Bvar(i,r,:)));
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
                    for ROCiter=1:IN
                        testFullName = [testName int2str(ROCiter)];
                        dir_spm = [dir1 filesep testFullName filesep dir_stat];
                        SPM = [];
                        load(fullfile(dir_spm,'SPM.mat'));
                        beta(ROCiter,:,:)   = SPM.beta;
                        Bcov(ROCiter,:,:)   = SPM.xX.Bcov;
                        ResSS(ROCiter,:)  = SPM.ResSS;
                        trRV(ROCiter)   = SPM.xX.trRV;
                        trRVRV(ROCiter) = SPM.xX.trRVRV;
                        erdf(ROCiter)   = trRV(ROCiter)^2/trRVRV(ROCiter);
                    end
            % %         %to get a t-stat, for example for regressor r and channel c,
            % %         %and iteration i
            % %         r = 1; c = 8; i =1;
            % %         t = beta(i,r,c)/sqrt(ResSS(i,c) * Bcov(i,r,r)/trRV(i));
            % %         %get t stat for each iteration
            % %         r = 1; c = 7; t = zeros(IN,1);
            % %         for i=1:IN
            % %             t(i) = beta(i,r,c)/sqrt(ResSS(i,c) * Bcov(i,r,r)/trRV(i));
            % %         end

                    %get t stat for each iteration
                    r = 1; t = zeros(IN,sz_beta(2));
                    for i=1:IN
                        t(i,:) = squeeze(beta(i,r,:)) ./sqrt(squeeze(ResSS(i,:)'*Bcov(i,r,r)/trRV(i)));
                    end
                end
                
                %Add bonferroni and choice of t-stat value.
                alpha_unc = 0.05; %uncorrected threshold

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Sensitivity = true positives
                %Will be on channels 7 and 8
                NcTP = 2; %Number of true positive channels
                TPn = IN; %NcTP*IN; %Number of data points for true positives and false negatives
                alpha_bonf_TPn = alpha_unc/TPn;

                %[TPb{Jidx,Idx} TPu{Jidx,Idx}] = count_TP_FP(IN,[1:2],t,alpha_bonf_TPn,alpha_unc,erdf,true,true); 
                [TPb{Jidx,Idx} TPu{Jidx,Idx}] = count_TP_FP(IN,[11:20],...
                    t,alpha_bonf_TPn,alpha_unc,erdf,true,true,false); 

                %Specificity = 1 - false positives
                %Will be on all other channels than 7 and 8
                FPn = IN; %(sz_beta(2)-NcTP)*IN; %Number of data points for false positives 
                %and true negatives
                alpha_bonf_FPn = alpha_unc/FPn;

                [FPb{Jidx,Idx} FPu{Jidx,Idx}] = count_TP_FP(IN,[1:10],...
                    t,alpha_bonf_FPn,alpha_unc,erdf,false,true,false); 
                
                %tFPb = FPb; tFPu = FPu; tTPb = TPb; tTPu = TPu; 
                %FPb = FPb{Jidx,:}; tFPu = FPu{Jidx,:}; tTPb = TPb{Jidx,:}; tTPu = TPu{Jidx,:};                  
                save(fullfile(dirROC,['ROC' testName int2str(Jidx) '.mat']),...
                    'TPb','FPb','TPu','FPu');
                %FPb = tFPb; FPu = tFPu; TPb = tTPb; TPu = tTPu;
                
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
figure;
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
linespec{4,1} = '*';
linespec{5,1} = 'x';
linespec{6,1} = 's';
linespec{7,1} = 'd';
linespec{8,1} = '^';


for Idx=1:nSubj    
    for Jidx= 1:nJob
        plot(mean(FPb{Jidx,Idx}),mean(TPb{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
             mean(FPu{Jidx,Idx}),mean(TPu{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
    end   
end
hold off
%max(reshape(cell2mat(TPu),10,[]))
%br = squeeze(beta(:,2,11:20))./squeeze(beta(:,1,11:20));
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

