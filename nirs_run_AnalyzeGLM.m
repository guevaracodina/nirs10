function out = nirs_run_AnalyzeGLM(job)
%specify some of the directory structure
dir_dataSPM = 'dataSPM';
nSubj = size(job.NIRSmat,1);
nJob = size(job.ROCLoopJob,1);

%Loop over subjects
for Idx=1:nSubj
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir1,dummy,dummy2] = fileparts(job.NIRSmat{Idx,1});
        %x = findstr(filesep,dir1);
        %dirROC = dir1(1:x(end));
        
        dir1 = [dir1 filesep dir_dataSPM];
        
        %loop over jobs
        for Jidx=1:nJob
            %loop over the specified job
            clear LoopJob
            %Load/read job
            LoopJob = load(job.ROCLoopJob{Jidx});
            %replace whatever is specified in the LoopJob NIRSmat by
            %subject Idx
            %LoopJob.matlabbatch{1}.spm.tools.nirs10.

            %build folder name
            statdir = LoopJob.matlabbatch{1}.spm.tools.nirs10. ...
                model_specify.wls_bglm_specify.dir1;
            dir_spm = [dir1 filesep statdir];
            
            SPM = [];
            load(fullfile(dir_spm,'SPM.mat'));
            NIRS = [];
            load(fullfile(dir_spm,'NIRS.mat'));
            NC = NIRS.Cf.H.C.N;
           
            %loop over sessions - old ROCiter becomes here the number of
            %sessions
            r = 2;
            nSess = length(SPM.xXn);
            clear t
            for f = 1:nSess
                t(f,:) = SPM.xXn{f}.t(r,:);
            end
            
            try 
                SPM.xX.Bvar;
                %very approximate value for erdf - assume enough
                %degrees of freedom to approximate t-stat by Z-score
                erdf   = 1000*ones(nSess,1);
            catch
               
                erdf   = zeros(nSess,1);
                for f=1:nSess                
                    trRV   = SPM.xXn{f}.trRV;
                    trRVRV = SPM.xXn{f}.trRVRV;
                    erdf(f)   = trRV^2/trRVRV;
                end
            end
            %Add bonferroni and choice of t-stat value.
            alpha_unc = 0.05; %uncorrected threshold

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Sensitivity = true positives
            nc = NC/2; 
            %channel names
            cn{Idx} = NIRS.Cf.H.C.n;
            alpha_bonf_TPn = alpha_unc/nc;

            [HbOb{Jidx,Idx} HbOu{Jidx,Idx}] = count_TP_FP(nSess,1:nc,...
                t,alpha_bonf_TPn,alpha_unc,erdf,false,true,true); 

            [HbRb{Jidx,Idx} HbRu{Jidx,Idx}] = count_TP_FP(nSess,nc+1:NC,...
                t,alpha_bonf_TPn,alpha_unc,erdf,false,false,true); 
            
            %to display a design matrix:
            %spm_DesRep('DesRepUI',SPM)

            
        end %end for Jidx
    catch exception
        disp(['Could not perform AnalyzeGLM iterations subject ' int2str(Idx)]);
        disp(exception.identifier);
    end   
end %end for Idx

linespec{1,1} = '.';
linespec{2,1} = '+';
linespec{3,1} = '*';
linespec{1,2} = 'o';
linespec{2,2} = 'x';
linespec{3,2} = 's';
linespec{1,3} = 'd';
linespec{2,3} = '^';
linespec{3,3} = 'v';
linespec{1,4} = '<';
linespec{2,4} = '>';
linespec{3,4} = 'p';
linespec{1,5} = 'h';

linecolor{1} = 'b';
linecolor{2} = 'r';
linecolor{3} = 'k';
linecolor{4} = 'g';
linecolor{5} = 'y';


% linespec{1,1} = '.';
% linespec{2,1} = '+';
% linespec{3,1} = 'o';
% linespec{4,1} = '*';
% linespec{5,1} = 'x';
% linespec{6,1} = 's';
% linespec{7,1} = 'd';
% linespec{8,1} = '^';
% linespec{9,1} = 'p';
% linespec{10,1} = 'h';


figure;
for Idx=1:nSubj    
    for Jidx= 1:nJob
        plot(sum(HbOu{Jidx,Idx},2),[linespec{Jidx,Idx} linecolor{Idx}]); hold on
    end   
end
hold off
figure;
for Idx=1:nSubj    
    for Jidx= 1:nJob
        plot(sum(HbRu{Jidx,Idx},2),[linespec{Jidx,Idx} linecolor{Idx}]); hold on
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

