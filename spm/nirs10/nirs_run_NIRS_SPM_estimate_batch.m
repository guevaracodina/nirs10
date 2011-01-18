function out = nirs_run_NIRS_SPM_estimate_batch(job)
%NIRS_SPM GLM estimation - first level
which_GLM = job.NIRS_SPM_which_GLM;

%Loop over all subjects
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        NC = NIRS.Cf.H.C.N;
        try
            switch which_GLM
                case 1 %first
                    fGLM = NIRS.SPM(1);
                case 2
                    %need to loop
                    fGLM = NIRS.SPM;
                case 3
                    fGLM = NIRS.SPM(end);
            end
        catch
            try 
                fGLM = NIRS.SPM(1);
            catch
                disp(['Could not find a GLM for subject ' int2str(Idx)]);
            end
        end
        
        %loop over GLMs to estimate for a given subject and set of sessions
        %(usually only one such GLM)
        for g=1:length(fGLM)
            SPM = [];
            iSPM = 0; %count number of subsessions
            try
                load(fullfile(fGLM{g},'SPM.mat'));
            end
                    
            try
                SPM.xY.P;
            catch
                disp(['Could not find data file for subject ' int2str(Idx)]); 
            end
            
            %Prepare SPM matrix -- if fields are present / removing field
            try
                if isfield(SPM.xVi, 'V') == 1
                    SPM = rmfield(SPM, 'xVi');
                    precolor = 1;                   
                elseif isfield(SPM_nirs.xVi, 'V') == 0                   
                    precolor = 0;
                end
            catch
                disp('No V field present - GLM already estimated. Estimating it again!');
                precolor = 1;
            end
            %loop over sessions
            nsess = length(SPM.xY.P);
            for s=1:nsess
                d = fopen_NIR(SPM.xY.P{s},NC);
                %loop over subsessions, defined as period in between
                %movement intervals -- to do later
                
                %Strategy is that precoloring and prewhitening code from
                %NIRS_SPM are designed to work with single sessions, and
                %also long sessions are inefficient due to large covariance
                %matrix
                
                %Extract SPM_NIRS K filter
                
                %Extract
                %Need to transpose
                Y = d';
                
                %carefully extract SPM info required to pass to precoloring
                %or prewhitening
                tSPM = [];
                tSPM.Sess = SPM.Sess(s);
                %tSPM = SPM;
                svec = SPM.xX.K(s).K.row;
                %find elements of X for session s
                nbeta = size(SPM.xX.X,2);
                nbetaS = (nbeta-nsess)/nsess;
                %last entry is the constant regressor
                beta = [(s-1)*nbetaS+1:nbetaS nbeta-nsess+s];
                
                tSPM.xX.X = SPM.xX.X(svec,beta);
                tSPM.xX.K = SPM.xX.K(s).K;
                tSPM.xX.K.row = 1:length(svec);
                
                if precolor
                    tSPM = precoloring_batch(tSPM, Y);
                else
                    %not done yet
                    tSPM = prewhitening(tSPM, Y);
                end
                %Add piece of SPM to the whole SPM               
                iSPM = iSPM + 1;
                SPM.xXn(iSPM) = tSPM.xX;
            end 
            %Reconstruct SPM
            SPM.xX.trRV   = 0;
            SPM.xX.trRVRV = 0;
            SPM.xX.Bcov   = [];
            SPM.beta      = 0;
            SPM.ResSS     = 0;
            for n1=1:length(SPM.xXn)
                SPM.xX.trRV   = SPM.xX.trRV   + SPM.xXn(iSPM).trRV;
                SPM.xX.trRVRV = SPM.xX.trRVRV + SPM.xXn(iSPM).trRVRV;
                SPM.xX.Bcov   = blkdiag(SPM.xX.Bcov,sparse(SPM.xXn(iSPM).Bcov));
                SPM.beta      = SPM.beta+ SPM.xXn(iSPM).beta; %incorrect
                SPM.ResSS     = SPM.ResSS + SPM.xXn(iSPM).ResSS;
            end
            SPM.xX.erdf = (SPM.xX.trRV)^2/SPM.xX.trRVRV;
            %
            %SPM.beta = tSPM.beta
            save(fullfile(fGLM{g},'SPM.mat'), 'SPM');
        end
    catch
        disp(['Could not estimate GLM for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;