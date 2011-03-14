function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

outNIRSmat ={};

%Big loop over all subjects
sN = size(job.subj,2);
for is=1:sN
    age = job.subj(1,is).age1;
    sDtp = job.subj(1,is).subj_path{:};
    
    %Reinitialize NIRS matrix for each subject
    NIRS = [];
    NIRS.Dt.s.age = age;
    NIRS.Dt.s.p = sDtp;
    UN = size(job.subj(1,is).nirs_files,1);
    
    % Protocol
    if ~isempty(job.protocol{:})
        
        %Ignore parametric modulations - cf spm_run_fmri_design.m
        P.name = 'none';
        P.h    = 0;
        
        for iU=1:UN
            load(job.protocol{:});
            for kk = 1:size(names, 2)
                NIRS.Dt.fir.Sess(iU).U(kk).name = names(kk);
                NIRS.Dt.fir.Sess(iU).U(kk).ons = onsets{kk};
                NIRS.Dt.fir.Sess(iU).U(kk).dur = durations{kk};
                NIRS.Dt.fir.Sess(iU).U(kk).P = P;
            end
        end
    end
    
    % Anatomical image
    if ~isempty(job.subj(1,is).anatT1{:})
        NIRS.Dt.ana.T1 = job.subj(1,is).anatT1{:};
    end
    
    % Helmet
    if isfield(job.subj(1,is).helmet,'text_brainsight')
        staxp =job.subj(1,is).helmet.text_brainsight{:};
        NIRS.Dt.fir.stax.n = 'Brainsight(c)';
        NIRS.Dt.fir.stax.p{1} = job.subj(1,is).helmet.text_brainsight{:};
        if ~isempty(strfind(job.subj(1,is).helmet.text_brainsight{:},'template'))
            %%% CB : etude Said.... a mettre coherent /////
            % coordinates
            load(fullfile(dir_nt,'Hcoregistered.mat'));
            NIRS.Cf.H = Hcoregistered;
        end
    elseif isfield(job.subj(1,is).helmet,'T1_vitamins')
        NIRS.Dt.fir.stax.n = 'T1_vitamins';
        %read nirs file if already specified
        %         try catch
    elseif isfield(job.subj(1,is).helmet,'no_helmet')
        NIRS.Dt.fir.stax.n = 'no_helmet';
    end
    
    % Topo Data
    if ~isempty(job.subj(1,is).TopoData{:})
        helmetdone=1;
        % nirs_run_coreg has already been executed to generated once for all the
        % TopoData matrix
    else
        helmetdone=0;
    end
    
    save(fullfile(sDtp, 'NIRS.mat'),'NIRS');
    NIRS =[];
    
    if ~helmetdone
        jobH.subj.sDtp = sDtp;
        jobH.subj.helmet.staxp = staxp;
        outH = nirs_criugm_getHelmet(jobH);% get helmet configuration (S,D,P,Q) from Brainsight
    end
    fig=findall(0,'name','Get positions from Brainsight (clbon)');
    waitfor(fig,'BeingDeleted','On');
    
    load(fullfile(sDtp, 'NIRS.mat'));
    %Loop over all sessions
    for iU=1:UN
        fp = job.subj(1,is).nirs_files(iU,1);
        f = load(fp{:},'-mat');
        
        NIRS.Dt.fir.pp(1).p{iU,1} = fp{:};
        NIRS.Dt.fir.pp(1).m{iU,1} = job.subj(1,is).baseline_method;
        NIRS.Dt.fir.pp(1).pre = 'readCriugm';
        NIRS.Dt.fir.pp(1).job = job;
        
        save(fullfile(sDtp, 'NIRS.mat'),'NIRS');
        
        % test sur la machine utilisee
        job1.system = job.subj(1,is).CWsystem;
        job1.nirs_file = f;
        job1.sDtp = sDtp;
        out = nirs_criugm_readtechen(job1);% get C configuration from nirs files
        clear f
    end
    outNIRSmat = [outNIRSmat; fullfile(sDtp,'NIRS.mat')];
end
out.NIRSmat = outNIRSmat;
end