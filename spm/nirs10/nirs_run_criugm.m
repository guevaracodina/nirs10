function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

%Big loop over all subjects
outNIRSmat = {};
sN = size(job.subj,2);
for is=1:sN
    age = job.subj(1,is).age1;
    sDtp = job.subj(1,is).subj_path{:};
    staxp = job.subj(1,is).helmet.text_brainsight{:};
    
    %Reinitialize NIRS matrix for each subject
    NIRS = [];
    NIRS.Dt.s.age = age;
    NIRS.Dt.s.p = sDtp;
    NIRS.Dt.fir.stax.n = 'Brainsight(c)';
    NIRS.Dt.fir.stax.p{1} = staxp;
    save(fullfile(sDtp, 'NIRS.mat'),'NIRS');
    
    if ~strcmp(staxp,'')%if job.LESCA==0 %if LESCA==1 then done for each nirs_file
        % One helmet for all the selected acquisitions
        jobH.subj.sDtp = sDtp;
        jobH.subj.helmet.staxp = staxp;
        outH = nirs_criugm_getHelmet(jobH);% from Brainsight

        fig=findall(0,'name','Get positions from Brainsight (clbon)');
        waitfor(fig,'BeingDeleted','On');       
    end
    
    NIRS = [];
    load(fullfile(sDtp, 'NIRS.mat'),'NIRS');
    
    %Loop over all acquisitions
    aN = size(job.subj(1,is).acquisition.nirs_file,1);
    for ia=1:aN
        fp = job.subj(1,is).acquisition.nirs_file(ia,1);
        f = load(fp{:},'-mat');
        
        NIRS.Dt.fir.pp(1).p{ia,1} = fp{:};
        %NIRS.Dt.fir.pp(1).m{ia,1} Should not be required anymore?
        NIRS.Dt.fir.pp(1).m{ia,1} = job.subj(1,is).baseline_method;
        NIRS.Dt.fir.pp(1).pre = 'readCriugm';
        NIRS.Dt.fir.pp(1).job = job;
        
        if strcmp(staxp,'')%job.LESCA==1
            % sources
            NIRS.Cf.H.S.N = f.SD.nSrcs;
            for i = 1:f.SD.nSrcs
                n{1,i} = ['S',int2str(i)];
            end
            NIRS.Cf.H.S.n = n;
        end
        save(fullfile(sDtp,'NIRS.mat'),'NIRS');
        
        % CW system (1 for all aquisitions)
        CWsystem = job.subj(1,is).acquisition.CWsystem;
        if CWsystem == 5
        elseif CWsystem == 6
            %%%% attention la matrice NIRS est prevue pour ne prendre qu'une
            %%%% seule valeur pour les longueurs d'onde et pour la
            %%%% frequence ce qui signifie une seule machine pour toutes
            %%%% les expériences lancées dans un même job....
            job1.nirs_file = f;
            job1.sDtp = sDtp;
            out = nirs_criugm_readtechenCW6(job1);
        end
        clear f
    end
    outNIRSmat = [outNIRSmat; fullfile(sDtp,'NIRS.mat')];
end
out.NIRSmat = outNIRSmat;
end