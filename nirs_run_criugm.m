function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

outNIRSmat ={};

% template directory
dir_nt = 'D:\Users\Clément\Projets_CRIUGM\nirs10_templates';

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
        for iU=1:UN
            NIRS.Dt.ana.rend{iU,1} = job.protocol{:};
        end
    end
    
    % Anatomical image
    if ~isempty(job.subj(1,is).anatT1{:})
        NIRS.Dt.ana.T1 = job.subj(1,is).anatT1{:};
    end
    
    % Helmet
    if ~isempty(job.subj(1,is).text_brainsight{:})
        NIRS.Dt.fir.stax.n = 'Brainsight(c)';
        NIRS.Dt.fir.stax.p{1} = staxp;
    else
        NIRS.Dt.fir.stax.n = 'Template LIOM'; % template
        [DirSPM,~,~] = fileparts(which('nirs10'));
        staxp = fullfile(DirSPM,'nirs10_templates','Brainsight(c).txt');
        NIRS.Dt.fir.stax.p{1} = staxp;
        % coordinates
        load(fullfile(dir_nt,'Hcoregistered.mat'));
        NIRS.Cf.H = Hcoregistered;
    end
    
    % Topo Data
    if isempty(job.subj(1,is).TopoData{:})
        % TopoData
        load(fullfile(dir_nt,'TopoData.mat'));
        NIRS.Dt.ana.rend = fullfile(dir_nt,'TopoData.mat');
        save(fullfile(NIRS.Dt.s.p,'TopoData.mat'),'rendered_MNI');
        disp('Inutile de faire la coregistration !');
        helmetdone =1;
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