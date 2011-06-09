function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

outNIRSmat ={};

% Study configuration
if isfield(job.study_cfg,'choose_path')
    mkdir(job.study_cfg.study_path.choose_path{:});
    study_p = job.study_cfg.study_path.choose_path{:};
else
    study_p = job.study_cfg.study_path.existing_study{:};
end

% Loop over all subjects
sN = size(job.subj,2);
for is=1:sN
    age = job.subj(1,is).age1;
    
    if isfield(job.subj(1,is),'subj_id')
        s_id = str2double(job.subj(1,is).subj_id);
        if s_id < 10,
            str0 = '00';
        else
            if s_id < 100
                str0 = '0';
            else
                str0 = '';
            end
        end
        sDtp = fullfile(study_p,['S' str0 int2str(s_id)]);
    end
    
    % Reinitialize NIRS matrix for each subject
    NIRS = [];
    NIRS.Dt.s.age = age;
    NIRS.Dt.s.p = sDtp;
    % on ecrase tout sujet deja existant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist(sDtp,'dir')==7, rmdir(sDtp,'s'); end
    mkdir(sDtp);
    mkdir(fullfile(sDtp,'T1'));
    mkdir(fullfile(sDtp,'fir'));
              
    % BOLD mask
    if ~isempty(job.subj(1,is).boldmask{:})
        NIRS.Cm.bold = job.subj(1,is).boldmask{:};
    end
    
    % Helmet
    if isfield(job.subj(1,is).helmet,'text_brainsight')% Reading subject-specific setup from BrainSight file
        staxp = job.subj(1,is).helmet.text_brainsight{:};
        NIRS.Dt.fir.stax.n = 'Brainsight(c)';
        NIRS.Dt.fir.stax.p{1} = job.subj(1,is).helmet.text_brainsight{:};
        
        if ~job.subj(1,is).allSD_autosave % if group analysis you may want to run all the subjects without having to confirm at each time...
            jobH.subj.sDtp = sDtp;
            jobH.subj.helmet.staxp = staxp;
            nirs_criugm_getHelmet(jobH); % get helmet configuration (S,D,P,Q) from Brainsight
            
            fig=gcf; %findall(0,'name','Get positions from Brainsight (%%%%%)');
            waitfor(fig,'BeingDeleted','On');
        else % all sources and all detectors are selected
            jobH.subj.sDtp = sDtp;
            jobH.subj.helmet.staxp = staxp;
            nirs_criugm_getHelmet_allSD_autosave(jobH);
        end
        
        % CB: NOUVELLE VERSION : LA partie HELMET DOIT enregistrer dans le
        % dossier du sujet une structure NIRS_Cf qui est lue ici puis effacee .
        NIRS_Cf = load(fullfile(sDtp, 'NIRS_Cf.mat'));
        NIRS.Cf = NIRS_Cf.NIRS.Cf;
        clear NIRS_Cf
        delete(fullfile(sDtp, 'NIRS_Cf.mat'));
        
    elseif isfield(job.subj(1,is).helmet,'T1_vitamins')
        NIRS.Dt.fir.stax.n = 'T1_vitamins';
        
    elseif isfield(job.subj(1,is).helmet,'no_helmet')
        NIRS.Dt.fir.stax.n = 'no_helmet';
        
    elseif  isfield(job.subj(1,is).helmet,'custom')
        NIRS.Dt.fir.stax.n = 'custom';
        % H : NIRS.Cf.H already done, just have to save it in the NIRS.mat
        load(fullfile(job.subj(1,is).helmet.custom{:},'Hcoregistered.mat'));
        NIRS.Cf.H = Hcoregistered;
        % Topo Data : nirs_run_coreg has already been executed to generate
        % once and for all the subjects the TopoData matrix.
        NIRS.Dt.ana.rend = fullfile(job.subj(1,is).helmet.custom{:},'TopoData.mat');
        % save template T1 as ana T1
        anatT1 = fullfile(fileparts(which('spm')),'templates','T1.nii');
        job.subj(1,is).anatT1 = {anatT1};
    end
    
    % Anatomical image
    if ~isempty(job.subj(1,is).anatT1{:})
        [ana,ana_nam,ana_ext] = fileparts(job.subj(1,is).anatT1{:});
        if ~strcmp(fullfile(sDtp,'T1'),ana) %%%%% a changer en ana, a terme !!!!!
            copyfile(job.subj(1,is).anatT1{:},fullfile(sDtp,'T1',[ana_nam ana_ext]));
        end
        NIRS.Dt.ana.T1 = job.subj(1,is).anatT1{:};
    end
    
    % nirs files
    if ~isempty(job.subj(1,is).nirs_files{1,1})
        
        for fi=1:size(job.subj(1,is).nirs_files,1)
            [dummy1,namef,extf] = fileparts(job.subj(1,is).nirs_files{fi,:});
            nirs_files{fi,:} = fullfile(sDtp,'fir',[namef extf]); 
            copyfile(job.subj(1,is).nirs_files{fi,:},nirs_files{fi,:});
        end
        % Read setup information from nirs file
        % System used for acquisition
        job1.system = job.subj(1,is).CWsystem;
        % Read only nirs file from first session
        job1.nirs_file = load(nirs_files{1,1},'-mat');
        job1.sDtp = sDtp;
        job1.coregType = NIRS.Dt.fir.stax.n;
        job1.NIRS = NIRS;
        % The function will update the NIRS matrix
        NIRS = nirs_criugm_readtechen(job1);% get C configuration from nirs files
        
        NU=[];% number of session EST ON BIEN SUR QUE C EST PAS SESS ??????????? NOTATIONS PAS CONSISTANTES
        if ~strcmp(nirs_files,''), NU = size(nirs_files,1); end
        % Loop over all sessions
        for iU=1:NU % # of data files
            fp = nirs_files(iU,1);
            %clear f
            %             f = load(fp{:},'-mat');
            
            NIRS.Dt.fir.pp(1).p{iU,1} = fp{:};
            % ON ne melange pas les inputs des codes !!!            NIRS.Dt.fir.pp(1).m{iU,1} = job.subj(1,is).baseline_method;
            NIRS.Dt.fir.pp(1).pre = 'readCriugm';
            NIRS.Dt.fir.pp(1).job = job;
            
            %Ignore parametric modulations - cf spm_run_fmri_design.m
            P.name = 'none';
            P.h    = 0;
            
            % Protocol
            if ~isempty(job.study_cfg.protocol{1,1}) && iU <= size(job.study_cfg.protocol,1)
                % Read "multiple conditions" file (.mat)
                load(job.study_cfg.protocol{iU,1},'-mat');
                for kk = 1:size(names, 2)
                    NIRS.Dt.fir.Sess(iU).U(kk).name = names(kk);
                    NIRS.Dt.fir.Sess(iU).U(kk).ons = onsets{kk};
                    NIRS.Dt.fir.Sess(iU).U(kk).dur = durations{kk};
                    NIRS.Dt.fir.Sess(iU).U(kk).P = P;
                end
            else
                NIRS.Dt.fir.Sess(iU).U.name = {};
                NIRS.Dt.fir.Sess(iU).U.ons = [];
                NIRS.Dt.fir.Sess(iU).U.dur = [];
                NIRS.Dt.fir.Sess(iU).U.P = P;
            end
        end
    end
    
    save(fullfile(sDtp, 'NIRS.mat'),'NIRS');
    outNIRSmat = [outNIRSmat; fullfile(sDtp,'NIRS.mat')];
end

out.NIRSmat = outNIRSmat;

end
