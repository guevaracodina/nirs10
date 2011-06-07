function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

outNIRSmat ={};

% Loop over all subjects
sN = size(job.subj,2);
for is=1:sN
    age = job.subj(1,is).age1;
    sDtp = job.subj(1,is).subj_path{:};
    
    % Reinitialize NIRS matrix for each subject
    NIRS = [];
    NIRS.Dt.s.age = age;
    NIRS.Dt.s.p = sDtp;
    
    % Anatomical image
    if ~isempty(job.subj(1,is).anatT1{:})
        NIRS.Dt.ana.T1 = job.subj(1,is).anatT1{:};
    end
    
    % BOLD mask
    if ~isempty(job.subj(1,is).boldmask{:})
        NIRS.Cm.bold = job.subj(1,is).boldmask{:};
    end
    
    % Helmet
    if isfield(job.subj(1,is).helmet,'text_brainsight')
        staxp = job.subj(1,is).helmet.text_brainsight{:};
        NIRS.Dt.fir.stax.n = 'Brainsight(c)';
        NIRS.Dt.fir.stax.p{1} = job.subj(1,is).helmet.text_brainsight{:};
        helmetdone = 0;
        if ~isempty(strfind(job.subj(1,is).helmet.text_brainsight{:},'template'))
            %%% CB : etude Said.... a mettre coherent /////
            % coordinates
            load(fullfile(fileparts(which('nirs10')),'nirs10_templates','Hcoregistered.mat'));
            NIRS.Cf.H = Hcoregistered;
            
            % Topo Data
            % CLEMENT: Vérifier Hcoregistered et TopoData (lien?)
            if ~isempty(job.subj(1,is).TopoData{:})
                helmetdone = 1;
                % If nirs_run_coreg has already been executed to generate once and
                % for all the TopoData matrix.
            else
                helmetdone = 0;
            end
            
        else% Reading subject-specific setup from BrainSight file
            if ~helmetdone
                jobH.subj.sDtp = sDtp;
                jobH.subj.helmet.staxp = staxp;
                outH = nirs_criugm_getHelmet(jobH); % get helmet configuration (S,D,P,Q) from Brainsight
            end
            fig=gcf; %findall(0,'name','Get positions from Brainsight (clbon)');
            waitfor(fig,'BeingDeleted','On');
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
    end
    
    % nirs files
    if ~isempty(job.subj(1,is).nirs_files{1,1})
        % Read setup information from nirs file
        % System used for acquisition
        job1.system = job.subj(1,is).CWsystem;
        % Read only nirs file from first session
        job1.nirs_file = load(job.subj(1,is).nirs_files{1,1},'-mat');
        job1.sDtp = sDtp;
        job1.coregType = NIRS.Dt.fir.stax.n;
        job1.NIRS = NIRS;
        % The function will update the NIRS matrix
        NIRS = nirs_criugm_readtechen(job1);% get C configuration from nirs files
        
        NU=[];% number of session EST ON BIEN SUR QUE C EST PAS SESS ??????????? NOTATIONS PAS CONSISTANTES
        if ~strcmp(job.subj(1,is).nirs_files,''), NU = size(job.subj(1,is).nirs_files,1); end
        % Loop over all sessions
        for iU=1:NU % # of data files
            fp = job.subj(1,is).nirs_files(iU,1);
            %clear f
            f = load(fp{:},'-mat');
            
            NIRS.Dt.fir.pp(1).p{iU,1} = fp{:};
            NIRS.Dt.fir.pp(1).m{iU,1} = job.subj(1,is).baseline_method;
            NIRS.Dt.fir.pp(1).pre = 'readCriugm';
            NIRS.Dt.fir.pp(1).job = job;
            
            %Ignore parametric modulations - cf spm_run_fmri_design.m
            P.name = 'none';
            P.h    = 0;
            
            % Protocol
            if ~isempty(job.subj(1,is).protocol{1,1}) && iU <= size(job.subj(1,is).protocol,1)
                % Read "multiple conditions" file (.mat)
                load(job.subj(1,is).protocol{iU,1},'-mat');
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
