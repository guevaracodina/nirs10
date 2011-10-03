function out = nirs_run_criugm(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-11-05

outNIRSmat ={};
% Study configuration
study_p = job.study_cfg.study_path{:};
%     if isfield(job.study_cfg.study_path,'choose_path')
%         study_p = job.study_cfg.study_path.choose_path;
%         if ~exist(study_p,'dir'),
%             mkdir(study_p);
%         end
%     else
%         study_p = job.study_cfg.study_path.existing_study{:};
%     end

% Loop over all subjects
sN = size(job.subj,2);
for is=1:sN
    try
        if isfield(job.study_cfg.indvdata,'template_chosen')
            % save template T1 as anatT1
            anatT1 = fullfile(fileparts(which('spm')),'templates','T1.nii');
            job.subj(1,is).anatT1 = {anatT1};
            % on saute ensuite les segmentations et coregistrations
            template4all=1;
            NIRS.Cf.H.template4all =template4all;
        else
            template4all=0;
            NIRS.Cf.H.template4all = template4all;
        end
        
        age = job.subj(1,is).age1;
        
        if isfield(job.subj(1,is),'subj_id')
            %%% verifier qu on a bien un nombre...
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
        %if exist(sDtp,'dir')==7, rmdir(sDtp,'s'); end
        if ~exist(sDtp,'dir'), mkdir(sDtp); end
        T1dir = fullfile(sDtp,'T1');
        if ~exist(T1dir,'dir'), mkdir(T1dir); end
        firdir = fullfile(sDtp,'fir');
        if ~exist(firdir,'dir'), mkdir(firdir); end
        
        % BOLD mask
        if ~isempty(job.subj(1,is).boldmask{:})
            NIRS.Cm.bold = job.subj(1,is).boldmask{:};
        end
        
        % Helmet
        if isfield(job.subj(1,is).helmet,'text_brainsight') || template4all% Reading subject-specific setup from BrainSight file
            if template4all
                staxp = job.study_cfg.indvdata.template_chosen.text_brainsight{:};
                NIRS.Dt.fir.stax.n = 'Brainsight(c)';
                NIRS.Dt.fir.stax.nota = 'Brainsight template';
            else
                staxp = job.subj(1,is).helmet.text_brainsight{:};
                NIRS.Dt.fir.stax.n = 'Brainsight(c)';
            end
            NIRS.Dt.fir.stax.p{1} = staxp;
            
            try % try to get Helmet automatically
                jobH.subj.sDtp = sDtp;
                jobH.subj.helmet.staxp = staxp;
                nirs_criugm_getHelmet_autosave(jobH);
            catch % launch GUI for the user to select the right points
                jobH.subj.sDtp = sDtp;
                jobH.subj.helmet.staxp = staxp;
                nirs_criugm_getHelmet2(jobH); % get helmet configuration (S,D,P,Q) from Brainsight
                
                fig=gcf;
                waitfor(fig,'BeingDeleted','On');
            end
            
            % CB: NOUVELLE VERSION : LA partie HELMET DOIT enregistrer dans le
            % dossier du sujet une structure NIRS_Cf qui est lue ici puis effacee .
            NIRS_Cf = load(fullfile(sDtp, 'NIRS_Cf.mat'));
            NIRS.Cf = NIRS_Cf.NIRS.Cf;
            clear NIRS_Cf
            delete(fullfile(sDtp, 'NIRS_Cf.mat'));
            
        elseif isfield(job.subj(1,is).helmet,'T1_vitamins') && ~template4all
            NIRS.Dt.fir.stax.n = 'T1_vitamins';
            
        elseif isfield(job.subj(1,is).helmet,'no_helmet') && ~template4all
            NIRS.Dt.fir.stax.n = 'no_helmet';% no so much interest
            
        elseif isfield(job.subj(1,is).helmet,'helm_temp') &&...
                ~isempty(job.subj(1,is).helmet.helm_temp) && ~template4all%%%% je pense au4il fqut enlever ce ~template4all
            %%% on a l anatomique et un template pour le helmet...
            % il faut coller les deux dans l espace normalise et la
            % coregistration necessaire n est que partielle !!!!
            NIRS.Dt.fir.stax.n = 'Brainsight(c)';
            NIRS.Dt.fir.stax.nota = 'Helmet template';
            NIRS_ht = load(job.subj(1,is).helmet.helm_temp{:});
            NIRS.Cf.H.P.w = NIRS_ht.NIRS.Cf.H.P.w;
            NIRS.Cf.H.C = NIRS_ht.NIRS.Cf.H.C;
            NIRS.Cf.H.p = NIRS_ht.NIRS.Cf.H.p;
            NIRS.Cf.H.n = NIRS_ht.NIRS.Cf.H.n;
            
            f = load(job.subj(1,is).nirs_files{:},'-mat');
            NIRS.Cf.H.S.n = f.SD.SrcNam;
            NIRS.Cf.H.S.N = size(f.SD.SrcPos,1);
            NIRS.Cf.H.D.n = f.SD.SrcNam;
            NIRS.Cf.H.D.N = size(f.SD.DetPos,1);
            NIRS.Cf.H.P.N = NIRS.Cf.H.S.N + NIRS.Cf.H.D.N;
            
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.NIRSmat = {fullfile(study_p,'NIRS4all.mat')};
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.DelPreviousData = 0;
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.anatT1 = {[fullfile(fileparts(which('spm')),'templates','T1.nii') ',1']};
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.segT1_4fit = {[fullfile(fileparts(which('spm')),'toolbox','nirs10','nirs10_templates','00044_segmented_T1.nii') ',1']};
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.anatT1_template = {[fullfile(fileparts(which('spm')),'templates','T1.nii') ',1']};
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.fid_in_subject_MNI = 0;
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.nasion_wMNI = [0 84 -48];
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.AL_wMNI = [-83 -19 -38];
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.AR_wMNI = [83 -19 -38];
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.GenDataTopo = 1;
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.render_choice.render_template = struct([]);
            %             matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg3.View6Projections = 0;
            %             spm_jobman('run_nogui',matlabbatch);
        end
        
        % Anatomical image
        if ~isempty(job.subj(1,is).anatT1{:})
            [ana,ana_nam,ana_ext] = fileparts(job.subj(1,is).anatT1{:});
            ana_ext = ana_ext(1:4);
            if ~strcmp(fullfile(sDtp,'T1'),ana) %%%%% a changer en ana, a terme !!!!!
                copyfile(fullfile(ana,[ana_nam ana_ext]),fullfile(sDtp,'T1',[ana_nam ana_ext]));
                if strcmp(ana_ext,'.img')
                    copyfile(fullfile(ana,[ana_nam '.hdr']),fullfile(sDtp,'T1',[ana_nam '.hdr']));
                end
            end
            NIRS.Dt.ana.T1 = fullfile(sDtp,'T1',[ana_nam ana_ext]);
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
            Nons=[];
            if ~strcmp(job.subj(1,is).protocol,''), Nons = size(job.subj(1,is).protocol,1); end
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
                
                
                %%%%%
                %%% MODIFIER POUR INTEGRER LES ONSETS DIFFEREMMENT %%%
                % Protocol
                if iU<=Nons && ~isempty(job.subj(1,is).protocol{iU,:})
                    % Read "multiple conditions" file (.mat)
                    clear names onsets durations;
                    load(job.subj(1,is).protocol{iU,:},'-mat');
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
        
        %%% on lance la coregistration pour le premier sujet uniquement
        %%% avec les donnees du sujet 0
        if template4all && is==1
            save(fullfile(study_p, 'NIRS4all.mat'),'NIRS');
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.NIRSmat = {fullfile(study_p,'NIRS4all.mat')};
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.DelPreviousData = 0;
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.NewDirCopyNIRS.CreateNIRSCopy_false = struct([]);
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.anatT1 = NIRS.Dt.ana.T1;
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.segT1_4fit = {[fullfile(fileparts(which('spm')),'toolbox','nirs10','nirs10_templates','00044_segmented_T1.nii') ',1']};
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.anatT1_template = {[fullfile(fileparts(which('spm')),'templates','T1.nii') ',1']};
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.fid_in_subject_MNI = 0;
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.nasion_wMNI = [0 84 -48];
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.AL_wMNI = [-83 -19 -38];
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.AR_wMNI = [83 -19 -38];
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.GenDataTopo = 1;
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.render_choice.render_template = struct([]);
            matlabbatch{1}.spm.tools.nirs10.coregNIRS.coreg2.View6Projections = 0;
            spm_jobman('run_nogui',matlabbatch);
        end
        
        if template4all
            NIRS4all = load(fullfile(study_p, 'NIRS4all.mat'));
            NIRS.Cf.H = NIRS4all.NIRS.Cf.H;
            NIRS.Dt.ana = NIRS4all.NIRS.Dt.ana;
            %             NIRS.Dt.pro.errValofCoreg_mm2 = NIRS4all.NIRS.Dt.pro.errValofCoreg_mm2;
            copyfile(fullfile(study_p, 'TopoData.mat'),fullfile(sDtp, 'TopoData.mat'));
        end
        
        save(fullfile(sDtp, 'NIRS.mat'),'NIRS');
        outNIRSmat = [outNIRSmat; fullfile(sDtp,'NIRS.mat')];
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
    end
end
out.NIRSmat = outNIRSmat;
end
