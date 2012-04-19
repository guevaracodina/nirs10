function out = nirs_run_liom_HDM(job) 
% User interface for hemodynamic model estimation
% Based on SPM HDM


nameHDM = job.nameHDM;
Model_Choice = job.Model_Choice;
%subjects = job.which_subjects; % path to SPM.mat files

TE = job.echo_time;    
Stimuli = job.which_condition;
sessions = job.which_session;
ROIs = job.whichROI;

dp_start = job.dp_start; %points to remove due to filter set-up
dp_end = job.dp_end; %points to remove due to filter set-up

removeWhitening = job.removeWhitening;

save_figures = job.save_figures;
generate_figures = job.generate_figures;


%Modalities:
%1: BOLD
%2: BOLD&flow from stat maps (not coded up)
%3: flow (not coded up)
%4: HbO HbR (not implemented in this version)
%5: BOLD&flow v2 from data

%Step 1: checking which modality to run
if isfield(job.xSPM_Modalities,'xSPM_BOLD')
    modal = 1;
    subjects_bold = job.xSPM_Modalities.xSPM_BOLD.which_subjects_bold;
    %sBOLD = job.xSPM_Modalities.xSPM_BOLD;
else
    if isfield(job.xSPM_Modalities,'xSPM_BOLD_ASL')
        modal = 2;
        subjects_bold = job.xSPM_Modalities.xSPM_BOLD_ASL.which_subjects_bold;                
        subjects_flow = job.xSPM_Modalities.xSPM_BOLD_ASL.which_subjects_flow;
        %sBOLD = job.xSPM_Modalities.xSPM_BOLD_ASL;
        %currently not used
        %sASL = job.xSPM_Modalities.xSPM_BOLD_ASL;
        %subjectsASL = sASL.which_subjects_ASL;
    else
        if isfield(job.xSPM_Modalities,'xSPM_ASL')
            modal = 3;
            subjects_flow = job.xSPM_Modalities.xSPM_ASL.which_subjects_flow;
            %sASL = job.xSPM_Modalities.xSPM_ASL;
        else
            if isfield(job.xSPM_Modalities,'xSPM_BOLD_ASL_V2')
                modal = 5;
                subjects_bold = job.xSPM_Modalities.xSPM_BOLD_ASL_V2.which_subjects_bold;
                subjects_flow = job.xSPM_Modalities.xSPM_BOLD_ASL_V2.which_subjects_flow;
                %sBOLD = job.xSPM_Modalities.xSPM_BOLD_ASL_V2;
            end
        end
    end
end

% Number of subjects
switch modal
    case {1,2,5}
        nSubj = length(subjects_bold);
    case 3
        nSubj = length(subjects_flow);
end

switch modal
    case {2,4,5}
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        plot_algebraic_CMRO2 = 1; %Boolean
    case {1,3}
        plot_algebraic_CMRO2 = 0;
end

if isfield(job,'priorFile')
    priorFile = job.priorFile{1}; % spm_hdm_prior.m filename (to set priors for model inversion)
else
    priorFile = '';
end

if isfield(job.simuOn,'simuYes')
    S.simuOn = 1;
    S.simuS     = job.simuOn.simuYes.simuS; %Stimuli types to include
    S.simuIt    = job.simuOn.simuYes.simuIt; %Number of random iterations
    S.simuA     = job.simuOn.simuYes.simuA; %Signal amplitude, as % of BOLD signal
    S.simuP     = job.simuOn.simuYes.simuP; %Parameters to vary
    if isfield(job.simuOn.simuYes.simuParamDistr,'distr_uniform')
        S.simuPriorDistr = 1;
        S.simuPrior = job.simuOn.simuYes.simuParamDistr.distr_uniform.simuPrior; %Priors to use
        S.simuR     = job.simuOn.simuYes.simuParamDistr.distr_uniform.simuR; %Range to sample
    end
    if isfield(job.simuOn.simuYes.simuParamDistr,'distr_bimodal')
        S.simuPriorDistr = 2;
        S.simuPrior1 = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuMean1; %Priors to use
        S.simuR1     = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuR1; %Range to sample
        S.simuDiffPrior = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuMean21; %Priors to use
        S.simuR2     = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuR2; %Range to sample
    end

    S.simuUpsample = job.simuOn.simuYes.simuUpsample; %Upsampling factor on data
    try % for back-compatibility
        S.simuInterp = job.simuOn.simuYes.simuInterp; %Interpolation on simulated data (after decimating)
    catch
        S.simuInterp = 1;
    end
    simuNoise1 = job.simuOn.simuYes.simuNoise; %Yes to include background noise based on restscans
    if isfield(simuNoise1,'noiseYes')
        S.simuNoise = 1;
        switch modal
            case 1
                restscans = job.simuOn.simuYes.simuNoise.noiseYes.restscans; %Rest scans to add signal to
            case {2,5}
                restscans_BOLD = job.simuOn.simuYes.simuNoise.noiseYes.restscans_BOLD;
                restscans_ASL = job.simuOn.simuYes.simuNoise.noiseYes.restscans_ASL;
            case 3
                restscans_ASL = job.simuOn.simuYes.simuNoise.noiseYes.restscans_ASL;
        end
    else
        S.simuNoise = 0;
    end
    
else
    S.simuOn = 0;
    
end

%EM parameters
EM.spm_integrator = job.EM_parameters.spm_integrator;
EM.Niterations = job.EM_parameters.Niterations;
EM.dFcriterion = job.EM_parameters.dFcriterion;
EM.LogAscentRate = job.EM_parameters.LogAscentRate;
EM.MaxLogAscentRate = job.EM_parameters.MaxLogAscentRate;
EM.Mstep_iterations = job.EM_parameters.Mstep_iterations;

%---------------------------------------------------------------------------
try
    
    
    % inputs
    %==========================================================================
    try 
        spmpath_tmp = subjects_bold{1};
    catch
        spmpath_tmp = subjects_flow{1};
    end
        
    tmp = strfind(spmpath_tmp,filesep);
    rootDir = spmpath_tmp(1:tmp(end-2));
    
    if ~isempty(nameHDM)
        figDir = fullfile(rootDir,['HDM_' nameHDM]);
    else
        figDir = fullfile(rootDir,'HDM');
    end
    if save_figures
        if ~exist(figDir,'dir'),mkdir(figDir);end
        HDMfile = fullfile(figDir,'HDM.mat');
    else
        HDMfile = fullfile(rootDir,'HDM.mat');
    end
    
    %loop over each session
    for s1=1:length(sessions)
        %current session
        cs1 = sessions(s1);
        %loop over each ROI
        for r1=1:length(ROIs)
            if isempty(ROIs(r1).nameROI)
                ROIs(r1).nameROI = int2str(r1);
            end
            %loop over each subject
            for SubjIdx = 1:nSubj
                fullfigDir = fullfile(figDir,['S' int2str(cs1) '_' ROIs(r1).nameROI]);

                try
                    
                    % Data files
                    clear SPM Y
                    if exist('subjects_bold')
                        fBOLD = fullfile(subjects_bold{SubjIdx},'SPM.mat');
                        fBOLD_old = fBOLD;
                    end
                    if exist('subjects_flow')
                        fASL = fullfile(subjects_flow{SubjIdx},'SPM.mat');
                    end
%                     
%                     % TO IMPLEMENT: REPLACE WITH A SEPARATE INPUT "FLOW SPM
%                     % FOLDERS"
%                     switch modal
%                         case {2,3,5}
%                             % Very inelegant HARD-CODED FOR MDEIE-P2
%                             tmpidx3 = strfind(fBOLD,'UR3');
%                             tmpidx1 = strfind(fBOLD,'UR1');
%                             if ~isempty(tmpidx1)
%                                 fASL = fBOLD;
%                             elseif isempty(tmpidx1) && ~isempty(tmpidx3)
%                                 fASL = fBOLD;
%                                 fASL(tmpidx3:tmpidx3+2) = 'UR1';
%                             else
%                                 tmpidxb = strfind(fBOLD,'bold');
%                                 tmpidxf = strfind(fBOLD,'flow');
%                                 if ~isempty(tmpidxf)
%                                     fASL = fBOLD;
%                                 elseif isempty(tmpidxf) && ~isempty(tmpidxb)
%                                     fASL = fBOLD;
%                                     fASL(tmpidxb:tmpidxb+3) = 'flow';
%                                 end
%                             end
%                     end
                    

                    % And for BOLD    
                    % CHECK THIS...??
                    if removeWhitening
                        % For flow...
                        switch modal
                            case {2,3,5}
                                load(fASL)
                                %quick fix to use only one session
                                 nSess = size(SPM.Sess,2);
                                 SPM.xX.W = speye(size(SPM.xX.W,1)./nSess);
                                 SPM.xX.K = SPM.xX.K(cs1);
                                 %Remove HPF entirely for now
                                 SPM.xX.K.X0 = zeros(size(SPM.xX.K.X0));
                                 nScans = size(SPM.xY.VY,1)/nSess;
                                 SPM.xY.VY = SPM.xY.VY((1:nScans) + (nScans*(cs1-1)) );
                                 
                                 %save it as a temporary file
                                [dir0 fil0 ext0] = fileparts(fASL);
                                tmpdir = fullfile(dir0,'tmp_noWhitening');
                                if ~exist(tmpdir,'dir')
                                    mkdir(tmpdir)
                                    fASL = fullfile(tmpdir,[fil0 ext0]);
                                    save(fASL,'SPM');
                                end
                        end
                        
                        switch modal
                            case {1,2,5}
                                load(fBOLD);
                                %quick fix to use only one session
                                 nSess = size(SPM.Sess,2);
                                 SPM.xX.W = speye(size(SPM.xX.W,1)./nSess);
                                 SPM.xX.K = SPM.xX.K(cs1);
                                 %Remove HPF entirely for now
                                 SPM.xX.K.X0 = zeros(size(SPM.xX.K.X0));
                                 nScans = size(SPM.xY.VY,1)/nSess;
                                 SPM.xY.VY = SPM.xY.VY((1:nScans) + (nScans*(cs1-1)) );
                        end
                        %save it as a temporary file
                        [dir0 fil0 ext0] = fileparts(fBOLD);
                        tmpdir = fullfile(dir0,'tmp_noWhitening');
                        if ~exist(tmpdir,'dir')
                            mkdir(tmpdir)
                            fBOLD = fullfile(tmpdir,[fil0 ext0]);
                            save(fBOLD,'SPM');
                        end
                        
                    else % if not removing filter
                        switch modal
                            case 3
                                load(fASL);
                            case {1,2,5}
                                load(fBOLD);
                        end
                        
                    end
                    
                    [Sess s] = get_session(SPM,cs1);
                    U = get_causes(Sess,Stimuli);
                    
                    %if S.simuOn, replace the data with the rest data
                    if S.simuOn 
                        %                        if S.simuNoise
                        switch modal
                            case 1
                                if S.simuNoise
                                    SPM.xY.P     = char(restscans{:});
                                    SPM.xY.VY = spm_vol(SPM.xY.P);
                                end
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fBOLD);
                                tmpdir = fullfile(dir00,'tmp_simu');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                fBOLD_old = fBOLD;
                                fBOLD = fullfile(tmpdir,[fil00 ext00]);
                                save(fBOLD,'SPM');
                                
                            case 3
                                 if S.simuNoise
                                    SPM.xY.P     = char(restscans_ASL{:});
                                    SPM.xY.VY = spm_vol(SPM.xY.P);
                                end
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fASL);
                                tmpdir = fullfile(dir00,'tmp_simu_ASL');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                %fASL_old = fBOLD;
                                fASL = fullfile(tmpdir,[fil00 ext00]);
                                save(fASL,'SPM');
                                
                            case 5
                                if S.simuNoise
                                    SPM.xY.P     = char(restscans_BOLD{:});
                                    SPM.xY.VY = spm_vol(SPM.xY.P);
                                end
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fBOLD);
                                tmpdir = fullfile(dir00,'tmp_simu_BOLD');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                fBOLD_old = fBOLD;
                                fBOLD = fullfile(tmpdir,[fil00 ext00]);
                                save(fBOLD,'SPM');
                                
                                if S.simuNoise
                                    SPM.xY.P     = char(restscans_ASL{:});
                                    SPM.xY.VY = spm_vol(SPM.xY.P);
                                end
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fBOLD_old);
                                tmpdir = fullfile(dir00,'tmp_simu_ASL');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                %fASL_old = fBOLD;
                                fASL = fullfile(tmpdir,[fil00 ext00]);
                                save(fASL,'SPM');
                                %                            end
                        end
                    end

                    
                    
                    % Extract VOI data
                    switch modal
                        case 1
                            %extract VOI for BOLD
                            VOI = get_VOI(fBOLD,ROIs(r1),cs1);                 
                            y1 = VOI.Y;                  
                            y1 = y1((1+dp_start):end-dp_end);                   
                            Y.y(:,1) = y1./mean(y1);
                        case {2,5}
                            %extract VOI for BOLD
                            VOI = get_VOI(fBOLD,ROIs(r1),cs1);                 
                            y1 = VOI.Y;                  
                            y1 = y1((1+dp_start):end-dp_end);                   
                            Y.y(:,1) = y1./mean(y1);
                            %extract VOI for ASL
                            VOI_ASL = get_VOI(fASL,ROIs(r1),cs1);
                            y2    = VOI_ASL.Y((1+dp_start):end-dp_end); %y/100;
                            %yf2 = ButterHPF(1/SPM.xY.RT,1/128,3,y2);
                            Y.y(:,2) = y2/mean(y2);
                        case 3
                            VOI = get_VOI(fASL,ROIs(r1),cs1);
                            y2    = VOI.Y((1+dp_start):end-dp_end); %y/100;
                            %yf2 = ButterHPF(1/SPM.xY.RT,1/128,3,y2);
                            Y.y(:,1) = y2/mean(y2);
                    end
                    
           
                    %-place response and confounds in response structure
                    %--------------------------------------------------------------------------
                    switch modal
                        case {1,2,3,5}
                            %y      = VOI.xY.u;
                            VOI.xY.RT = SPM.xY.RT;
                            Y.dt   = VOI.xY.RT; %SPM.xY.RT;
                            Y.X0   = VOI.xY.X0((1+dp_start):(end-dp_end),:);
                    end
                    %remove start and end points from stimuli too to match the data
                    U.u = U.u( (1+dp_start*(Y.dt/U.dt)):(end-dp_end*(Y.dt/U.dt)),:);
                    
                    if generate_figures || save_figures
                        spm('Pointer','Watch')
                        spm('FigName','Estimation in progress');
                    end                                       
                    
                    % Priors
                    %--------------------------------------------------------------------------
                    pE = []; pC = [];
                    %number of inputs
                    m       = size(U.u,2);
                    [pE,pC] = set_priors(priorFile,m);        
                    
                    % Generate random parameter distributions
                    if S.simuOn
                        pA = generate_random_param(S,pE);
                    end
                    
                    % Set model
                    %--------------------------------------------------------------------------
                    M.modal = modal;
                    M.Model_Choice = Model_Choice;                    
                    M.N     = 64;
                    M.dt    = 24/M.N; %24/M.N;
                    M.TE    = TE;
                    M.pE    = pE;
                    M.pC    = pC;
                    M.m     = m;
                    M.EM    = EM;
                    if generate_figures || save_figures
                        M.nograph = 0;
                    else
                        M.nograph = 1;
                    end
                    
                    M = set_model(M);
   
                    if ~S.simuOn
                        S.simuIt = 1;
                    else
                        Y0 = Y;
                        if S.simuS == 0 % Stimuli type(s) to includes
                            S.simuS = 1:size(U.u,2);
                        end
                    end
                    
                    %--------------------------------------------------------------------------
                    warning('off','MATLAB:nearlySingularMatrix');
                    warning('off','MATLAB:singularMatrix');
                    
                    for it1=1:S.simuIt
                        try
                            if S.simuOn  
                                Y = Y0;      
                                [Y,P] = simulate_data(pA,modal,M,U,Y,S,it1);
                            end


                            % INVERT MODEL %
                            % nonlinear system identification
                            %--------------------------------------------------------------------------
                            hfigevolution = []; hfigevolution2 = [];
                            if save_figures
                                fullfigDir1 = fullfigDir;
                            else
                                fullfigDir1 = '';
                            end
                            
                            % Plot data for debugging 
                            tt = (0:size(Y.y,1)-1)*(Y.dt);
                            figure('Units','normalized','Position',[0.1, 0.35, 0.6, 0.5])
                            plot(tt,Y.y(:,1)*100,'.-b')
                            %ylim([-0.3 1.8])
                            xlim([0 200])
                            xlabel('Time (s)')
                            switch modal
                                case 1
                                    ylabel('BOLD (%)')
                                    ylim(100*[min(Y.y(:)) max(Y.y(:))]);
                                case 3
                                    ylabel('Flow (%)')
                                    ylim(100*[min(Y.y(:)) max(Y.y(:))]);
                                case {2,5}
                                    hold on, plot(tt,Y.y(:,2)*100,'.-r')
                                    %ylim([-3 18])
                                    ylabel('ASL (BOLD & flow) (%)')
                                    ylim(100*[min(Y.y(:)) max(Y.y(:))]);
                            end
                            try 
                                title(['Upsampling: ' num2str(S.simuUpsample) ' ; Interpolation: ' num2str(S.simuInterp)])
                            end
                            
                            if ~exist(fullfigDir1,'dir'), mkdir(fullfigDir1); end
                            saveas(gcf,fullfile(fullfigDir1,['Yy_' gen_num_str(it1,3)]),'fig');
                            print(gcf, '-dtiffn', fullfile(fullfigDir1,['Yy_' gen_num_str(it1,3)]));
                            close(gcf);
                                                    
                            if S.simuOn
                                SubjIdx0 = SubjIdx; %should not be used
                                SubjIdx = it1; %changed inside a for loop but it is restored later
                            end
                            %[Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F,hfigevolution,hfigevolution2] = nirs_nlsi(M,U,Y);
                            [Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F] = nirs_nlsi(M,U,Y,fullfigDir1,SubjIdx);
                            %for simulations: store results in place of subjects

                            if S.simuOn
                                HDM{SubjIdx,r1}{s1}.EpS = P;
                            end
                            %Store results
                            if it1 == 1
                                %information independent of simulation iteration number
                                HDM{SubjIdx,r1}{s1}.M = M;
                                HDM{SubjIdx,r1}{s1}.pE = pE;
                                HDM{SubjIdx,r1}{s1}.job = job;                            
                            end
                            HDM{SubjIdx,r1}{s1}.Ep = Ep;
                            HDM{SubjIdx,r1}{s1}.Cp = Cp;
                            HDM{SubjIdx,r1}{s1}.F = F;
                            HDM{SubjIdx,r1}{s1}.name = ROIs(r1).nameROI;
                            HDM{SubjIdx,r1}{s1}.subj = fBOLD; 
                            HDM{SubjIdx,r1}{s1}.session = cs1;
                            HDM{SubjIdx,r1}{s1}.Y = Y;
                            HDM{SubjIdx,r1}{s1}.U = U;
                            
                            % Volterra kernels of states
                            %--------------------------------------------------------------------------
                            [dummy,H1] = spm_kernels(M0,M1,M.N,M.dt);
                            HDM{SubjIdx,r1}{s1}.H1 = H1;
                            %save - overwriting previous
                            save(HDMfile,'HDM');
                            
                            %-display results
                            %==========================================================================
                            if generate_figures || save_figures
                                display_results(fullfigDir,Model_Choice,pE,Ep,Cp,U,m,M,H1,K1,K2,...
                                    modal,plot_algebraic_CMRO2,save_figures,SubjIdx,hfigevolution,hfigevolution2);
                            end

                            if S.simuOn 
                                SubjIdx = SubjIdx0; % not used...
                            end
                        
                        catch exception
                            disp(['Model estimation or Simu #' int2str(it1) ' failed']);
                            disp(exception.identifier);
                            disp(exception.stack(1));
                        end
                        
                    end
                    
                    out = [];
                    
                catch exception
                    disp(exception.identifier);
                    disp(exception.stack(1));
                    out = [0];
                end
            end
        end
    end
    
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end

function [Sess s] = get_session(SPM,sessions)
% which session?
%---------------------------------------------------------------------------
s    = length(SPM.Sess);
if s > 1
    if sessions <= s
        if sessions
            s = sessions;
        end
    end
end
Sess = SPM.Sess(s);

end

function  U = get_causes(Sess,Stimuli)
% 'causes' or imputs U
%---------------------------------------------------------------------------
u_list = Stimuli; %1; %LFCP
u = length(u_list);
flagu = 1;

U.dt = Sess.U(1).dt;

if u == 1 && length(Sess.U(1).name) == 1
    U.name = Sess.U(1).name;
    U.u    = Sess.U(1).u(33:end,1);
else
    U.name = {};
    U.u    = [];
    for  i = 1:u
        if flagu
            U.u             = [U.u Sess.U(u_list(i)).u(33:end,1)];
            U.name{end + 1} = Sess.U(u_list(i)).name{1};
        else
            for  j = 1:length(Sess.U(i).name)
                str   = ['include ' Sess.U(i).name{j} '?'];
                if spm_input(str,2,'y/n',[1 0])
                    U.u             = [U.u Sess.U(i).u(33:end,j)];
                    U.name{end + 1} = Sess.U(i).name{j};
                end
            end
        end
    end
end
end

function VOI = get_VOI(SPMfile,ROI,s)
clear matlabbatch
matlabbatch{1}.spm.util.voi.spmmat = {SPMfile};
matlabbatch{1}.spm.util.voi.adjust = 0;
matlabbatch{1}.spm.util.voi.session = s;
matlabbatch{1}.spm.util.voi.name = ROI.nameROI;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = ROI.coordinateROI;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = ROI.radiusROI;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';
spm_jobman('run',matlabbatch);
%load VOI
[dir1 dummy] = fileparts(SPMfile);
load(fullfile(dir1,['VOI_' ROI.nameROI '_' int2str(s) '.mat']));
VOI.Y = Y;
VOI.xY = xY;
end

function [pE,pC] = set_priors(priorFile,m)
if ~isempty(priorFile)
    dirnow = pwd;
    [pdir1 pfil1 pext1] = fileparts(priorFile);
    cd(pdir1);
    eval(['[pE,pC] = ' pfil1 '(m,3);']);
    cd(dirnow);
end
if isempty(pE) || isempty(pC)
    switch Model
        case 0 %Buxton-Friston
            %[pE,pC] = spm_hdm_priors_YO(m);%PP pour Michèle.
            [pE,pC] = spm_hdm_priors(m,3);
        case 1 %Zheng-Mayhew
            [pE,pC] = nirs_hdm_priors_ZM(m,5);   %MODIFIER
        case 2 %Huppert1
            [pE,pC] = nirs_hdm_priors_Huppert1(m,3);
        otherwise
            [pE,pC] = spm_hdm_priors(m,3);
    end
end
end

function pA = generate_random_param(S,pE)
                        
    pA = repmat(pE',[S.simuIt 1]); %zeros(S.simuIt,length(pE));
    ct = 0;
    for pE1=1:length(pE)
        if any(S.simuP == 0) || any(pE1 == S.simuP)
            ct = ct+1;
            % Initialize the stream
            mtstream = RandStream('mt19937ar','Seed',pE1);
            RandStream.setDefaultStream(mtstream);
            % Generate the random numbers

            switch S.simuPriorDistr                                 
                case 1 % Uniform distribution
                    for it1=1:S.simuIt % for each simulation
                        % User-specified priors
                        if ~isempty(S.simuPrior)
                            pA(it1,pE1) = S.simuPrior(ct);
                        end
                        tpA = pA(it1,pE1);
                        % User-specified range of simulated
                        % values for parameters
                        if length(S.simuR) == 1
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR/100),tpA*(1+S.simuR/100));
                        else
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR(ct)/100),tpA*(1+S.simuR(ct)/100));
                        end
                    end

                case 2 % 2 Gaussians
                    for it1=1:round(S.simuIt/2) % for 1st half of simulations                                                                                      
                        % User-specified priors
                        if ~isempty(S.simuPrior1)
                            pA(it1,pE1) = S.simuPrior1(ct);
                        end
                        tpA = pA(it1,pE1);
                        % User-specified range of simulated
                        % values for parameters
                        if length(S.simuR1) == 1
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR1/100),tpA*(1+S.simuR1/100));
                            %normrnd(tpA,S.simuR1/100);
                        else
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR1(ct)/100),tpA*(1+S.simuR1(ct)/100));
                            %pA(it1,pE1) = normrnd(tpA,S.simuR1(ct)/100);
                        end
                    end

                    for it1=round(S.simuIt/2)+1:S.simuIt % for 2nd half of simulations                                                                                      
                        % User-specified priors
                        if ~isempty(S.simuDiffPrior)
                            pA(it1,pE1) = pA(it1-round(S.simuIt/2),pE1)*(1+S.simuDiffPrior(ct)/100);
                        end
                        tpA = pA(it1,pE1);
                        % User-specified range of simulated
                        % values for parameters
                        if length(S.simuR2) == 1
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR2/100),tpA*(1+S.simuR2/100));
                            %pA(it1,pE1) = normrnd(tpA,S.simuR2);
                        else
                            pA(it1,pE1) = unifrnd(tpA*(1-S.simuR2(ct)/100),tpA*(1+S.simuR2(ct)/100));
                            %pA(it1,pE1) = normrnd(tpA,S.simuR2(ct));
                        end
                    end

            end                                    
        end
    end

end

function M = set_model(M)
    switch M.Model_Choice
        case 0 %Buxton-Friston
            M.f     = 'spm_fx_hdm';
            M.g     = 'nirs_gx_hdm';
            M.x     = [0 0 0 0]';
            M.n     = 4;
        case 1 %Zheng-Mayhew
            M.f     = 'nirs_fx_hdm_ZM_BOLD';  %MODIFIER
            M.g     = 'nirs_gx_hdm_ZM_BOLD';  %MODIFIER
            M.x     = [0 0 0 0 0]';      %MODIFIER
            M.n     = 5; %MODIFIER
        case 2 %Huppert1
            M.f     = 'nirs_fx_hdm_Huppert1';  %MODIFIER
            M.g     = 'nirs_gx_hdm';  %MODIFIER
            M.x     = [0 0 0 0 0 0 0 0]';      %MODIFIER
            M.n     = 8; %MODIFIER
        case 3 % Buxton-Friston, part 1 (neural->flow, linear, eq. 1-2)
            M.f     = 'nirs_fx_hdm_part1';
            M.g     = 'nirs_gx_hdm';
            M.x     = [0 0]';
            M.n     = 2; % What is that?? number of state equations/variables??
        case 4 % Buxton-Friston, part 2 (flow->bold, non-linear, eq. 3-4)
            M.f     = 'spm_fx_hdm';%'spm_fx_hdm_part2';
            M.g     = 'nirs_gx_hdm';%'nirs_gx_hdm_part2';
            M.x     = [0 0 0 0]';
            M.n     = 4; % What is that?? umber of state equations/variables??    
        otherwise %Buxton-Friston
            M.f     = 'spm_fx_hdm';
            M.g     = 'nirs_gx_hdm';
            M.x     = [0 0 0 0]';
            M.n     = 4;
    end
    % NB: resting value/expansion point of x(1) is -Inf in log space; this is
    % taken into account in spm_fx_hdm.

    switch M.modal
        case {1,3}
            M.l = 1; %BOLD or flow only
        case 2
            M.l = 2; % BOLD & flow
        case 4 % HbO, HbR
            switch Model_Choice
                case 3 % BF part 1 (flow)
                    M.l = 1;
                case {0,1,4} %Buxton, Zheng-Mayhew - HbO, HbR
                    M.l = 2;
                case 2 %Huppert1 - HbO, HbR, HbT??
                    M.l = 3;
            end
        case 5
            M.l = 2; %BOLD & flow
    end
end

function [Y,P] = simulate_data(pA,modal,M,U,Y,S,it1)                                
                                
%Buxton-Friston balloon and other models
P = pA(it1,:)';
%efficacies
P(end-size(U.u,2)+S.simuS) = 0.1;                               
ys = spm_int(P,M,U);                           
% Upsample
if S.simuUpsample > 1
    for iX0=1:size(Y.y,2)
        tmp0(:,iX0) = interp(Y.y(:,iX0),S.simuUpsample);
    end
    Y.y = tmp0; clear tmp0

    for iX0=1:size(Y.X0,2)
        tmp0(:,iX0) = interp(Y.X0(:,iX0),S.simuUpsample);
    end
    Y.X0 = tmp0; clear tmp0
end

% PP: 
% %get rid of drifts for now
% Y.X0 = [];
% % Michele : why??

ns = size(Y.y,1);
Y.dt = Y.dt/S.simuUpsample;
% Downsampling of simulated data (16 time bins
% per TR -> experimental TR (or upsampled TR))
% for iModal = 1:size(ys,2)
%     tmppp(:,iModal) = decimate(ys(:,iModal),Y.dt/U.dt);
% end
% ys = tmppp; clear tmppp
ys = ys(( (0:round(Y.dt/U.dt):(size(ys,1)-1)) )+1,:); %
%ys = ys(round((0:(ns - 1))*Y.dt/U.dt)+1,:); %

if S.simuNoise % User-specified response amplitude
    if S.simuA>100 % then reduce the noise level
        Y.y = Y.y/(S.simuA/100) + ys;
    else
        Y.y = Y.y + ys*S.simuA/100;
    end
else
    Y.y = ys*S.simuA/100;
end
% Filter included in function decimate

% Now interpolate simulated or real data to
% (16)
% time bins per TR, to help model estimation
if S.simuInterp>1
    tBase = 0:Y.dt:(Y.dt*(size(Y.y,1)-1));
    tInterp = 0:Y.dt/S.simuInterp:(Y.dt*(size(Y.y,1)-1));
    for iModal = 1:size(ys,2)
        %tmppp(:,iModal) = interp(Y.y(:,iModal),S.simuInterp);
        tmppp(:,iModal) = interp1(tBase,Y.y(:,iModal),tInterp,'spline');
    end
    Y.y = tmppp; clear tmppp
    if ~isempty(Y.X0)
        for iX0=1:size(Y.X0,2)
            tmp0(:,iX0) = interp1(tBase,Y.X0(:,iX0),tInterp,'spline');
            %tmp0(:,iX0) = interp(Y.X0(:,iX0),S.simuInterp);
        end
        Y.X0 = tmp0; clear tmp0
    end
    Y.dt = Y.dt/S.simuInterp;
end

end

function display_results(fullfigDir,Model_Choice,pE,Ep,Cp,U,m,M,H1,K1,K2,...
    modal,plot_algebraic_CMRO2,save_figures,SubjIdx,hfigevolution,hfigevolution2)

    Finter = spm_figure('GetWin','Interactive');
    header = get(Finter,'Name');
    set(Finter,'Name','Hemodynamic modelling')

    t       = [1:M.N]*M.dt;
    Fhdm    = spm_figure;
    set(Fhdm,'name','Hemodynamic Modeling')

    % display input parameters
    %--------------------------------------------------------------------------
    subplot(2,2,1)
    switch Model_Choice
        case 0 %Buxton-Friston
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
        case 1 %Zheng-Mayhew
            P     = Ep(9:end);                  %MODIFIER
            C     = diag(Cp(9:end,9:end));      %MODIFIER
        case 2 %Huppert1
            P     = Ep(11:end);                 %MODIFIER
            C     = diag(Cp(11:end,11:end));    %MODIFIER
        case 3 % BF part 1
            P     = Ep(3:end);
            C     = diag(Cp(3:end,3:end));
        case 4 % BF part 2
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
        otherwise
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
    end

    [dummy, j] = max(abs(P));
    spm_barh(P,C)
    axis square
    title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
    switch Model_Choice
        case {0,1,3,4} %Buxton-Friston
            set(gca,'Ytick',[1:m],'YTickLabel',U.name,'FontSize',8)
            str = {};
            for i = 1:m
                str{end + 1} = U.name{i};
                str{end + 1} = sprintf('mean = %0.2f',P(i));
                str{end + 1} = '';
            end
            set(gca,'Ytick',[1:m*3]/3 + 1/2,'YTickLabel',str)
        case 2
            m2 = 2*m;
            for m3=1:m
                temp_labels{m3} = ['CMRO2' U.name{m3}];
            end
            labels2 = [U.name temp_labels];
            set(gca,'Ytick',[1:m2],'YTickLabel',labels2,'FontSize',8)
            str = {};
            for i = 1:m2
                str{end + 1} = labels2{i}; %U.name{i};
                str{end + 1} = sprintf('mean = %0.2f',P(i));
                str{end + 1} = '';
            end
            set(gca,'Ytick',[1:m2*3]/3 + 1/2,'YTickLabel',str)
        otherwise
    end
    xlabel('relative efficacy per event/sec')


    % display hemodynamic parameters
    %---------------------------------------------------------------------------
    subplot(2,2,3)
    switch Model_Choice
        case {0,4} %Buxton-Friston
            P     = Ep(1:6);
            pE    = pE(1:6);
            C     = diag(Cp(1:6,1:6));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)
        case 1 %Zheng-Mayhew
            P     = Ep(1:8);            %MODIFIER
            pE    = pE(1:8);            %MODIFIER
            C     = diag(Cp(1:8,1:8));  %MODIFIER
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18+6]/3 + 1/2)   %MODIFIER? was 18
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'VASCULAR TONE',...                                %MODIFIER
                sprintf('%0.2f *10 seconds',P(6)),'',...              %MODIFIER
                'GAIN PARAMETER',...                               %MODIFIER
                sprintf('%0.3f *10 seconds',P(7)),'',...              %MODIFIER
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(8),'%'),''},'FontSize',8)   %MODIFIER
        case 2 %Huppert1
            P     = Ep(1:10);
            pE    = pE(1:10);
            C     = diag(Cp(1:10,1:10));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18+12]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'',...
                'CMRO2 SIGNAL decay',...
                sprintf('%0.2f per sec',P(7)),'',...
                'CMRO2 FEEDBACK',...
                sprintf('%0.2f per sec',P(8)),'',...
                'Volume fraction',...
                sprintf('%0.2f',P(9)),'',...
                'O2 Diffusion K',...
                sprintf('%0.2f per sec',P(10)),'','FontSize',8)
         case 3 %Buxton-Friston part 1
            P     = Ep(1:2);
            pE    = pE(1:2);
            C     = diag(Cp(1:2,1:2));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:6]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                ''},'FontSize',8)
        otherwise
            P     = Ep(1:6);
            pE    = pE(1:6);
            C     = diag(Cp(1:6,1:6));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                sprintf('%0.2f per sec',P(1)),'',...
                'FEEDBACK',...
                sprintf('%0.2f per sec',P(2)),'',...
                'TRANSIT TIME',...
                sprintf('%0.2f seconds',P(3)),'',...
                'EXPONENT',...
                sprintf('%0.2f',P(4)),'',...
                'EXTRACTION',...
                sprintf('%0.0f %s',P(5)*100,'%'),'',...
                'log SIGNAL RATIO',...
                sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)
    end


% get display state kernels (i.e. state dynamics)
%==========================================================================

    subplot(3,2,2)
    if plot_algebraic_CMRO2
        tmp_H1 = exp(H1(:,:,j));
        %Algebraic relation for m = CMRO2, in arbitrary units
        %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
        tmp_ma = tmp_H1(:,2) .* tmp_H1(:,4) ./ tmp_H1(:,3);
        tmp_H1 = [tmp_H1 tmp_ma];
        plot(t,tmp_H1)
        axis square
        title({['1st order kernels for ' U.name{j}];...
            'state variables'},'FontSize',9)
        ylabel('normalized values')
        switch Model_Choice
            case {0,4} %Buxton-Friston
                legend('s','f','v','q','ma',0);
            case 1 %Zheng-Mayhew
                legend('s','f','v','q','w','ma',0); %MODIFIER
            case 2 %Huppert1
                legend('s','f','v','q','s2','m','Ct','Cv','ma',0)
            case 3 %Buxton-Friston part 1
                legend('s','f',0);
            otherwise
                legend('s','f','v','q','ma',0);
        end
        grid on
    else
        plot(t,exp(H1(:,:,j)))
        axis square
        title({['1st order kernels for ' U.name{j}];...
            'state variables'},'FontSize',9)
        ylabel('normalized values')
        switch Model_Choice
            case {0,4} %Buxton-Friston
                legend('s','f','v','q',0);
            case 1 %Zheng-Mayhew
                legend('s','f','v','q','w', 0); %MODIFIER
            case 2 %Huppert1
                legend('s','f','v','q','s2','m','Ct','Cv',0)
            case 3 %Buxton-Friston part 1
                legend('s','f',0);
            otherwise
                legend('s','f','v','q',0);
        end
        grid on
    end

    % display output kernels (i.e. BOLD response)
    %--------------------------------------------------------------------------
    subplot(3,2,4)
    plot(t,K1(:,:,j))
    axis square
    switch modal
        case 1
            title({ '1st order kernel';'output: BOLD'},'FontSize',9)
        case 2
            title({ '1st order kernel';'output: BOLD, flow'},'FontSize',9)
        case 3
            title({ '1st order kernel';'output: flow'},'FontSize',9)
        case 4
            switch Model_Choice
                case {0,1,3,4}
                    title({ '1st order kernel';'output: HbT, HbR'},'FontSize',9)
                case 2
                    title({ '1st order kernel';'output: HbT, HbR, HbO'},'FontSize',9)
            end
    end
    %ylabel('normalized flow signal')
    switch modal
        case {1,3}
            ylabel('normalized measure response')
        case {2,4}
            ylabel('normalized measure responses')
    end
    grid on

    subplot(3,2,6)
    axis square
    switch modal
        case {1,2}
            imagesc(t,t,K2(:,:,1,j,j))
            title({ '2nd order kernel';'output: BOLD'},'FontSize',9)
        case 3
            imagesc(t,t,K2(:,:,1,j,j))
            title({ '2nd order kernel';'output: flow'},'FontSize',9)
        case 4
            imagesc(t,t,K2(:,:,2,j,j))
            title({ '2nd order kernel';'output: HbR'},'FontSize',9)
    end
    xlabel({'time {seconds} for'; U.name{j}})
    grid on

    %-Reset title
    %--------------------------------------------------------------------------
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input('Thank you',1,'d')
    if save_figures
        if ~exist(fullfigDir,'dir'), mkdir(fullfigDir); end
        %Save figure
        filen1 = fullfile(fullfigDir,['HDM_' gen_num_str(SubjIdx,3)]);
        fullfigDir2 = fullfile(fullfigDir,'fig');
        if ~exist(fullfigDir2,'dir'),mkdir(fullfigDir2); end
        filen3 = fullfile(fullfigDir2,['HDM_' gen_num_str(SubjIdx,3)]);
        saveas(Fhdm,filen3,'fig');
        print(Fhdm, '-dtiffn', filen1);
        Fsi = spm_figure('GetWin','SI');
        filen2 = fullfile(fullfigDir,['EHDM_' gen_num_str(SubjIdx,3)]);
        filen4 = fullfile(fullfigDir2,['EHDM_' gen_num_str(SubjIdx,3)]);
        saveas(Fsi,filen4,'fig');
        print(Fsi, '-dtiffn', filen2);
        
        try
            close(Fhdm);
        end
    end

end



% function VOI = get_VOI5(SPMfile,ROI,s)
% clear matlabbatch
% matlabbatch{1}.spm.util.voi.spmmat = {SPMfile};
% matlabbatch{1}.spm.util.voi.adjust = 0;
% matlabbatch{1}.spm.util.voi.session = s;
% matlabbatch{1}.spm.util.voi.name = ROI.nameROI;
% matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = ROI.coordinateROI;
% matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = ROI.radiusROI;
% matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
% matlabbatch{1}.spm.util.voi.expression = 'i1';
% spm_jobman('run',matlabbatch);
% %load VOI
% [dir1 dummy] = fileparts(SPMfile);
% load(fullfile(dir1,['VOI_' ROI.nameROI '_' int2str(s) '.mat']));
% VOI.Y = Y;
% VOI.xY = xY;
% end