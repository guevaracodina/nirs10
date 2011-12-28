function out = nirs_run_liom_HDM(job) %(xSPM,SPM,hReg)
% user interface for hemodynamic model estimation
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_ui.m 3666 2010-01-10 17:54:31Z klaas $
save_figures = job.save_figures;
nameHDM = job.nameHDM;
TE = job.echo_time;
dp_start = job.dp_start; %points to remove due to filter set-up
dp_end = job.dp_end; %points to remove due to filter set-up
removeWhitening = job.removeWhitening;
if isfield(job.simuOn,'simuYes')
    simuOn = 1;
    simuS     = job.simuOn.simuYes.simuS; %Stimuli types to include
    simuIt    = job.simuOn.simuYes.simuIt; %Number of random iterations
    simuA     = job.simuOn.simuYes.simuA; %Signal amplitude, as % of BOLD signal
    simuP     = job.simuOn.simuYes.simuP; %Parameters to vary
    simuPrior = job.simuOn.simuYes.simuPrior; %Priors to use
    simuR     = job.simuOn.simuYes.simuR; %Range to sample
    simuUpsample = job.simuOn.simuYes.simuUpsample; %Upsampling factor on data
    simuNoise = job.simuOn.simuYes.simuNoise; %Yes to include background noise based on restscans
    restscans = job.simuOn.simuYes.restscans; %Rest scans to add signal to
    restscans_BOLD = job.simuOn.simuYes.restscans_BOLD; %Rest scans to add signal to
    restscans_ASL = job.simuOn.simuYes.restscans_ASL; %Rest scans to add signal to
else
    simuOn = 0;
end
%Individual figure display
HDMdisplay = job.HDMdisplay;
%---------------------------------------------------------------------------
try
    if HDMdisplay || save_figures
        Finter = spm_figure('GetWin','Interactive');
        header = get(Finter,'Name');
        set(Finter,'Name','Hemodynamic modelling')
    end
    %Algebraic relation for m = CMRO2, in arbitrary units
    %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
    plot_algebraic_CMRO2 = 1; %Boolean
    
    % inputs
    %==========================================================================
    
    %Modalities:
    %1: BOLD
    %2: BOLD+ASL
    %3: ASL (not coded up)
    %4: HbO HbR (not implemented in this version)
    %5: BOLD+ASL v2 (not coded up)
    
    %Step 1: checking which modality to run
    if isfield(job.xSPM_Modalities,'xSPM_BOLD')
        modal = 1;
        sBOLD = job.xSPM_Modalities.xSPM_BOLD;
    else
        if isfield(job.xSPM_Modalities,'xSPM_BOLD_ASL')
            modal = 2;
            sBOLD = job.xSPM_Modalities.xSPM_BOLD_ASL;
            %currently not used
            sASL = job.xSPM_Modalities.xSPM_BOLD_ASL;
            subjectsASL = sASL.which_subjects_ASL;
        else
            if isfield(job.xSPM_Modalities,'xSPM_ASL')
                modal = 3;
                sASL = job.xSPM_Modalities.xSPM_ASL;
            else
                if isfield(job.xSPM_Modalities,'xSPM_BOLD_ASL_V2')
                    modal = 5;
                    sBOLD = job.xSPM_Modalities.xSPM_BOLD_ASL_V2;
                end
            end
        end
    end
    %load SPM and xSPM
    clear SPM xSPM
    switch modal
        case {1,2,5}
            spmmat = sBOLD.spmmat{1};
            xSPMmat = sBOLD.xSPM_spmmat{1};
            load(spmmat);
            SPM_BOLD = SPM;
            load(xSPMmat);
            xSPM_BOLD = xSPM;
    end
    clear SPM xSPM
    switch modal
        case {3,2}
            spmmat = sASL.spmmat_ASL{1};
            load(spmmat);
            SPM_ASL = SPM;
            xSPMmat = sASL.xSPM_spmmat_ASL{1};
            load(xSPMmat);
            xSPM_ASL = xSPM;
    end
    clear SPM xSPM
    %get rootDir
    tmp = strfind(spmmat,filesep);
    rootDir = spmmat(1:tmp(end-2));
    if save_figures
        if ~isempty(nameHDM)
            figDir = fullfile(rootDir,['HDM_' nameHDM]);
        else
            figDir = fullfile(rootDir,'HDM');
        end
        if ~exist(figDir,'dir'),mkdir(figDir);end
        HDMfile = fullfile(figDir,'HDM.mat');
    else
        HDMfile = fullfile(rootDir,'HDM.mat');
    end
    
    %Various inputs:
    %Choice of model:
    Model = job.Model_Choice;
    %Stimuli:
    Stimuli = job.Stimuli;
    %StimuliSign = job.StimuliSign;
    %Sessions:
    sessions = job.which_session;
    %subjects:
    subjects = job.which_subjects;
    %ROIs
    ROIs = job.whichROI;
    
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
            for SubjIdx=1:length(subjects)
                try
                    %load data
                    clear SPM Y
                    fBOLD= fullfile(subjects{SubjIdx},'SPM.mat');
                    fBOLD_old = fBOLD;
                    load(fBOLD);
                    if removeWhitening
                        switch modal
                            case 1
                        SPM.xX.W = speye(size(SPM.xX.W,1));
                            case 5
                                %quick fix to use only one session
                                 %SPM.xX.W = speye(size(SPM.xX.W,1)/2);
                                 %need to subtract 1 image too
                                 SPM.xX.W = speye(size(SPM.xX.W,1)/2-1);
                                 SPM.xX.K = SPM.xX.K(1);
                                 SPM.xX.K.row = SPM.xX.K.row(1:end-1);
                                 SPM.xX.K.X0 = SPM.xX.K.X0(1:end-1,:);
                                 %Remove HPF entirely for now
                                 SPM.xX.K.X0 =zeros(size(SPM.xX.K.X0));
                        end
                        %save it as a temporary file
                        [dir0 fil0 ext0] = fileparts(fBOLD);
                        tmpdir = fullfile(dir0,'tmp_noWhitening');
                        if ~exist(tmpdir,'dir')
                            mkdir(tmpdir)
                            fBOLD = fullfile(tmpdir,[fil0 ext0]);
                            save(fBOLD,'SPM');
                        end
                    end
                    switch modal
                        case 2
                            fASL = fullfile(subjectsASL{SubjIdx},'SPM.mat');
                    end
                    
                    %if r1 == 1 && SubjIdx == 1 %some generic information to get
                    [Sess s] = get_session(SPM,cs1);
                    U = get_causes(Sess,Stimuli);
                    
                    %if simuOn, replace the data with the rest data
                    if simuOn
                        %                        if simuNoise
                        switch modal
                            case 1
                                SPM.xY.P     = char(restscans{:});
                                SPM.xY.VY = spm_vol(SPM.xY.P);
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fBOLD);
                                tmpdir = fullfile(dir00,'tmp_simu');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                fBOLD_old = fBOLD;
                                fBOLD = fullfile(tmpdir,[fil00 ext00]);
                                save(fBOLD,'SPM');
                            case 5
                                SPM.xY.P     = char(restscans_BOLD{:});
                                SPM.xY.VY = spm_vol(SPM.xY.P);
                                %save SPM in a temporary location
                                [dir00 fil00 ext00] = fileparts(fBOLD);
                                tmpdir = fullfile(dir00,'tmp_simu_BOLD');
                                if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
                                fBOLD_old = fBOLD;
                                fBOLD = fullfile(tmpdir,[fil00 ext00]);
                                save(fBOLD,'SPM');
                                
                                SPM.xY.P     = char(restscans_ASL{:});
                                SPM.xY.VY = spm_vol(SPM.xY.P);
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
                    %cwd = pwd;
                    %extract VOI for BOLD
                    switch modal
                        case {1,2}
                            VOI = get_VOI(fBOLD,ROIs(r1),cs1);
                        case 5
                            VOI = get_VOI5(fBOLD,ROIs(r1),cs1);
                    end
                    %extract VOI for ASL
                    switch modal
                        case 2                            
                            VOI_ASL = get_VOI(fASL,ROIs(r1),cs1);
                        case 5
                            VOI_ASL = get_VOI5(fASL,ROIs(r1),cs1);
                    end
                    %cd(cwd);
                    %xY = VOI.xY;
                    %y = VOI.Y; %/100;
                    switch modal
                        case {2,3,5}
                            Y.y(:,2)    = VOI_ASL.Y((1+dp_start):end-dp_end); %y/100;
                            %Specific to Michèle`s project:
                            %calibration factor for ASL
                            load(fASL);
                            
                            [dir1 fil1] = fileparts(SPM.xY.P(1,:));
                            switch modal
                                case {2,3}
                                    
                                    load(fullfile(dir1,'CBFcalibrFactor.mat'));
                                case 5
                                    load(fullfile(dir1,'CBFcalibrFactor.mat'));
                                    %load(fullfile(dir1,'flow','CBFcalibrFactor.mat'));
                            end
                            M.CBFcalibrFactor = calibrFactor;
                    end
                    
                    %Add high pass filter
                    %Y.y(:,2) = ButterHPF(1/SPM.xY.RT,1/128,2,Y.y(:,2));
                    %only rescale BOLD signal, not ASL:
                    y1 = VOI.Y;
                    
                    y1 = y1((1+dp_start):end-dp_end);
                    %y1 = ButterHPF(1/SPM.xY.RT,1/128,2,y1);
                    %y1(1:4) = mean(y1(5:10)); y1(end-4:end) = mean(y1(end-10:end-5));
                    
                    %y1 = 100*(y1 - repmat(median(y1),[size(y1) 1]))./repmat(median(y1),[size(y1,1) 1]); % repmat(std(y1),[size(y1,1) 1]);
                    tmp_m = repmat(mean(y1),[size(y1,1) 1]);
                    y1 = (y1 - tmp_m)./std(y1); %rescale to zero mean and unit standard deviation %tmp_m;
                    Y.y(:,1) = y1;
                    %rescale to unit variance and zero mean
                    %Y.y = 100*(Y.y - repmat(median(Y.y),[size(Y.y,1) 1]))./ repmat(std(Y.y),[size(Y.y,1) 1]);
                    
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
                    % estimate
                    %===========================================================================
                    if HDMdisplay || save_figures
                        spm('Pointer','Watch')
                        spm('FigName','Estimation in progress');
                    end
                    %number of inputs
                    m       = size(U.u,2);
                    switch Model
                        case 0 %Buxton-Friston
                            [pE,pC] = spm_hdm_priors_YO(m); %,3); %PP pour Michèle.
                        case 1 %Zheng-Mayhew
                            [pE,pC] = nirs_hdm_priors_ZM(m,5);   %MODIFIER
                        case 2 %Huppert1
                            [pE,pC] = nirs_hdm_priors_Huppert1(m,3);
                        otherwise
                            [pE,pC] = spm_hdm_priors(m,3);
                    end
                    if simuOn
                        pA = repmat(pE',[simuIt 1]); %zeros(simuIt,length(pE));
                        ct = 0;
                        for pE1=1:length(pE)
                            if any(simuP == 0) || any(pE1 == simuP)
                                ct = ct+1;
                                %initialize the stream
                                mtstream = RandStream('mt19937ar','Seed',pE1);
                                RandStream.setDefaultStream(mtstream);
                                for it1=1:simuIt
                                    %generate the random numbers
                                    if ~isempty(simuPrior)
                                        pA(it1,pE1) = simuPrior(ct);
                                    end
                                    tpA = pA(it1,pE1);
                                    if length(simuR) == 1
                                        pA(it1,pE1) = unifrnd(tpA*(1-simuR/100),tpA*(1+simuR/100));
                                    else
                                        pA(it1,pE1) = unifrnd(tpA*(1-simuR(ct)/100),tpA*(1+simuR(ct)/100));
                                    end
                                end
                            end
                        end
                    end
                    % model
                    %--------------------------------------------------------------------------
                    M.modal = modal;
                    M.Model_Choice = job.Model_Choice;
                    switch job.Model_Choice
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
                        otherwise %Buxton-Friston
                            M.f     = 'spm_fx_hdm';
                            M.g     = 'nirs_gx_hdm';
                            M.x     = [0 0 0 0]';
                            M.n     = 4;
                    end
                    if HDMdisplay || save_figures
                        M.nograph = 0;
                    else
                        M.nograph = 1;
                    end
                    % NB: resting value/expansion point of x(1) is -Inf in log space; this is
                    % taken into account in spm_fx_hdm.
                    M.pE    = pE;
                    M.pC    = pC;
                    M.m     = m;
                    
                    switch modal
                        case {1,3}
                            M.l = 1; %BOLD
                        case 2
                            M.l = 2;
                        case 4
                            switch job.Model_Choice
                                case {0,1} %Buxton, Zheng-Mayhew
                                    M.l = 2;
                                case 2 %Huppert1
                                    M.l = 3;
                            end
                        case 5
                            M.l = 2; %BOLD, ASL
                    end
                    
                    M.N     = 64;
                    M.dt    = 24/M.N; %24/M.N;
                    M.TE    = TE;
                    
                    if ~simuOn
                        simuIt = 1;
                    else
                        Y0 = Y;
                        if simuS == 0
                            simuS = 1:size(U.u,2);
                        end
                    end
                    %--------------------------------------------------------------------------
                    warning('off','MATLAB:nearlySingularMatrix');
                    warning('off','MATLAB:singularMatrix');
                    for it1=1:simuIt
                        if simuOn
                            %Buxton-Friston balloon and other models
                            P = pA(it1,:)';
                            %efficacies
                            P(end-size(U.u,2)+simuS) = 1;
                            %tic
                            %Careful! must use the same integrator for both the direct and the inverse model
                            %otherwise, there are systematic biases in parameter estimations
                            %ys = spm_int_J(P,M,U); %6.7 times slower than spm_int_D, but produces a very different result
                            ys = spm_int(P,M,U); 
                            %ys = spm_int_D(P,M,U);
                            %toc
                            %set to zero mean and rescale to unit standard deviation
                            % %                             tmp_m = repmat(mean(ys),[size(ys,1) 1]);
                            % %                             ys = (ys - tmp_m)./std(ys);
                            %                             %Canonical response
                            %                             ns = size(Y.y,1);
                            %                             xBF.T = 10;
                            %                             xBF.T0 = 0;
                            %                             xBF.dt = Y.dt/xBF.T; % - time bin length {seconds}
                            %                             xBF.name = 'hrf'; %description of basis functions specified
                            %                             xBF = spm_get_bf(xBF);
                            %                             %unit area
                            %                             xBF.bf = xBF.bf/sum(xBF.bf); %normalize
                            %                             %convolve stimuli U with basis functions
                            %                             [X,Xn,Fc] = spm_Volterra(U,xBF.bf,1); %
                            %                             X = X((0:(ns - 1))*xBF.T + xBF.T0 + 32,:);
                            %                             %add response to rest data -- direct model
                            Y = Y0;
                            %Y.y = Y0.y + sum(X(:,simuS),2)*simuA/100;
                            %Upsample
                            if simuUpsample > 1
                                switch modal 
                                    case 1
                                        Y.y = interp(Y.y,simuUpsample);
                                       
                                    case {2,5}
                                        for iX0=1:size(Y.y,2)
                                            tmp0(:,iX0) = interp(Y.y(:,iX0),simuUpsample);
                                        end
                                        Y.y = tmp0; clear tmp0
                                end
                                
                                for iX0=1:size(Y.X0,2)
                                    tmp0(:,iX0) = interp(Y.X0(:,iX0),simuUpsample);
                                end
                                Y.X0 = tmp0; clear tmp0
                                
                            end
                            %get rid of drifts for now
                            Y.X0 = [];
                            ns = size(Y.y,1);
                            Y.dt = Y.dt/simuUpsample;
                            ys = ys(round((0:(ns - 1))*Y.dt/U.dt)+1,:); %????
                            
                            if simuNoise
                                Y.y = Y.y + ys*simuA/100;
                            else
                                Y.y = ys*simuA/100;
                            end
                            %Low pass filtering of the data after downsampling -- otherwise there will be aliasing
                            if simuUpsample < 16 && simuUpsample > 1 %here Y0 might get LPF twice...
                                Y.y = ButterLPF(1/Y.dt,0.95*1/(2*Y.dt),3,Y.y);
                            end
                        end
                        % nonlinear system identification
                        %--------------------------------------------------------------------------
                        [Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F] = nirs_nlsi(M,U,Y);
                        %for simulations: store results in place of subjects
                        if simuOn
                            SubjIdx0 = SubjIdx; %should not be used
                            SubjIdx = it1; %changed inside a for loop but it is restored later
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
                        HDM{SubjIdx,r1}{s1}.subj = fBOLD; %to remove any doubt
                        HDM{SubjIdx,r1}{s1}.session = cs1;
                        %-display results
                        %==========================================================================
                        t       = [1:M.N]*M.dt;
                        if HDMdisplay || save_figures
                            Fhdm    = spm_figure;
                            set(Fhdm,'name','Hemodynamic Modeling')
                            
                            
                            % display input parameters
                            %--------------------------------------------------------------------------
                            subplot(2,2,1)
                            switch job.Model_Choice
                                case 0 %Buxton-Friston
                                    P     = Ep(7:end);
                                    C     = diag(Cp(7:end,7:end));
                                case 1 %Zheng-Mayhew
                                    P     = Ep(9:end);                  %MODIFIER
                                    C     = diag(Cp(9:end,9:end));      %MODIFIER
                                case 2 %Huppert1
                                    P     = Ep(11:end);                 %MODIFIER
                                    C     = diag(Cp(11:end,11:end));    %MODIFIER
                                otherwise
                                    P     = Ep(7:end);
                                    C     = diag(Cp(7:end,7:end));
                            end
                            
                            [dummy, j] = max(abs(P));
                            spm_barh(P,C)
                            axis square
                            title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
                            switch job.Model_Choice
                                case {0,1} %Buxton-Friston
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
                            switch job.Model_Choice
                                case 0 %Buxton-Friston
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
                            
                        end
                        % get display state kernels (i.e. state dynamics)
                        %==========================================================================
                        
                        % Volterra kernels of states
                        %--------------------------------------------------------------------------
                        [dummy,H1] = spm_kernels(M0,M1,M.N,M.dt);
                        HDM{SubjIdx,r1}{s1}.H1 = H1;
                        %save - overwriting previous
                        save(HDMfile,'HDM');
                        
                        if HDMdisplay || save_figures
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
                                switch job.Model_Choice
                                    case 0 %Buxton-Friston
                                        legend('s','f','v','q','ma',0);
                                    case 1 %Zheng-Mayhew
                                        legend('s','f','v','q','w','ma',0); %MODIFIER
                                    case 2 %Huppert1
                                        legend('s','f','v','q','s2','m','Ct','Cv','ma',0)
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
                                switch job.Model_Choice
                                    case 0 %Buxton-Friston
                                        legend('s','f','v','q',0);
                                    case 1 %Zheng-Mayhew
                                        legend('s','f','v','q','w', 0); %MODIFIER
                                    case 2 %Huppert1
                                        legend('s','f','v','q','s2','m','Ct','Cv',0)
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
                                    title({ '1st order kernel';'output: BOLD, ASL'},'FontSize',9)
                                case 3
                                    title({ '1st order kernel';'output: ASL'},'FontSize',9)
                                case 4
                                    switch job.Model_Choice
                                        case {0,1}
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
                                    title({ '2nd order kernel';'output: ASL'},'FontSize',9)
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
                                fullfigDir = fullfile(figDir,['S' int2str(cs1) '_' ROIs(r1).nameROI]);
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
                        if simuOn
                            SubjIdx = SubjIdx0;
                        end
                    end
                    %                     %remove temp files -- careful, this is not coded correctly, so a good SPM.mat might get deleted by mistake
                    %                     if removeWhitening
                    %                         if simuOn
                    %                             delete(fBOLD_old);
                    %                         end
                    %                     end
                    %                     if simuOn
                    %                         delete(fBOLD);
                    %                     end
                    out = [];
                catch exception
                    disp(exception.identifier);
                    disp(exception.stack(1));
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


function VOI = get_VOI5(SPMfile,ROI,s)
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