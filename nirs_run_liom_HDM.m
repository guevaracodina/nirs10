function out = nirs_run_liom_HDM(job) %(xSPM,SPM,hReg)
% user interface for hemodynamic model estimation
% FORMAT [Ep,Cp,K1,K2] = spm_hdm_ui(xSPM,SPM,hReg);
%
% xSPM   - structure containing specific SPM details
% SPM    - structure containing generic  SPM details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Ep     - conditional expectations of the hemodynamic model parameters
% Cp     - conditional  covariance  of the hemodynamic model parameters
% K1     - 1st order kernels
% K2     - 2nd order kernels
%          (see main body of routine for details of model specification)
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_ui.m 3666 2010-01-10 17:54:31Z klaas $
try 
    save_figures = job.save_figures;
catch
    save_figures = 1;
end
try 
    nameHDM = job.nameHDM;
catch
    nameHDM = '';
end

TE    = 0.03; TE_ok = 1; %echo time
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
    end
    HDMfile = fullfile(figDir,'HDM.mat');
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
                    load(fBOLD);
                    %if r1 == 1 && SubjIdx == 1 %some generic information to get
                    [Sess s] = get_session(SPM,cs1);
                    U = get_causes(Sess,Stimuli);
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
                            fASL = fullfile(subjectsASL{SubjIdx},'SPM.mat');
                            VOI_ASL = get_VOI(fASL,ROIs(r1),cs1);
                    end
                    
                    %xY = VOI.xY;
                    %y = VOI.Y; %/100;
                    switch modal
                        case {2,3,5}
                            Y.y(:,2)    = VOI_ASL.Y; %y/100;
                            %Specific to Michèle`s project:
                            %calibration factor for ASL
                            load(fASL);
                            [dir1 fil1] = fileparts(SPM.xY.P(1,:));
                            load(fullfile(dir1,'CBFcalibrFactor.mat'));
                            M.CBFcalibrFactor = calibrFactor;
                    end
                    %Add high pass filter
                    %Y.y(:,2) = ButterHPF(1/SPM.xY.RT,1/128,2,Y.y(:,2));
                    %only rescale BOLD signal, not ASL:
                    y1 = VOI.Y;
                    %y1(1:4) = mean(y1(5:10)); y1(end-4:end) = mean(y1(end-10:end-5));
                    %y1 = ButterHPF(1/SPM.xY.RT,1/128,2,y1);
                    %y1 = 100*(y1 - repmat(median(y1),[size(y1) 1]))./repmat(median(y1),[size(y1,1) 1]); % repmat(std(y1),[size(y1,1) 1]);
                    y1 = (y1 - repmat(mean(y1),[size(y1) 1]))./repmat(mean(y1),[size(y1,1) 1]); % repmat(std(y1),[size(y1,1) 1]);
                    Y.y(:,1) = y1;
                    %rescale to unit variance and zero mean
                    %Y.y = 100*(Y.y - repmat(median(Y.y),[size(Y.y,1) 1]))./ repmat(std(Y.y),[size(Y.y,1) 1]);
                    
                    %-place response and confounds in response structure
                    %--------------------------------------------------------------------------
                    switch modal
                        case {1,2,3}
                            %y      = VOI.xY.u;
                            VOI.xY.RT = SPM.xY.RT;
                            Y.dt   = VOI.xY.RT; %SPM.xY.RT;
                            Y.X0   = VOI.xY.X0;
                    end
                    
                    % estimate
                    %===========================================================================
                    if HDMdisplay || save_figures
                        spm('Pointer','Watch')
                        spm('FigName','Estimation in progress');
                    end
                    
                    % Model specification: m input; 4 states; 1 outout; m + 6 parameters
                    %---------------------------------------------------------------------------
                    % u(m) - mth stimulus function     (u)
                    %
                    % x(1) - vascular signal           log(s)
                    % x(2) - rCBF                      log(f)
                    % x(3) - venous volume             log(v)
                    % x(4) - deoyxHb                   log(q)
                    % x(5) - normalized vascular tone  log(w)
                    %
                    % y(1) - BOLD                      (y)
                    %
                    % P(1)       - signal decay               d(ds/dt)/ds)      half-life = log(2)/P(1) ~ 1sec
                    % P(2)       - autoregulation             d(ds/dt)/df)      2*pi*sqrt(1/P(1)) ~ 10 sec
                    % P(3)       - transit time               (t0)              ~ 1 sec
                    % P(4)       - exponent for Fout(v)       (alpha)           c.f. Grubb's exponent (~ 0.38)
                    % P(5)       - resting oxygen extraction  (E0)              ~ range 20 - 50%
                    % P(6) - time constant of vascular tone w               (tau_w)   ~range: 0 - 30?
                    % P(7) - gain parameter b = b0 V0                       (b)       ~range: 0 - 30?
                    % P(8) - ratio of intra- to extra-vascular components   (epsilon)  ~range 0.5 - 2
                    %          of the gradient echo signal
                    %
                    % P(8 + 1:m) - input efficacies - d(ds/dt)/du)  ~0.3 per event
                    %--------------------------------------------------------------------------
                    
                    % priors (3 modes of hemodynamic variation)
                    %--------------------------------------------------------------------------
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
                            M.l = 1;
                        case 2
                            M.l = 2;
                        case 4
                            switch job.Model_Choice
                                case {0,1} %Buxton, Zheng-Mayhew
                                    M.l = 2;
                                case 2 %Huppert1
                                    M.l = 3;
                            end
                    end
                    
                    M.N     = 64;
                    M.dt    = 24/M.N; %24/M.N;
                    M.TE    = TE;
                    
                    %--------------------------------------------------------------------------
                    warning('off','MATLAB:nearlySingularMatrix');
                    warning('off','MATLAB:singularMatrix');
                    % nonlinear system identification
                    %--------------------------------------------------------------------------
                    [Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F] = spm_nlsi(M,U,Y);
                    %Store results
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
                            %Save figure
                            filen1 = fullfile(figDir,['HDM_' gen_num_str(SubjIdx,3)]);
                            figDir2 = fullfile(figDir,'fig');
                            if ~exist(figDir2,'dir'),mkdir(figDir2); end
                            filen3 = fullfile(figDir2,['HDM_' gen_num_str(SubjIdx,3)]);
                            saveas(Fhdm,filen3,'fig');
                            print(Fhdm, '-dtiffn', filen1);
                            Fsi = spm_figure('GetWin','SI');
                            filen2 = fullfile(figDir,['EHDM_' gen_num_str(SubjIdx,3)]);   
                            filen4 = fullfile(figDir2,['EHDM_' gen_num_str(SubjIdx,3)]);   
                            saveas(Fsi,filen4,'fig');
                            print(Fsi, '-dtiffn', filen2);
                            try
                                close(Fhdm);
                            end
                        end
                    end
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