function out = nirs_run_NIRS_HDM(job) %(xSPM,SPM,hReg) 
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

%addpath('D:\Users\Philippe Pouliot\spm8\toolbox\dcm_liom\')
%Usage, from command line: Call to function
%   [Ep,Cp,K1,K2] = liom_hdm_ui_ZM(xSPM,SPM,hReg);
% get figure handles


%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Hemodynamic modelling')

%Algebraic relation for m = CMRO2, in arbitrary units
%m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
plot_algebraic_CMRO2 = 1; %Boolean 
        
% inputs
%==========================================================================

try 
    spmmat = job.Modalities.BOLD.spmmat{1,1};
    modal = 1; %BOLD only
catch
    try
        spmmat = job.Modalities.BOLD_ASL.spmmat{1,1};
        spmmat_ASL = job.Modalities.BOLD_ASL.spmmat_ASL{1,1};
        modal = 2; %BOLD and ASL 
    catch
        try
            spmmat_ASL = job.Modalities.ASL.spmmat_ASL{1,1};
            modal = 3; %ASL only
        catch
            try
                ch_mode = job.Modalities.HbO_HbR.ch_mode;
                spmmat_HbO = job.Modalities.HbO_HbR.spmmat_HbO{1,1};
                spmmat_HbR = job.Modalities.HbO_HbR.spmmat_HbR{1,1};
                preproc_file =job.Modalities.HbO_HbR.NIRS_SPM_Coregistration_Channels{1,1};
                view = job.Modalities.HbO_HbR.view;
                rescalingFactor1 = job.Modalities.HbO_HbR.rescalingFactor1;
                downsamplingFactor1 = job.Modalities.HbO_HbR.downsamplingFactor1;
                modal = 4;
            catch
                disp('Unrecognized modality or missing inputs, aborting');
                %return
            end
        end
    end
end

try 
    modal;
catch
    modal = spm_input('1:BOLD only 2:BOLD+ASL 3:ASL only',1,'n1',3,3);
end
            
switch modal
    case {1,2}     
        try
            spmmat;
        catch
            spmmat = spm_select(1,'^SPM\.mat$','Select SPM.mat of BOLD estimation');
        end
        load(fullfile(spmmat));
        SPM_BOLD = SPM;
    case 4 %HbO+HbR
        cinterp_SPM_nirs = [];
        preproc_info = [];
        SPM_nirs = [];
        
        try
            load(preproc_file);
            %preproc_info;
            %give chpos.rchn and chpos.cchn (2D projection coordinates as row
            %and column)
            chpos = preproc_info.rend_ch_pos{view};
        catch
            disp('Could not load coregistration file');
        end
        try
            ch_num_avg = ch_mode.ch_avg_mode.ch_num_avg; 
        catch
            try
                map_file = ch_mode.stat_map_mode.map_file{1,1};
                load(map_file);
                try
                    %in case this map comes from a group estimation    
                    tstat = SPM_nirs.nirs.tstat;
                catch
                    try 
                        cinterp_SPM_nirs;
                        %need to calculate tstat: (without correction for tube formula
                        %-- see nirs_run_NIRS_SPM_contrast for code (which could be
                        %saved in cinterp_ files if desired
                        tstat = cinterp_SPM_nirs.cbeta ./sqrt(cinterp_SPM_nirs.ccovbeta);
                        tstat = reshape(tstat,cinterp_SPM_nirs.s1,cinterp_SPM_nirs.s2);
                    catch
                        disp('Content of map file is unrecognized.');
                    end
                end
            catch
                disp('Could not load statistical map file');
            end
        end
        try
            load(spmmat_HbO);
            SPM_HbO = SPM_nirs;
        catch
            disp('Could not load HbO data file');
        end
        try
            load(spmmat_HbR);
            SPM_HbR = SPM_nirs;
        catch
            disp('Could not load HbR data file');
        end
end

switch modal
    case {2,3}      
        try 
            spmmat_ASL;
        catch
            spmmat = spm_select(1,'^SPM\.mat$','Select SPM.mat of ASL estimation');
        end
        load(fullfile(spmmat_ASL));
        SPM_ASL = SPM;
        xSPM.swd  = SPM_ASL.swd;
        disp('Questions regarding ASL');
        [hReg_ASL,xSPM_ASL,SPM] = spm_results_ui('Setup',xSPM);

        [dummy,Q_ASL] = max(xSPM_ASL.Z); 
        xyz_ASL = xSPM_ASL.XYZ(:,Q_ASL);
        xSPM_ASL.xY.xyz = xyz_ASL;   
        [xyzmm,d] = spm_XYZreg('SetCoords',xSPM_ASL.XYZmm(:,Q_ASL),hReg_ASL);
        xSPM_ASL.xY.xyz = xyzmm;
end

switch modal
    case {1,2} 
        SPM = SPM_BOLD;
    case 3
        SPM = SPM_ASL;
    case 4
        SPM = SPM_HbO;
end

% which session?
%---------------------------------------------------------------------------
s    = length(SPM.Sess);
if s > 1
    s = spm_input('which session',1,'n1',s,s);
end
Sess = SPM.Sess(s);

% 'causes' or imputs U
%---------------------------------------------------------------------------
try 
    u_list = job.Stimuli; %1; %LFCP
    u = length(u_list);
    flagu = 1;
catch
    flagu = 0;
    u    = length(Sess.U);
    spm_input('Input specification:...  ',1,'d');
end

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


%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------

TE    = 0.03;
switch modal
    case {1,2,3}
        TE_ok = 0;
        while ~TE_ok
            TE = spm_input('Echo time, TE [s]', '+1', 'r', TE);
            if ~TE || (TE < 0) || (TE > 0.1)
                str = { 'Extreme value for TE or TE undefined.',...
                    'Please re-enter TE (in seconds!)'};
                spm_input(str,'+1','bd','OK',[1],1);
            else
                TE_ok = 1;
            end
        end
    case 4
        TE_ok = 1; %no fMRI
end
% system outputs
%===========================================================================

% enforce adjustment w.r.t. all effects
%---------------------------------------------------------------------------
xY     = struct(    'Ic'        ,1,...  
            'name'      ,'HDM',...
            'Sess'      ,s);

switch modal
%     case 1 %BOLD only
%         [y xY] = spm_regions(xSPM,SPM,hReg,xY);
%         y = xY.u;
%         Y.y    = y;
    case 2 %BOLD and ASL
        xY2 = xY;
        [y xY2] = spm_regions(xSPM_ASL,SPM_ASL,hReg_ASL,xY2);
        y = xY2.u;
        Y.y(:,2)    = y/100;
    case 3 %ASL only
        [y xY] = spm_regions(xSPM_ASL,SPM_ASL,hReg_ASL,xY);
        y = xY.u;
        Y.y    = y/100; %Scaling factor --should be 100? why is this required??? 
    case 4
        %Specify channel to look at
        %Find the channel with the minimum distance to the map's maximum
        try 
            %mode stat_map
            ch_mode.stat_map_mode.map_file{1,1};
            [r1 c1] = find(tstat == max(max(tstat)));
            d1 = zeros(size(chpos.rchn,1),1);
            for i1=1:size(chpos.rchn,1)
                d1(i1) = sum(([r1; c1]-[chpos.rchn(i1); chpos.cchn(i1)]).^2).^0.5;
            end
            [vd1 id1] = min(d1);
            disp(['Channel is: ' int2str(id1) ' and distance to max is: ' num2str(vd1)]); 
            %id1 is now the index of the channel closest to the maximum of the stat map
            load(SPM_HbO.nirs.fname);
            y(:,1) = SPM_nirs.nirs.KY(:,id1); %nirs_data.oxyData(:,id1); %HbO
            %load(SPM_HbR.nirs.fname); should really be coming from the same
            %file
            load(SPM_HbR.nirs.fname)
            y(:,2) = SPM_nirs.nirs.KY(:,id1); %nirs_data.dxyData(:,id1); %HbR
        catch
            try
                %mode channel average
                ch_mode.ch_avg_mode.ch_num_avg;
                load(SPM_HbO.nirs.fname);                
                %average of selected channels
                y(:,1) = sum(SPM_nirs.nirs.KY(:,ch_num_avg),2)/length(ch_num_avg);
                load(SPM_HbR.nirs.fname);
                y(:,2) = sum(SPM_nirs.nirs.KY(:,ch_num_avg),2)/length(ch_num_avg); 
            catch
                disp('Unrecognized channel mode');
            end
        end
            
end

switch modal
    case {1,2} 
        xSPM = [];
        xSPM.swd  = SPM_BOLD.swd;
        SPM = SPM_BOLD;
        disp('Questions regarding BOLD');
        [hReg_BOLD,xSPM_BOLD,SPM] = spm_results_ui('Setup',xSPM);
    case 4
    otherwise
end

switch modal
    case 1
        [dummy,Q_BOLD] = max(xSPM_BOLD.Z);
    case 2
        Q_BOLD = spm_XYZreg('FindXYZ',xyz_ASL,xSPM_BOLD.XYZ);
        if isempty(Q_BOLD) 
            [xyz_BOLD,Q_BOLD,d_BOLD] = spm_XYZreg('NearestXYZ',xyz_ASL,xSPM_BOLD.XYZ);
        end
    case 4
    otherwise
end

switch modal
    case {1,2}
        [xyzmm,d] = spm_XYZreg('SetCoords',xSPM_BOLD.XYZmm(:,Q_BOLD),hReg_BOLD);
        xSPM_BOLD.xY.xyz  = xyzmm;   
        % get region stucture for BOLD
        %---------------------------------------------------------------------------
        [y xY] = spm_regions(xSPM_BOLD,SPM_BOLD,hReg_BOLD,xY);   
    case 4
    otherwise
end

switch modal
    case 1
        Y.y      = y;
    case 2
        Y.y(:,1) = y;
    case 4
        Y.y      = y;
    otherwise
end

%-place response and confounds in response structure
%--------------------------------------------------------------------------
%-
switch modal
    case {1,2,3}
        y      = xY.u;
        Y.dt   = SPM.xY.RT;
        Y.X0   = xY.X0;
    case 4
        % KY = spm_filter_HPF_LPF_WMDL(SPM_HbR.xX.K.KL, nirs_data.oxyData(:,99));
        %y(:,1) = y(:,1) - spm_filter_HPF_LPF_WMDL(SPM_HbO.xX.K.KL, y(:,1));
        %y(:,2) = y(:,2) - spm_filter_HPF_LPF_WMDL(SPM_HbR.xX.K.KL, y(:,2));
        %y(:,1) = spm_filter_HPF_LPF_WMDL(SPM_HbO.xX.K.KL, y(:,1));
        %y(:,2) = spm_filter_HPF_LPF_WMDL(SPM_HbR.xX.K.KL, y(:,2));
        switch job.Model_Choice
            case {0,1} %Buxton, Zheng-Mayhew
                y(:,1) = y(:,1) + y(:,2); %Assume v = HbT = HbO+HbR
            case 2
                y(:,3) = y(:,1); %HbO
                y(:,1) = y(:,1) + y(:,2); %Assume v = HbT = HbO+HbR
            otherwise
        end
                
        Y.dt = SPM_HbO.xY.RT;
        %Y.X0 = []; %ones(size(y,1),1); %Constant?
        
        %Lower the sampling rate
        %Resample = rescalingFactor15; %long (15 to 30 mins) estimation if Resample < 20, e.g. tested at 5.
        Y.dt = Y.dt*downsamplingFactor1;
        y = y(1:downsamplingFactor1:end,:);
        Y.y = y*rescalingFactor1;  %Scaling factor - required for convergence
        %U.dt = U.dt * Resample;
        %U.u = U.u(1:Resample:end);
        %enlarge a bit U.u so that spm_int doesn't break
        %U.u = sparse([1:size(U.u,1)],1,U.u,175800,1);
    otherwise         
end

% estimate
%===========================================================================
spm('Pointer','Watch')
spm('FigName','Estimation in progress');


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
switch job.Model_Choice
    case 0 %Buxton-Friston
        [pE,pC] = spm_hdm_priors(m,3);
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
M.dt    = 24/M.N;
M.TE    = TE;

%--------------------------------------------------------------------------

% nonlinear system identification
%--------------------------------------------------------------------------
[Ep,Cp,dummy,dummy2,K1,K2,M0,M1] = nirs_nlsi(M,U,Y);

%-display results
%==========================================================================
t       = [1:M.N]*M.dt;
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


% get display state kernels (i.e. state dynamics) 
%==========================================================================

% Volterra kernels of states
%--------------------------------------------------------------------------
[dummy,H1] = spm_kernels(M0,M1,M.N,M.dt);

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

out = [];
end