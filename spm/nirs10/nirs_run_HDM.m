function out = nirs_run_HDM(job)
% User interface for hemodynamic model estimation
% Based on SPM HDM
Model_Choice = job.Model_Choice;
Stimuli = job.which_condition;
S = nirs_get_simu(job);
EM = nirs_get_EM(job);
O = nirs_get_common_O_HDM_SCKS(job);
DO = nirs_get_common_hdm_display_options(job);

% Loop over subjects
for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        clear NIRS HDM0
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'HDM_OK') || job.force_redo)
            
            Sess = NIRS.Dt.fir.Sess;
            %loop over each session
            for s1=1:length(Sess)
                if all_sessions || sum(s1==selected_sessions)
                    %get onsets in right format for HDM
                    U = get_causes(Sess,Stimuli);
                    %get data - selected channels        
                    
                    %Loop over channels -- could this be done in parallel?        
                            Y.dt
                            Y.X0
                            
                            
                                    
                            if generate_figures || save_figures
                                spm('Pointer','Watch')
                                spm('FigName','Estimation in progress');
                            end
                            
                            % Priors
                            %--------------------------------------------------------------------------
                            %number of inputs
                            m       = size(U.u,2);
                            [pE,pC] = set_priors(priorFile,m);
                            
                            % Generate random parameter distributions
                            if S.simuOn
                                pA = nirs_generate_random_param(S,pE);
                            end
                            
                            
                            HDM_str = ['S' gen_num_str(s1,2) '_ROI' gen_num_str(r1,3)];
                                        HDMfname = fullfile(dir1,['HDM_' HDM_str '.mat']);
                                        IOI.HDM{s1,r1}.HDMfname = HDMfname;
                                        if ~O.only_display
                                            HDM0 = [];
                                            %Various options
                                            HDM0.O = O;
                                            HDM0.DO = DO;
                                            HDM0.EM = EM;
                                            HDM0.dir1 = dir1;
                                            HDM0.HDM_str = HDM_str;
                                            %Simulation option
                                            HDM0.S = S;
                                            %onsets: they are now specified in the module create_onsets
                                            name = IOI.sess_res{s1}.names{1};
                                            ons = IOI.sess_res{s1}.onsets{1};
                                            dur = IOI.sess_res{s1}.durations{1};
                                            if O.use_onset_amplitudes
                                                amp = IOI.sess_res{s1}.parameters{1};
                                                %normalize -- but what if amp takes extreme values?
                                                amp = amp/mean(amp);
                                            else
                                                amp = [];
                                            end
                                            bases.hrf.derivs = [0 0]; %not used
                                            [dummyX U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,1); %call only to get U                                           
                                            U.u = U.u(33:end); %?
                                            HDM0.U = U;
                                            HDM0.dt = IOI.dev.TR; 
                                            HDM0.N = 20/IOI.dev.TR; %too large?
                                            %Filters
                                            HDM0.HPF = HPF; %High pass filter on data
                                            HDM0.LPF = LPF;
                                            %data specification - which modalities to include:
                                            HDM0=ioi_get_data(ROI,HDM0,r1,s1);
                                            HDM0=ioi_set_physiomodel(HDM0);
                                            %choose priors
                                            HDM0=ioi_set_priors(HDM0);
                                            %scaling factors on covariance ... do we need that?
                                            %HDM0.pC = 10*HDM0.pC;
                                            %HDM0.pC(end,end) = 10*HDM0.pC(end,end);
                                            if S.simuOn
                                                HDM0 = ioi_simu_gen_parameters(HDM0);
                                                simuIt = HDM0.S.simuIt; 
                                                HDM0.Y0 = HDM0.Y; %background
                                                HDM0.HDM_str0 = HDM0.HDM_str; %save each figure of fit
                                                SHDM = []; %simulation HDM structure
                                            else                                              
                                                simuIt = 1;
                                            end
                                            warning('off','MATLAB:nearlySingularMatrix');
                                            warning('off','MATLAB:singularMatrix');
                                            % big loop over simulation iterations
                                            %--------------------------------------------------------------------------
                                            for it1=1:simuIt
                                                if S.simuOn
                                                    HDM0 = ioi_set_simu(HDM0,it1);
                                                end
                                                % nonlinear system identification
                                                %--------------------------------------------------------------------------                                               
                                                HDM0 = ioi_nlsi(HDM0);
                                                if S.simuOn
                                                    SHDM{it1} = ioi_save_simu(HDM0,it1); 
                                                end
                                            end
                                            warning('on','MATLAB:nearlySingularMatrix');
                                            warning('on','MATLAB:singularMatrix');
                                        else
                                            try
                                                load(HDMfname);
                                                HDM0 = HDM{r1,s1};
                                            catch
                                                disp(['HDM.mat not found for Session ' int2str(s1) ' and ROI ' int2str(r1)]);
                                            end
                                        end
                                        if ~S.simuOn
                                            HDM{r1,s1} = HDM0;
                                            save(HDMfname,'HDM');
                                            ioi_HDM_display(HDM0);
                                          
                                            %Store some of the information in IOI structure
                                            IOI.HDM{s1,r1}.Ep = HDM0.Ep;
                                            IOI.HDM{s1,r1}.Cp = HDM0.Cp;
                                            IOI.HDM{s1,r1}.K1 = HDM0.K1;
                                            IOI.HDM{s1,r1}.H1 = HDM0.H1;
                                            
                                            IOI.HDM{s1,r1}.F  = HDM0.F;
                                        else
                                            save(HDMfname,'SHDM');
                                        end
                                    end
                                end
                                disp(['HDM for session ' int2str(s1) ' completed']);
                                
        end
    end
end
% end
%                                 
%                             
%                             % Set model
%                             %--------------------------------------------------------------------------
%                             M.modal = modal;
%                             M.Model_Choice = Model_Choice;
%                             M.N     = 64;
%                             M.dt    = 24/M.N; %24/M.N;
%                             M.TE    = TE;
%                             M.pE    = pE;
%                             M.pC    = pC;
%                             M.m     = m;
%                             M.EM    = EM;
%                             if generate_figures || save_figures
%                                 M.nograph = 0;
%                             else
%                                 M.nograph = 1;
%                             end
%                             
%                             M = set_model(M);
%                             
%                             if ~S.simuOn
%                                 S.simuIt = 1;
%                             else
%                                 Y0 = Y;
%                                 if S.simuS == 0 % Stimuli type(s) to includes
%                                     S.simuS = 1:size(U.u,2);
%                                 end
%                             end
%                             
%                             %--------------------------------------------------------------------------
%                             warning('off','MATLAB:nearlySingularMatrix');
%                             warning('off','MATLAB:singularMatrix');
%                             
%                             for it1=1:S.simuIt
%                                 try
%                                     if S.simuOn
%                                         Y = Y0;
%                                         [Y,P] = simulate_data(pA,modal,M,U,Y,S,it1);
%                                     end
%                                     
%                                     
%                                     % INVERT MODEL %
%                                     % nonlinear system identification
%                                     %--------------------------------------------------------------------------
%                                     hfigevolution = []; hfigevolution2 = [];
%                                     if save_figures
%                                         fullfigDir1 = fullfigDir;
%                                     else
%                                         fullfigDir1 = '';
%                                     end
%                                     
%                                     % Plot data for debugging
%                                     tt = (0:size(Y.y,1)-1)*(Y.dt);
%                                     hfigData = figure('Units','normalized','Position',[0.1, 0.35, 0.6, 0.5]);
%                                     plot(tt,Y.y(:,1)*100,'.-b')
%                                     %ylim([-0.3 1.8])
%                                     xlim([0 200])
%                                     xlabel('Time (s)')
%                                     switch modal
%                                         case 1
%                                             ylabel('BOLD (%)')
%                                             ylim(100*[min(Y.y(:)) max(Y.y(:))]);
%                                         case 3
%                                             ylabel('Flow (%)')
%                                             ylim(100*[min(Y.y(:)) max(Y.y(:))]);
%                                         case {2,5}
%                                             hold on, plot(tt,Y.y(:,2)*100,'.-r')
%                                             %ylim([-3 18])
%                                             ylabel('ASL (BOLD & flow) (%)')
%                                             ylim(100*[min(Y.y(:)) max(Y.y(:))]);
%                                     end
%                                     try
%                                         title(['Upsampling: ' num2str(S.simuUpsample) ' ; Interpolation: ' num2str(S.simuInterp)])
%                                     end
%                                     
%                                     if save_figures
%                                         if ~exist(fullfigDir1,'dir'), mkdir(fullfigDir1); end
%                                         saveas(hfigData,fullfile(fullfigDir1,['Yy_' gen_num_str(it1,3)]),'fig');
%                                         print(hfigData, '-dtiffn', fullfile(fullfigDir1,['Yy_' gen_num_str(it1,3)]));
%                                     end
%                                     close(hfigData);
%                                     
%                                     if S.simuOn
%                                         SubjIdx0 = SubjIdx; %should not be used
%                                         SubjIdx = it1; %changed inside a for loop but it is restored later
%                                     end
%                                     %[Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F,hfigevolution,hfigevolution2] = nirs_nlsi(M,U,Y);
%                                     [Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2,F] = nirs_nlsi(M,U,Y,fullfigDir1,SubjIdx);
%                                     %for simulations: store results in place of subjects
%                                     
%                                     if S.simuOn
%                                         HDM{SubjIdx,r1}{s1}.EpS = P;
%                                     end
%                                     %Store results
%                                     if it1 == 1
%                                         %information independent of simulation iteration number
%                                         HDM{SubjIdx,r1}{s1}.M = M;
%                                         HDM{SubjIdx,r1}{s1}.pE = pE;
%                                         HDM{SubjIdx,r1}{s1}.job = job;
%                                     end
%                                     HDM{SubjIdx,r1}{s1}.Ep = Ep;
%                                     HDM{SubjIdx,r1}{s1}.Cp = Cp;
%                                     HDM{SubjIdx,r1}{s1}.F = F;
%                                     HDM{SubjIdx,r1}{s1}.name = ROIs(r1).nameROI;
%                                     HDM{SubjIdx,r1}{s1}.subj = fBOLD;
%                                     HDM{SubjIdx,r1}{s1}.session = cs1;
%                                     HDM{SubjIdx,r1}{s1}.Y = Y;
%                                     HDM{SubjIdx,r1}{s1}.U = U;
%                                     
%                                     % Volterra kernels of states
%                                     %--------------------------------------------------------------------------
%                                     [dummy,H1] = spm_kernels(M0,M1,M.N,M.dt);
%                                     HDM{SubjIdx,r1}{s1}.H1 = H1;
%                                     %save - overwriting previous
%                                     save(HDMfile,'HDM');
%                                     
%                                     %-display results
%                                     %==========================================================================
%                                     if generate_figures || save_figures
%                                         display_results(fullfigDir,Model_Choice,pE,Ep,Cp,U,m,M,H1,K1,K2,...
%                                             modal,plot_algebraic_CMRO2,save_figures,SubjIdx,hfigevolution,hfigevolution2);
%                                     end
%                                     
%                                     if S.simuOn
%                                         SubjIdx = SubjIdx0; % not used...
%                                     end
%                                     
%                                 catch exception
%                                     disp(['Model estimation or Simu #' int2str(it1) ' failed']);
%                                     disp(exception.identifier);
%                                     disp(exception.stack(1));
%                                 end
%                                 
%                             end
%                             
%                             out = [];
%                             
%                         catch exception
%                             disp(exception.identifier);
%                             disp(exception.stack(1));
%                             out = [0];
%                         end
%                     end
%                 end
%             end
%             
%             catch exception
%                 disp(exception.identifier);
%                 disp(exception.stack(1));
%         end
%     end
%     
    