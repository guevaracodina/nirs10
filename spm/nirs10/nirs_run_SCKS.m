function out = nirs_run_SCKS(job)
% User interface for hemodynamic model estimation -- based on SPM HDM
Stimuli = job.which_condition;
S = nirs_get_simu_SCKS(job);
O = nirs_get_common_O_HDM_SCKS(job);
DO = nirs_get_common_SCKS_display_options(job);
IC = nirs_get_colors(job);
[all_sessions selected_sessions] = nirs_get_sessions(job);
[all_channels selected_channels] = nirs_get_channels(job);

% Loop over subjects
for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        clear NIRS SCKS
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && ((~isfield(NIRS,'flags') || ~isfield(NIRS.flags,'SCKS_OK')) || job.force_redo)
            if DO.generate_figures || DO.save_figures
                spm('Pointer','Watch')
                spm('FigName','Estimation in progress');
            end
            Sess = NIRS.Dt.fir.Sess;
            dt = 1/NIRS.Cf.dev.fs;
            NC = NIRS.Cf.H.C.N; %number of HbO and HbR channels
            [dir1 fil1] = fileparts(newNIRSlocation);
            SCKSfname = fullfile(dir1,'SCKS.mat');                        
            NIRS.SCKS.SCKSfname = SCKSfname;
            %loop over each session
            for s1=1:length(Sess)
                if all_sessions || sum(s1==selected_sessions)
                    %get data - selected channels
                    dname = NIRS.Dt.fir.pp(end).p{s1};
                    d = fopen_NIR(dname,NC);
                    %get onsets in right format for SCKS?
                    U = Sess(s1).U(Stimuli);
                    if O.use_onset_amplitudes
                        U.amp = NIRS.sess_res{s1}.parameters{1}; %to be calculated
                        %normalize -- but what if amp takes extreme values?
                        U.amp = U.amp/mean(amp);
                    else
                        U.amp = [];
                    end
                    U = nirs_get_U(dt,size(d,2),U); %call only to get U
                    %Loop over channels -- could this be done in parallel?
                    for c1=1:NC/2
                        if all_channels || sum(c1==selected_channels)
                            Y.y = [];
                            if IC.include_HbR
                                Y.y(:,1) = d(c1+NC/2,:)'; %HbR
                            end
                            if IC.include_HbT || IC.include_HbO
                                Y.y = [Y.y d(c1,:)'+d(c1+NC/2,:)'];
                            end
                            Y.dt = dt;
                            Y = nirs_SCKS_filter(Y,O.LPF,O.HPF); %creates Y.X0
                            
                            if ~DO.only_display
                                %SCKS on each channel
                                SCKS0 = nirs_initial_SCKS_setup(O,DO,IC,S,dt,dir1);
                                SCKS0.subj_id = NIRS.Dt.s.subj_id;
                                SCKS0.s1 = s1; %session
                                SCKS0.c1 = c1; %channel
                                SCKS0.SCKSparams = job.SCKSparams;
                                SCKS0.HPF = O.HPF; %High pass filter on data
                                SCKS0.LPF = O.LPF;
                                SCKS0.Y = Y;
                                %SCKS0 = SCKS_get_data(SCKS0,Y);
                                SCKS0 = nirs_set_physiomodel(SCKS0);
                                SCKS0 = nirs_set_SCKS_priors(SCKS0);
                                if S.simuOn
                                    SCKS0 = nirs_simu_gen_parameters_SCKS(SCKS0);
                                    simuIt = SCKS0.S.simuIt;
                                else
                                    simuIt = 1;
                                end
                                SCKS0 = nirs_SCKS_set_SCKS(SCKS0,U);
                                for it1=1:simuIt
                                    if S.simuOn
                                        HDM0 = nirs_set_simu_SCKS(HDM0,it1); %not coded up yet
                                    end
                                    
                                    if SCKS0.SCKSparams.SCKSnoise
                                        SCKS0 = nirs_SCKS(SCKS0,DO.generate_figures); %noise and parameters promoted to states
                                    else
                                        SCKS0 = nirs_SCKS2(SCKS0,DO.generate_figures); %reduced version, noise and parameters are time independent
                                    end
                                    if S.simuOn
                                        SSCKS{it1} = nirs_save_simu_SCKS(SCKS0,it1); %not coded up yet
                                        disp(['SCKS for simulation ' int2str(it1) ' completed']);
                                        %must think about saving SSCKS
                                    end
                                end
                                %Store and save results obtained so far
                                SCKS{c1,s1} = SCKS0;
                                save(SCKSfname,'SCKS');
                            else
                                try
                                    if ~exist('SCKS','var')
                                        load(SCKSfname);
                                    end
                                    SCKS0 = SCKS{c1,s1};
                                catch
                                    disp(['SCKS.mat not found for Session ' int2str(s1) ' and channel ' int2str(c1)]);
                                end
                            end
                            
                            if ~S.simuOn
%                                 SCKS{s1,c1} = SCKS;
%                                 save(SCKSfname,'SCKS');
                                nirs_SCKS_display(SCKS);                            
                                %Store some of the information in NIRS structure
                                %NIRS.SCKS{s1,c1}.Ep = SCKS.Ep;
%                                 NIRS.SCKS{s1,c1}.Cp = SCKS.Cp;
%                                 NIRS.SCKS{s1,c1}.K1 = SCKS.K1;
%                                 NIRS.SCKS{s1,c1}.H1 = SCKS.H1;                                
%                                 NIRS.SCKS{s1,c1}.F  = SCKS.F;
%                             else
%                                 save(SCKSfname,'SCKS');
                            end
                            disp(['SCKS for channel ' int2str(c1) ' completed']);
                        end
                    end
                    disp(['SCKS for session ' int2str(s1) ' completed']);
                end                
            end
            NIRS.flags.SCKS_OK = 1;
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
    end
end
out.NIRSmat = job.NIRSmat;

