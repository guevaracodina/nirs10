function out = nirs_run_HDM(job)
% User interface for hemodynamic model estimation -- based on SPM HDM
Stimuli = job.which_condition;
S = nirs_get_simu(job);
EM = nirs_get_EM(job);
O = nirs_get_common_O_HDM_SCKS(job);
DO = nirs_get_common_hdm_display_options(job);
IC = nirs_get_colors(job);
[all_sessions selected_sessions] = nirs_get_sessions(job);
[all_channels selected_channels] = nirs_get_channels(job);

% Loop over subjects
for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        clear NIRS HDM0
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'HDM_OK') || job.force_redo)
            if DO.generate_figures || DO.save_figures
                spm('Pointer','Watch')
                spm('FigName','Estimation in progress');
            end
            Sess = NIRS.Dt.fir.Sess;
            dt = 1/NIRS.Cf.dev.fs;
            NC = NIRS.Cf.H.C.N; %number of HbO and HbR channels
            [dir1 fil1] = fileparts(newNIRSlocation);
            HDMfname = fullfile(dir1,'HDM.mat');
            NIRS.HDM.HDMfname = HDMfname;
            %set up common HDM structure
            HDM0 = nirs_initial_hdm_setup(O,DO,EM,IC,S,dt,dir1);
            HDM0 = nirs_set_physiomodel(HDM0);
            %choose priors
            HDM0 = nirs_set_hdm_priors(HDM0);
            if S.simuOn
                HDM0 = nirs_simu_gen_parameters(HDM0);
                simuIt = HDM0.S.simuIt;
            else
                simuIt = 1;
            end
            %loop over each session
            for s1=1:length(Sess)
                if all_sessions || sum(s1==selected_sessions)
                    %get data - selected channels
                    dname = NIRS.Dt.fir.pp(end).p{s1};
                    d = fopen_NIR(dname,NC);
                    %get onsets in right format for HDM
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
                            Y.X0 = []; %confounds                            
                            
                            if ~DO.only_display
                                HDM = HDM0;
                                HDM.U = U;
                                HDM.HDM_str = ['S' gen_num_str(s1,2) 'C' gen_num_str(c1,3)];                                                   
                                %data specification - which modalities to include:
                                HDM=nirs_get_hdm_data(HDM,Y,c1,s1);                               
                                %scaling factors on covariance ... do we need that?
                                %HDM0.pC = 10*HDM0.pC;
                                %HDM0.pC(end,end) = 10*HDM0.pC(end,end);
                                if S.simuOn
                                    HDM.Y0 = HDM.Y; %background
                                    HDM.HDM_str0 = HDM.HDM_str; %save each figure of fit
                                    SHDM = []; %simulation HDM structure
                                end
                                warning('off','MATLAB:nearlySingularMatrix');
                                warning('off','MATLAB:singularMatrix');
                                % big loop over simulation iterations
                                for it1=1:simuIt
                                    if S.simuOn
                                        HDM = nirs_set_simu(HDM,it1);
                                    end
                                    % nonlinear system identification
                                    HDM = nirs_nlsi_V2(HDM);
                                    if S.simuOn
                                        SHDM{it1} = nirs_save_simu(HDM,it1);
                                    end
                                end
                                warning('on','MATLAB:nearlySingularMatrix');
                                warning('on','MATLAB:singularMatrix');
                            else
                                try
                                    load(HDMfname);
                                    HDM = HDM{s1,c1};
                                catch
                                    disp(['HDM.mat not found for Session ' int2str(s1) ' and channel ' int2str(c1)]);
                                end
                            end
                            if ~S.simuOn
                                HDM{s1,c1} = HDM;
                                save(HDMfname,'HDM');
                                nirs_HDM_display(HDM);                            
                                %Store some of the information in NIRS structure
                                NIRS.HDM{s1,c1}.Ep = HDM.Ep;
                                NIRS.HDM{s1,c1}.Cp = HDM.Cp;
                                NIRS.HDM{s1,c1}.K1 = HDM.K1;
                                NIRS.HDM{s1,c1}.H1 = HDM.H1;                                
                                NIRS.HDM{s1,c1}.F  = HDM.F;
                            else
                                save(HDMfname,'SHDM');
                            end
                        end
                    end
                end
            end
            disp(['HDM for session ' int2str(s1) ' completed']);
            NIRS.flags.HDM_OK = 1;
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
    end
end
out.NIRSmat = job.NIRSmat;

