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
if isfield(job.target_sampling_rate,'specified_sampling_rate')
    target_sampling_rate = 1;
    sampling_rate = job.target_sampling_rate.specified_sampling_rate.sampling_rate;
else
    target_sampling_rate = 0;
end
DO.verbose = 0;
% Loop over subjects
for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        %Various HDM structures
        %HDM0: generic HDM template, to use for all sessions and channels
        %HDM: for results specific to one channel and session
        %HDMc: 2D cell of HDM structures, to be saved
        %SHDM: HDM for simulations
        clear NIRS HDM0
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && ((~isfield(NIRS,'flags') || ~isfield(NIRS.flags,'HDM_OK')) || job.force_redo)
            if DO.generate_figures || DO.save_figures
                spm('Pointer','Watch')
                spm('FigName','Estimation in progress');
            end
            Sess = NIRS.Dt.fir.Sess;
            dt = 1/NIRS.Cf.dev.fs;
            dt0= dt;
            if target_sampling_rate
                decimate_factor = round(1/(sampling_rate*dt));
                dt = 1/sampling_rate;
            end
            NC = NIRS.Cf.H.C.N; %number of HbO and HbR channels
            [dir1 fil1] = fileparts(newNIRSlocation);
            HDMfname = fullfile(dir1,'HDM.mat');
            NIRS.HDM.HDMfname = HDMfname;
            %set up common HDM structure
            HDM0 = nirs_initial_hdm_setup(O,DO,EM,IC,S,dt,dir1);
            HDM0.subj_id = NIRS.Dt.s.subj_id;
            HDM0.dir1 = dir1;
            HDM0.N = round(EM.kernel_window/dt);
            HDM0 = nirs_set_physiomodel(HDM0);
            %choose priors
            HDM0 = nirs_set_hdm_priors(HDM0);
            if S.simuOn
                HDM0 = nirs_simu_gen_parameters(HDM0);
                simuIt = HDM0.S.simuIt;
            else
                simuIt = 1;
            end
            HDMloaded = 0;
            %loop over each session
            for s1=1:length(Sess)
                if all_sessions || sum(s1==selected_sessions)
                    %get data - selected channels
                    dname = NIRS.Dt.fir.pp(end).p{s1};
                    d = fopen_NIR(dname,NC);   
                    if O.LPF.lpf_gauss_On || target_sampling_rate
                        K = get_K(1:size(d,2),O.LPF.fwhm1,dt0);
                        d = nirs_filter_HPF_LPF_WMDL(K,d')';
                    end                   
                    if target_sampling_rate
                        %downsample d
                        d = d(:,round(decimate_factor/2):decimate_factor:end);
                    end
                    
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
                                HDM.c1 = c1;
                                HDM.s1 = s1;
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
%                                     clear Fin0 dataOLD
%                                     global Fin0 dataOLD
                                    HDM = nirs_nlsi_V2(HDM);
                                    if S.simuOn
                                        SHDM{it1} = nirs_save_simu(HDM,it1);
                                    end
                                end
                                warning('on','MATLAB:nearlySingularMatrix');
                                warning('on','MATLAB:singularMatrix');
                                HDM.HDM_OK = 1;
                            else
                                try
                                    if ~HDMloaded
                                        load(HDMfname);
                                        HDMloaded = 1;
                                    end
                                    HDM = HDMc{s1,c1};
                                catch
                                    disp(['HDM.mat not found for Session ' int2str(s1) ' and channel ' int2str(c1)]);
                                end
                            end
                            if ~S.simuOn
                                HDMc{s1,c1} = HDM;
                                save(HDMfname,'HDMc');
                                nirs_HDM_display(HDM);
                                %Store some of the information in NIRS structure
                                %                                 NIRS.HDMc{s1,c1}.Ep = HDM.Ep;
                                %                                 NIRS.HDMc{s1,c1}.Cp = HDM.Cp;
                                %                                 NIRS.HDMc{s1,c1}.K1 = HDM.K1;
                                %                                 NIRS.HDMc{s1,c1}.H1 = HDM.H1;
                                %                                 NIRS.HDMc{s1,c1}.F  = HDM.F;
                            else
                                save(HDMfname,'SHDM');
                            end
                            disp(['HDM for channel ' int2str(c1) ' completed']);
                        end
                    end
                    disp(['HDM for session ' int2str(s1) ' completed']);
                end
            end
            NIRS.flags.HDM_OK = 1;
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
    end
end
out.NIRSmat = job.NIRSmat;

