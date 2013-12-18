function out = nirs_run_liom_tpHRF(job)
%==========================================================================
%Ke Peng, LIOM, Polytechnique
%2013-12-17, version 0.1, Function created
%==========================================================================

for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
        Topodatafile = NIRS.Dt.ana.rend;
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        sess = job.session_select;
        onsets = NIRS.Dt.fir.Sess(sess).U;
        [new_path,name,ext] = fileparts(newNIRSlocation);
        
        %******************************************************************
        %Topodata
        load(Topodatafile);
        brain_view = job.brain_view;
        rchn = rendered_MNI{brain_view}.rchn;
        cchn = rendered_MNI{brain_view}.cchn;
        %Two options
        W.AllowExtrapolation = job.AllowExtrapolation; %Option: 0: do not extrapolate
        W.no_interpolation = job.no_interpolation; %Option: 1: do not interpolate
        
        %find channels which are visible from this projection view
        W.index_ch = find(rchn ~= -1);
        %rendering surface
        brain = rendered_MNI{brain_view}.ren;
        W.s1 = size(brain, 1);
        W.s2 = size(brain, 2);
        %split into HbO and HbR interpolations
        W.ch = W.index_ch;
        W.rchn = rchn(W.index_ch);
        W.cchn = cchn(W.index_ch);
        %Generate interpolation matrix
        Q = interpolation_kernel_cine_simplified(W);
        %******************************************************************
        
        %******************************************************************
        %Time-course
        Data_dir = [NIRS.Dt.s.p '\dataSPM\coreg\Stat'];
        SessionTimeCourse = fopen_NIR(fullfile(Data_dir, ['Sess' int2str(sess) '.nir']), NC*3/2); %Load all time course of the session
        sample_interval = job.sample_interval;
        chromophore = job.chromophore_select;
        if isfield(job.HRFmethod, 'simple_averaging')
            time_before_spk = job.HRFmethod.simple_averaging.HRFstart;
            time_after_spk = job.HRFmethod.simple_averaging.HRFend;
            EvtofInterest = job.HRFmethod.simple_averaging.EvtInterest;
            if isfield(job.HRFmethod.simple_averaging.TCavoidance, 'avoidance_on')
                av.enabled = 1;
                av.EvtofNonInterest = job.HRFmethod.simple_averaging.TCavoidance.avoidance_on.EvtNonInterest;
            else
                av.enabled = 0;
            end
            estiHRF = average_HRF_timecourse(SessionTimeCourse, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, chromophore, av);
        end
        %******************************************************************
        
        %******************************************************************
        %Interpolation
        estiHRF = estiHRF(W.ch,:);
        tp_interval = round(sample_interval*fs);
        cp = 1;
        i0 = 0;
        while cp <= size(estiHRF,2)
            amplitude = estiHRF(:,cp)';
            amplitude_i = amplitude * Q.B;
            amplitude_i_reshaped = zeros(W.s1,W.s2);
            %amplitude_i_reshaped = reshape(amplitude_i,W.s1,W.s2);
            amplitude_i_reshaped(Q.index_mask) = amplitude_i(1,:);
            c = max(max(abs(estiHRF)));
            i0 = i0+1;
            h = figure;
            imagesc(amplitude_i_reshaped);
            colorbar
            caxis([-c c]);
            %Write markers
            if isfield(job.HRFmethod, 'simple_averaging')
                mk = [int2str(time_before_spk + (i0-1)*sample_interval) ' s'];
                text(300, 450, mk, 'Color', 'b', 'FontSize',18);
            end
            %save
            switch chromophore
                case 0
                    filestr = [new_path '\HbO_' int2str(i0)];
                case 1
                    filestr = [new_path '\HbR_' int2str(i0)];
                case 2
                    filestr = [new_path '\HbT_' int2str(i0)];
            end
            saveas(gcf,filestr,'png');
            pause(1);
            close(h);
            cp = cp + tp_interval;
        end
        %******************************************************************
        
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1));
        disp(['Could not estimate HRF for subject' int2str(iSubj) ' for ' job.NIRSmat{iSubj,1}]);
    end
end
outNIRSmat = NIRS;
out.NIRSmat = outNIRSmat;