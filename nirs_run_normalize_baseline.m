function out = nirs_run_normalize_baseline(job)
%Filename prefix for normalized files
prefix = 'b'; % for "baseline"
% User-specified options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bl_m = job.Normalize_OD; % {0='Median',1='Initial Value',2='Mean'};
group_consecutive_markers = 0; %boolean
add_or_mult = job.add_or_mult; % {1='Additive', 0='Multiplicative'} normalization
normalization_type = job.normalization_type; % {1='Global', 2='By bad point segments', 3='By stimuli'};
%Scaling factor, useful to rescale data, for example to view it in Brain
%Vision Analyzer 2
scaling_factor = job.Analyzer_sf; %job.scaling_factor; % Scaling factor applied to the amplitude of all channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div_factor0 = 1;
%**************************************************************************
%Introduce nirs data jump filling process
%Ke Peng
%**************************************************************************
%Should be done before normalization of data
if isfield(job, 'nirs_filling_jumps') && isfield(job.nirs_filling_jumps, 'nirs_fill_jumps_on')
    fill_jump_on = 1;
    num_standard_deviation = job.nirs_filling_jumps.nirs_fill_jumps_on.num_standard_deviation;
    num_points = job.nirs_filling_jumps.nirs_fill_jumps_on.num_points;
    size_gap = job.nirs_filling_jumps.nirs_fill_jumps_on.size_gap;
    
    if isfield(job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable, 'HPF_enable_on')
        HPF_enable_on = 1;
        hpf_butter_order = job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable.HPF_enable_on.hpf_butter_order;
        hpf_butter_freq = job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable.HPF_enable_on.hpf_butter_freq;
    else
        HPF_enable_on = 0;
    end
else
    if isfield(job, 'nirs_filling_jumps') && isfield(job.nirs_filling_jumps, 'nirs_new_fill_jumps')
        fill_jump_on = 2;
    else
        fill_jump_on = 0;
    end
end

if isfield(job.filters,'filters_on')
    FiltersOn = 1;
    filters_on = job.filters.filters_on;
    %PCA
    PCA = filters_on.channel_pca;
    NumPCAComponents = filters_on.NumPCAComponents;
    %LPF
%     lpf0 = filters_on.lpf;
%     if isfield(lpf0,'lpf_gauss')
%         FWHM = lpf0.lpf_gauss.fwhm1;
%         lpf = 'gaussian';
%     else
%         if isfield(lpf0,'lpf_hrf')
%             lpf = 'hrf';
%         else
%             if isfield(lpf0,'lpf_none')
%                 lpf = 'none';
%             else
%                 disp('Unrecognized low pass filter');
%             end
%         end
%     end

    %LPF
    lpf0 = filters_on.lpf;
    if isfield(lpf0,'lpf_butter')
        lpf_butter_freq = lpf0.lpf_butter.lpf_butter_freq;
        lpf_butter_order = lpf0.lpf_butter.lpf_butter_order;
        lpf = 1;
    else
        lpf = 0; %no low pass filter
    end    
    %HPF
    hpf0 = filters_on.hpf;
    if isfield(hpf0,'hpf_butter_On')
        hpf_butter_freq = hpf0.hpf_butter_On.hpf_butter_freq;
        hpf_butter_order = hpf0.hpf_butter_On.hpf_butter_order;
        hpf = 1;
    else
        if isfield(job.hpf_butter,'remove_linear')
            hpf = 2;
        else
            if isfield(job.hpf_butter,'GLM_remove_linear')
                hpf = 3; %no high pass filter
            else
                if isfield(job.hpf_butter,'SPM_cosine_filter')
                    hpf = 4; %no high pass filter
                    hpf_cos_freq = 1/300; %1/240; %4 minute
                    hpf_cos_order = 3; %not used
                else
                    hpf = 0; %no high pass filter
                end
            end
        end
    end
else
    FiltersOn = 0;
end
%**************************************************************************


% Loop over subjects
for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'normbaseOK') || job.force_redo)
            % Perform normalization on output of the last step of preprocessing
            % that has been performed
            lst = length(NIRS.Dt.fir.pp);
            rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
            NC = NIRS.Cf.H.C.N; % Number of channels
            fs = NIRS.Cf.dev.fs; % Sampling frequency
            % For each data file - loop over sessions
            for f=1:size(rDtp,1)
                % Read raw data (intensity)
                d = fopen_NIR(rDtp{f,1},NC);
                %reset negative values to the minimal positive value: not a good idea
                %Instead: add twice the minimum value to all the channels
                %treat each channel separately
                minimal_correction = 1;
                if minimal_correction
                    correction1 = 1.01;
                    correction2 = 0.01;
                else
                    correction1 = 2;
                    correction2 = 1;
                end
                
                for c0=1:size(d,1)
                    td = d(c0,:);
                    dMin = min(td);
                    if dMin < 0
                        d(c0,:) = td - correction1*dMin;
                    else
                        if dMin == 0
                            dMin = min(td(td>0));
                            d(c0,:) = td + correction2*dMin;
                        end
                    end
                end
                
                %**************************************************************************
                %Filling jumps in nirs data
                %Ke Peng
                %**************************************************************************
                if job.outputdatafigures
                    nirs_time_plots(d,fs,NC,f,newNIRSlocation,'rm_Neg',{NIRS.Cf.dev.wl(1) NIRS.Cf.dev.wl(2)});
                end
                
                %disp(newNIRSlocation)
                switch fill_jump_on
                    case 2
                        d = nirs_new_remove_jumps(d,newNIRSlocation);
                        if job.outputdatafigures
                            nirs_time_plots(d,fs,NC,f,newNIRSlocation,'rm_jumps',{NIRS.Cf.dev.wl(1) NIRS.Cf.dev.wl(2)});
                        end
                    case 1
                        OP.Sb = num_standard_deviation;
                        OP.Nr = num_points;
                        OP.Mp = size_gap;
                        OP.sf = fs;
                        OP.ubf = 0;
                        for i0 = 1:size(d,1)
                            d(i0,:) = nirs_remove_jumps_new(d(i0,:),OP);
                        end
                    otherwise
                end
                
                
                %**************************************************************************
                %Filters
                if FiltersOn
                    %PCA
                    if PCA
                        d = d';
                        %process each wavelength separately
                        which_channels = 1:(NC/2); %all HbOchannels
                        Y1 = nirs_makePcaV2(d(:,which_channels),which_channels,nComponents);
                        Y2 = nirs_makePcaV2(d(:,which_channels+NC/2),which_channels,nComponents);
                        d = [Y1 Y2];
                        d = d';
                        %Output graphs of PCA result
                        if job.outputdatafigures
                            nirs_time_plots(d,fs,NC,f,newNIRSlocation,'PCA',{'WL1' 'WL2'});
                        end
                    end
                    %HPF
                    if hpf
                        %to add DC value back
                        DC = mean(d,2);
                    end
                    switch hpf
                        case 1
                            d = ButterHPF(fs,hpf_butter_freq,hpf_butter_order,d')';
                        case 2
                            Y = d';
                            nS = size(Y,1);
                            mX = linspace(0,round(nS/fs),nS);
                            mX = [mX' ones(nS,1)];
                            pmX = pinv(mX);
                            Y = Y - mX * (pmX * Y);
                            d = Y';
                    end
                    if hpf
                        d = d + repmat(DC,[1 size(d,2)]);
                    end
                    if hpf && job.outputdatafigures
                        nirs_time_plots(d,fs,NC,f,newNIRSlocation,'HPF',{'WL1' 'WL2'});
                    end
                    %LPF
                    if lpf
                        d = ButterLPF(fs,lpf_butter_freq,lpf_butter_order,d')';
                        if job.outputdatafigures
                            nirs_time_plots(d,fs,NC,f,newNIRSlocation,'LPF',{'WL1' 'WL2'});
                        end
                    end
                end
                %**************************************************************************
                
                % Read markers for movement if available
                try
                    bpi = NIRS.Dt.fir.pp(lst).bpi{f,1}; %bad point indices
                    bpd = NIRS.Dt.fir.pp(lst).bpd{f,1}; %bad point durations
                    si = NIRS.Dt.fir.pp(lst).si{f,1};
                    ei = NIRS.Dt.fir.pp(lst).ei{f,1};
                    if isempty(bpi)
                        markers_available = 0;
                    else
                        markers_available = 1;
                    end
                catch
                    markers_available = 0;
                end
                
                % Normalization by bad point segments
                if markers_available && normalization_type == 2
                    % temporary for data
                    td = zeros(size(d));
                    % group consecutive markers, as an option
                    if group_consecutive_markers
                        [bpi bpd] = find_marker_durations(bpi);
                    end
                    for i=1:length(si)
                        try
                            e = d(:,si(i):ei(i));
                        catch
                            e = [];
                        end
                        win = 1:size(e,2);
                        [div_factor div_factor0] = get_div_factor(e,fs,bl_m,f,win,div_factor0);
                        
                        % Perform normalization
                        try
                            if ~isempty(e)
                                div_factor = div_factor*ones(1,size(e,2));
                                
                                if add_or_mult % Additive
                                    %normalize median to 75 uM for HbO and 25 uM for HbR
                                    wl = NIRS.Cf.dev.wl;
                                    HbO_like = [];
                                    HbR_like = [];
                                    for i=1:length(wl)
                                        if wl(i) > 750 %in nanometer
                                            %found a wavelength that is HbO-like
                                            %Note code at present will work for only one HbO-like
                                            %wavelength
                                            HbO_like = [HbO_like i];
                                        else
                                            HbR_like = [HbR_like i];
                                        end
                                    end
                                    %HbO channels
                                    ch = NIRS.Cf.H.C.wl== HbO_like;
                                    td(:,si(i):ei(i)) = (75+(e(ch,:)-div_factor(ch,:)))*scaling_factor;
                                    %HbR channels
                                    ch = NIRS.Cf.H.C.wl== HbR_like;
                                    td(:,si(i):ei(i)) = (25+(e(ch,:)-div_factor(ch,:)))*scaling_factor;
                                else
                                    td(:,si(i):ei(i)) = e./div_factor*scaling_factor;
                                end
                            end
                        catch
                            %do nothing
                        end
                        
                    end
                    %replace d - this sets intervals with movement to zero
                    d=td;
                    
                else
                    % Normalization by stimuli segment
                    if normalization_type == 3
                        %Normalize using window prior to stimuli
                        baseline_duration = round(job.baseline_duration*fs);
                        %take first stimulus - could generalize to loop over all
                        %stimuli
                        onsets = NIRS.Dt.fir.Sess(f).U(1).ons;
                        %loop over onsets
                        for i1=1:length(onsets)
                            %find window to normalize over
                            %end of window
                            wine = round(onsets(i1)*fs);
                            if baseline_duration < wine - 1
                                %start of baseline
                                wins = wine - baseline_duration;
                            else
                                wins = 1;
                            end
                            win = wins:wine;
                            
                            [div_factor div_factor0] = get_div_factor(d,fs,bl_m,f,win,div_factor0);
                            
                            % Normalize to the next stimulus only
                            if i1 < length(onsets)
                                wine2 = round(onsets(i1+1)*fs)-1;
                            else
                                wine2 = size(d,2);
                            end
                            win2 = wins:wine2;
                            div_factor = div_factor * ones(1,length(win2));
                        end
                        
                        % Global normalization (all time series)
                    else
                        %global normalization either requested or because
                        %no markers available - then normalize whole series
                        bpi = [];
                        bpd = [];
                        win = 1:size(d,2);
                        [div_factor div_factor0] = get_div_factor(d,fs,bl_m,f,win,div_factor0);
                        
                        div_factor = div_factor * ones(1,size(d,2));
                        
                        % Perform normalization
                        if add_or_mult % "Additive"
                            % Normalize median/mean/initial value to 75 uM for HbO and 25 uM for HbR
                            wl = NIRS.Cf.dev.wl;
                            HbO_like = [];
                            HbR_like = [];
                            for i=1:length(wl)
                                if wl(i) > 750 %in nanometer
                                    %found a wavelength that is HbO-like
                                    %Note code at present will work for only one HbO-like
                                    %wavelength
                                    HbO_like = [HbO_like i];
                                else
                                    HbR_like = [HbR_like i];
                                end
                            end
                            %HbO channels
                            ch = find(NIRS.Cf.H.C.wl== HbO_like);
                            d(ch,:) = (75+(d(ch,:)-div_factor(ch,:)))*scaling_factor;
                            %HbR channels
                            ch = find(NIRS.Cf.H.C.wl== HbR_like);
                            d(ch,:) = (25+(d(ch,:)-div_factor(ch,:)))*scaling_factor;
                        else % "Multiplicative" (I/I0)
                            d = d./div_factor*scaling_factor;
                        end
                    end
                end
                
                [dir1,fil1,ext1] = fileparts(rDtp{f});
                outfile = fullfile(dir1,[prefix fil1 ext1]);
                
                fwrite_NIR(outfile,d);
                % add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'normalize_baseline';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
                try
                    NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
                    NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
                    NIRS.Dt.fir.pp(lst+1).si{f,1} = si;
                    NIRS.Dt.fir.pp(lst+1).ei{f,1} = ei;
                catch
                end
                if job.outputdatafigures
                    if add_or_mult
                        nirs_time_plots(d,fs,NC,f,newNIRSlocation,'norm_add',{'HbO' 'HbR'});
                    else
                        nirs_time_plots(d,fs,NC,f,newNIRSlocation,'norm_mult',{NIRS.Cf.dev.wl(1) NIRS.Cf.dev.wl(2)});
                        nirs_time_plots(1./d,fs,NC,f,newNIRSlocation,'inverse_norm_mult',{NIRS.Cf.dev.wl(1) NIRS.Cf.dev.wl(2)});
                    end
                end
            end
            NIRS.flags.normbaseOK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not normalize to baseline for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;
