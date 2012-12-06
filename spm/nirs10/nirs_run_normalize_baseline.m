function out = nirs_run_normalize_baseline(job)
% Filename prefix for normalized files
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
            % Option to remove negative values - such channels are probably too noisy to be
            % useful anyway - but a better method should be found
            if strcmp(NIRS.Cf.dev.n,'ISS Imagent')
                legacy_option_to_remove_negative_values = 1;
            else
                legacy_option_to_remove_negative_values = 0;
            end
            % For each data file - loop over sessions
            for f=1:size(rDtp,1)
                % Read raw data (intensity)
                d = fopen_NIR(rDtp{f,1},NC);
%                 test_clear_protocol = 1;
%                 if test_clear_protocol
%                      %protocole 
%                      d_max = max(d,[],2);
%                      lp1 = [linspace(60*25,120*25-1,60*25) linspace(300*25,360*25-1,60*25)];
%                      d(29:56,lp1) = d(29:56,lp1) - repmat(d_max(29:56)/4,1,120*25);
%                      d(1:28,lp1) = d(1:28,lp1) + repmat(d_max(1:28)/4,1,120*25);
%                 end
                % Option to remove negative values (to improve!)
                if legacy_option_to_remove_negative_values
                    threshold = 0.1;
                    % Some values of the optical intensity d may be negative
                    % Compute the minimum of d for each channel
                    mind = min(d,[],2);
                    for i1=1:size(mind,1)
                        if mind(i1) > threshold
                            mind(i1) = 0;
                        else
                            NIRS.CAUTION = 'legacy_option_to_remove_negative_values';
                            if mind(i1) < threshold
                                mind(i1) = -mind(i1)+threshold;
                            end
                        end
                    end
                    % Subtract twice this minimum value for regularization - in
                    % anticipation of taking the log later
                    d = d + mind * ones(1,size(d,2));
                else
                    %reset negative values to the minimal positive value
                    dMin = min(d(d>0));
                    nL0 = sum(d(:)<0);
                    d(d<=0) = dMin;
                end
                
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
