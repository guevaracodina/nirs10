function out = nirs_run_normalize_baseline(job)
%filename prefix 
prefix = 'b'; %for "baseline"
DelPreviousData  = job.DelPreviousData;
try 
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end
bl_m = job.Normalize_OD;
group_consecutive_markers = 0; %boolean
add_or_mult = job.add_or_mult;

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
           
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        
        
        for f=1:size(rDtp,1)
            d = fopen_NIR(rDtp{f,1},NC);
            
            legacy_option_to_remove_negative_values = 1;
            if legacy_option_to_remove_negative_values
                threshold = 0.1;
                %Some values of the optical intensity d may be negative
                %Compute the minimum of d for each channel
                mind = min(d,[],2);
                for i1=1:size(mind,1)
                    if mind(i1) > threshold 
                        mind(i1) = 0;
                    else 
                        if mind(i1) > -threshold
                             mind(i1) = -threshold;
                        end
                    end
                end
                %Subtract twice this minimum value for regularization - in
                %anticipation of taking the log later
                d = d - 2 *  mind * ones(1,size(d,2));
            end
            
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
            if markers_available && job.normalization_type == 2
                %temporary for data
                td = zeros(size(d));
                %group consecutive markers, as an option
                if group_consecutive_markers
                    [bpi bpd] = find_marker_durations(bpi);
                end
                for i=1:length(si)
                    try
                        e = d(:,si(i):ei(i));
                    catch
                        e = [];
                    end
                    
                    try
                        %Normalization factor
                        switch bl_m
                            case 0 %0: Median;
                                div_factor = median(e,2);
                            case 1 %1: Initial value of that interval
                                div_factor = e(:,1);
                            case 2 %mean
                                div_factor = mean(e,2);
                            otherwise %take median
                                div_factor = median(e,2);
                        end
                    catch
                        div_factor = 1;
                    end
                    %normalize
                    try
                        if ~isempty(e)
                            div_factor = div_factor*ones(1,size(e,2));

                            if add_or_mult
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
                                td(:,si(i):ei(i)) = (75+(e(ch,:)-div_factor))*job.Analyzer_sf;
                                %HbR channels
                                ch = NIRS.Cf.H.C.wl== HbR_like;
                                td(:,si(i):ei(i)) = (25+(e(ch,:)-div_factor))*job.Analyzer_sf;
                            else
                                td(:,si(i):ei(i)) = e./div_factor*job.Analyzer_sf; 
                            end 
                        end
                    catch
                        %do nothing
                    end
                                        
                end
                %replace d - this sets intervals with movement to zero
                d=td;
            else
                if job.normalization_type == 3
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
                        try
                            %Normalization factor
                            switch bl_m
                                case 0 %0: Median;
                                    div_factor = median(d(:,win),2);
                                case 1 %1: Initial value - not too sensible - better last value
                                    div_factor = d(:,win(end));
                                case 2 %mean
                                    div_factor = mean(d(:,win),2);
                                otherwise %take median
                                    div_factor = median(d(:,win),2);
                            end
                        catch
                            div_factor = 1;
                        end
                        %Normalize to the next stimulus only
                        if i1 < length(onsets)
                            wine2 = round(onsets(i1+1)*fs)-1;                           
                        else
                            wine2 = size(d,2);                           
                        end
                        win2 = wins:wine2;
                        div_factor = div_factor * ones(1,length(win2));
                    end
                else
                    %global normalization either requested or because
                    %no markers available - then normalize whole series
                    bpi = [];
                    bpd = [];

                    try
                        %Normalization factor
                        switch bl_m
                            case 0 %0: Median;
                                div_factor = median(d,2);
                            case 1 %1: Initial value
                                div_factor = d(:,1);
                            case 2 %mean
                                div_factor = mean(d,2);
                            otherwise %take median
                                div_factor = median(d,2);
                        end
                    catch
                        div_factor = 1;
                    end
                    div_factor = div_factor * ones(1,size(d,2));
                    %normalize
                    if add_or_mult
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
                        ch = find(NIRS.Cf.H.C.wl== HbO_like);
                        d(ch,:) = (75+(d(ch,:)-div_factor))*job.Analyzer_sf;
                        %HbR channels
                        ch = find(NIRS.Cf.H.C.wl== HbR_like);
                        d(ch,:) = (25+(d(ch,:)-div_factor))*job.Analyzer_sf; 
                    else
                        d = d./div_factor*job.Analyzer_sf; 
                    end
                end
            end
           
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            if NewDirCopyNIRS
                dir2 = [dir1 filesep NewNIRSdir];
                if ~exist(dir2,'dir'), mkdir(dir2); end; 
                outfile = fullfile(dir2,[prefix fil1 ext1]);
            else
                outfile = fullfile(dir1,[prefix fil1 ext1]);
            end
            if DelPreviousData
                delete(rDtp{f,1});
            end
            fwrite_NIR(outfile,d);
            %add outfile name to NIRS
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
        if NewDirCopyNIRS
            save(fullfile(dir2,'NIRS.mat'),'NIRS');            
        else
            save(job.NIRSmat{Idx,1},'NIRS'); 
        end 
    catch
        disp(['Could not normalize to baseline for subject' int2str(Idx)]);
    end       
end
out.NIRSmat = job.NIRSmat;