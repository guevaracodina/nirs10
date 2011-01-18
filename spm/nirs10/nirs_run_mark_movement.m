function out = nirs_run_mark_movement(job)
%filename prefix 
prefix = 'm'; %for "movement"
win = job.mvt_window_length;
cutoff = job.mvt_cutoff/100;
MTh = job.sum_mvt_threshold/100;
group_consecutive_markers = 1; %boolean
try 
    job.min_session_duration;
catch
    job.min_session_duration = 60;
end

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
        winfs = floor(win*fs);
        
        for f=1:size(rDtp,1)
            d = fopen_NIR(rDtp{f,1},NC);
        
            %array to store very bad points
            pvb = zeros(size(d));
            
            %over each window, find channels that exceed the cutoff
            for i=1:size(d,2)-winfs 
                %bad channels, as a column vector
                cb = false(size(d,1),1);
                e = d(:,i:i+winfs);
                mx = max(e,[],2);
                mn = min(e,[],2);
                %md = median(e,2);
                %find channels
                cb( (mx./mn-1) > cutoff) = 1;
                %store bad points over the window
                pvb(cb,i:i+winfs) = 1;
            end
            
            %vector of bad points common to lp*NC channels; as a column
            bp = zeros(size(d,2),1);
            s = sum(pvb,1);
            bp(s>=floor(MTh*NC)) = 1;
            %indices of bad points
            bpi = find(bp);            
            bpd = ones(length(bpi),1); %all ones since duration = one data point
            
            %group consecutive markers, as an option
            if group_consecutive_markers
                [bpi bpd] = find_marker_durations(bpi);
            end
            
            %remove subsessions that are too short, less than 60s in length
            if ~isempty(bpi)
                tbpi = [];
                tbpd = [];
                %begin
                if bpi(1)-1 < fs*job.min_session_duration
                    %discard the first interval
                    next_bpi = 1;
                    next_bpd = bpi(1);
                else
                    next_bpi = bpi(1);
                    next_bpd = bpd(1);
                end
                %main loop
                for i=1:length(bpi)-1
                    if bpi(i+1)-bpi(i)-bpd(i) < fs*job.min_session_duration
                        %discard this subsession
                        next_bpd = next_bpd + bpi(i+1)-bpi(i);
                    else
                        %arrays of sufficiently spaced bad points
                        tbpi = [tbpi next_bpi];
                        tbpd = [tbpd next_bpd];
                        next_bpi = bpi(i+1);
                        next_bpd = bpd(i+1);
                    end
                end
                %conclude
                if size(d,2)-bpi(end)-bpd(end) < fs*job.min_session_duration 
                    %tbpi(end) doesn't change, but need to extend tpbd to the
                    %end
                    tbpi = [tbpi next_bpi];
                    tbpd = [tbpd size(d,2)-next_bpi-1];
                else
                    tbpi = [tbpi next_bpi];
                    tbpd = [tbpd next_bpd];
                end
                bpi = tbpi;
                bpd = tbpd;
                %catalog subsessions
                %start and end indices for subsessions
                if bpi(1) == 1
                   si = [];
                   ei = [];
                else
                    si = 1;
                    ei = bpi(1)-1;
                end

                for i=1:length(bpi)-1
                    si = [si bpi(i)+bpd(i)+1];
                    ei = [ei bpi(i+1)-1];
                end

                if bpi(end)+bpd(end)+1 == size(d,2)
                    %no more sessions
                else
                    si = [si bpi(end)+bpd(end)+1];
                    ei = [si size(d,2)];
                end
            else
                si = 1;
                ei = size(d,2);
            end
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            outfile = fullfile(dir1,[prefix fil1 ext1]);
            fwrite_NIR(outfile,d);
            %add outfile name to NIRS
            if f == 1
                NIRS.Dt.fir.pp(lst+1).pre = 'mark_movement';
                NIRS.Dt.fir.pp(lst+1).job = job;
            end
            NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
            NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
            NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
            NIRS.Dt.fir.pp(lst+1).si{f,1} = si;
            NIRS.Dt.fir.pp(lst+1).ei{f,1} = ei;
        end                 
        save(job.NIRSmat{Idx,1},'NIRS'); 
    catch
        disp(['Could not mark movement for subject' int2str(Idx)]);
    end   
end
out.NIRSmat = job.NIRSmat;