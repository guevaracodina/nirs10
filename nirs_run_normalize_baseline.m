function out = nirs_run_normalize_baseline(job)
%filename prefix 
prefix = 'b'; %for "baseline"
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
            if markers_available
                %temporary for data
                td = zeros(size(d));
                %group consecutive markers, as an option
                if group_consecutive_markers
                    [bpi bpd] = find_marker_durations(bpi);
                end
                for i=1:length(si)
                    try
                        e = d(:,si(i):ei(i));
%                         switch i
%                             case 1
%                                 if bpi(1) == 1
%                                     e = [];
%                                 else
%                                     %interval before first marker
%                                     e=d(:,1:bpi(1)-1);
%                                 end
%                             case length(bpi)+1
%                                 if size(d,2) == bpi(end)+bpd(end)+1
%                                     e = [];
%                                 else
%                                     %interval after last marker
%                                     e=d(:,bpi(end)+bpd(end)+1:end);
%                                 end
%                             otherwise
%                                 %all other intervals
%                                 e=d(:,bpi(i-1)+bpd(i-1)+1:bpi(i)-1);
%                         end
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
                                td(:,si(i):ei(i)) = (100+(e-div_factor))*job.Analyzer_sf;
                            else
                                td(:,si(i):ei(i)) = e./div_factor*job.Analyzer_sf; 
                            end
                            
%                             switch i
%                                 case 1
%                                     %interval before first marker
%                                     if add_or_mult
%                                         td(:,1:bpi(1)-1) = (100+(e-div_factor))*job.Analyzer_sf;
%                                     else
%                                         td(:,1:bpi(1)-1) = e./div_factor*job.Analyzer_sf;
%                                     end
%                                 case length(bpi)+1
%                                     if add_or_mult
%                                         %interval after last marker
%                                         td(:,bpi(end)+bpd(end)+1:end) = (100+(e-div_factor))*job.Analyzer_sf;
%                                     else
%                                         td(:,bpi(end)+bpd(end)+1:end) = e./div_factor*job.Analyzer_sf;
%                                     end
%                                 otherwise
%                                     if add_or_mult
%                                     %all other intervals
%                                         td(:,bpi(i-1)+bpd(i-1)+1:bpi(i)-1) = (100+(e-div_factor))*job.Analyzer_sf;
%                                     else
%                                         td(:,bpi(i-1)+bpd(i-1)+1:bpi(i)-1) = e./div_factor*job.Analyzer_sf;
%                                     end
%                             end
                        end
                    catch
                        %do nothing
                    end
                    %interpolate linearly over bad intervals??? or do
                    %nothing??
                    
                    
                end
                %replace d - this sets intervals with movement to zero
                d=td;
            else
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
                    %normalize median to 100 uM
                    
                    d = (100+(d-div_factor))*job.Analyzer_sf; 
                else
                    d = d./div_factor*job.Analyzer_sf; 
                end
            end
           
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            outfile = fullfile(dir1,[prefix fil1 ext1]);
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
        save(job.NIRSmat{Idx,1},'NIRS');   
    catch
        disp(['Could not normalize to baseline for subject' int2str(Idx)]);
    end       
end
out.NIRSmat = job.NIRSmat;