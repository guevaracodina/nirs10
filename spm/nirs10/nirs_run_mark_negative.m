function out = nirs_run_mark_negative(job)
%To do: also mark saturated channels
%filename prefix 
prefix = 'n'; %for "negative"
DelPreviousData  = job.DelPreviousData;
try 
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

STh = job.sum_neg_threshold/100;
group_consecutive_markers = 1; %boolean

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
           
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        
        for f=1:size(rDtp,1)
            d = fopen_NIR(rDtp{f,1},NC);
        
            %Some values of d are actually negative, even though this is an
            %optical intensity
            %Find negative values
            
            %channel indices: ch; which data points are bad: p
            [ch p] = find(d<=0);
            
            %array to store very bad points
            pvb = zeros(size(d));
            
            %interpolate to remove isolated bad points
            %loop over channels
            for j=1:NC
                %indices of bad points for channel j
                pb = p(ch==j);
                tpvb = [];
                %first bad point
                if pb, d(j,pb(1)) = d(j,pb(1)-1); end
                %loop over bad points, skipping first and last
                for i=2:length(pb)-1
                    if pb(i)-pb(i-1)>1 && pb(i+1)-pb(i)>1
                        %isolated bad point, interpolate
                        d(j,pb(i)) = (d(j,pb(i)-1)+d(j,pb(i)+1))/2;
                    else
                        %take value to the left, but mark as very bad point
                        d(j,pb(i)) = d(j,pb(i)-1);
                        tpvb = [tpvb pb(i)];
                    end
                end
                %last bad point
                if length(pb)>1, d(j,pb(end)) = d(j,pb(end)-1); end
                %store very bad points
                pvb(j,tpvb) = 1;
            end
            %vector of bad points common to lp*NC channels; as a column
            bp = zeros(size(d,2),1);
            s = sum(pvb,1);
            bp(s>=floor(STh*NC)) = 1;
            %indices of bad points
            bpi = find(bp);
            bpd = ones(length(bpi),1); 
            
            %group consecutive markers, as an option
            if group_consecutive_markers
                [bpi bpd] = find_marker_durations(bpi);
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
                NIRS.Dt.fir.pp(lst+1).pre = 'mark_negative';
                NIRS.Dt.fir.pp(lst+1).job = job;               
            end
            NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
            NIRS.Dt.fir.pp(lst+1).kept{f,1} = 1:NC;
            NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
            NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
        end  
        if NewDirCopyNIRS
            save(fullfile(dir2,'NIRS.mat'),'NIRS');            
        else
            save(job.NIRSmat{Idx,1},'NIRS'); 
        end
    catch
        disp(['Could not mark negative values for subject' int2str(Idx)]);
    end        
end
out.NIRSmat = job.NIRSmat;
