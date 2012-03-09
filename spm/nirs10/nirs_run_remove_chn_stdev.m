function out = nirs_run_remove_chn_stdev(job)
%filename prefix
prefix = 'd'; %for standard "deviation"

threshold_stdev = job.threshold_stdev;
win = job.window_stdev;

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'remov_std_OK') || job.force_redo)
            %use last step of preprocessing
            lst = length(NIRS.Dt.fir.pp);
            rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
            NC = NIRS.Cf.H.C.N;
            fs = NIRS.Cf.dev.fs;
            wl = NIRS.Cf.dev.wl;
            nc = NC/length(wl);
            gc = ones(nc,size(rDtp,1));
            npt = floor(fs * win);
            %first loop over all sessions to find all bad channels
            for f=1:size(rDtp,1)
                d = fopen_NIR(rDtp{f,1},NC);
                ns = size(d,2);
                %number of complete windows
                nwin = floor(ns/(fs * win));
                stdev = zeros(NC,nwin);
                med2  = zeros(NC,nwin);
                for i1=1:nwin
                    %standard deviation channelwise over
                    %non overlapping windows
                    if i1==nwin
                        lpt = min((i1*npt),size(d,2));
                    else
                        lpt = i1*npt;
                    end
                    %stdev(:,i1) = std(d(:,((i1-1)*npt+1):lpt),0,2);
                    d1 = d(:,((i1-1)*npt+1):lpt);
                    md2 =  mean(d1,2)*ones(1,size(d1,2));
                    stdev(:,i1) = std(d1./md2,0,2);
                    %med2(:,i1) = abs(median(d(:,i1:i1+floor(fs * win)),2));
                    %mean2(:,i1) = abs(mean(d(:,i1:i1+floor(fs * win)),2));
                end %end for i1
                
                %channelwise median of standard deviations
                %med1 = median(stdev./med2,2);
                med1 = median(stdev,2);
                %mean1 = mean(stdev./mean2,2);
                %identify good channels gc
                
                %k1 = [];
                for i2=1:nc
                    for i3=1:length(wl)
                        if med1(i2+(i3-1)*nc)> threshold_stdev
                            gc(i2,f) = 0; %bad channel
                        end
                    end
                    %if gc(i2)
                    %    k1 = [k1 i2]; %indices for first wavelength
                    %end
                end
                %k2 = k1;
                %channel indices for all wavelengths
                %for i3=1:length(wl)-1
                %    k2 = [k2 k2+nc];
                %end
            end
            k1 = [];
            gc = sum(gc,2);
            for i4=1:nc
                if gc(i4) == size(rDtp,1)
                    %this is a good channel
                    k1 = [k1 i4];
                end
            end
            k2 = k1;
            %channel indices for all wavelengths
            for i3=1:length(wl)-1
                k2 = [k2 k1+i3*nc];
            end
            %Now write the data session by session
            for f=1:size(rDtp,1)
                d = fopen_NIR(rDtp{f,1},NC);
                %ns = size(d,2);
                %Note that
                %find(1-gc)
                %gives the channels removed
                %and
                %length(find(1-gc))
                %is the number of channels that were removed
                d = d(k2,:);
                
                [dir1,fil1,ext1] = fileparts(rDtp{f});
                
                    outfile = fullfile(dir1,[prefix fil1 ext1]);
         
                fwrite_NIR(outfile,d);
                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'remove_chn_stdev';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
                NIRS.Dt.fir.pp(lst+1).kept{f,1} = k2; %kept channels
            end
            %update NIRS matrix
            NIRS.Cf.H.C.N = length(k2);
            try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(k2); end
            try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,k2); end
            try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(k2); end
            try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(k2); end
            try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(k2); end
            
            %keep copy of original NIRS structure
            [dir1 fil1 ext1] =fileparts(job.NIRSmat{Idx,1});
            %the old file will be dNIRS.mat
            %out1 = fullfile(dir1,[prefix fil1 ext1]);
            NIRS.flags.remov_std_OK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
        
        save(job.NIRSmat{Idx,1},'NIRS');
    catch
        disp(['Could not remove channels based on standard deviation for subject' int2str(Idx)]);
    end
    
end
out.NIRSmat = job.NIRSmat;