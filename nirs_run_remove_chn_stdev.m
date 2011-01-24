function out = nirs_run_remove_chn_stdev(job)
%filename prefix 
prefix = 'd'; %for standard "deviation"
DelPreviousData  = job.DelPreviousData;
try 
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

threshold_stdev = job.threshold_stdev;
win = job.window_stdev;

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
        wl = NIRS.Cf.dev.wl;
        
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
                stdev(:,i1) = std(d(:,i1:i1+floor(fs * win)),0,2);
                med2(:,i1) = median(d(:,i1:i1+floor(fs * win)),2);
            end %end for i1
            %channelwise median of standard deviations
            med1 = median(stdev./med2,2);
            %identify good channels gc
            nc = NC/length(wl);
            gc = ones(nc,1);
            k1 = [];
            for i2=1:nc
                for i3=1:length(wl)
                    if med1(i2+(i3-1)*nc)> threshold_stdev
                        gc(i2) = 0; %bad channel
                    end
                end
                if gc(i2)
                    k1 = [k1 i2]; %indices for first wavelength
                end
            end
            k2 = k1;
            %channel indices for all wavelengths
            for i3=1:length(wl)-1
                k2 = [k2 k2+nc];
            end
            %Note that 
            %find(1-gc)
            %gives the channels removed
            %and
            %length(find(1-gc))
            %is the number of channels that were removed
            d = d(k2,:);
            %update NIRS matrix
            NIRS.Cf.H.C.N = length(k2);
            try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(k2); end
            try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,k2); end
            try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(k2); end 
            try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(k2); end
            try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(k2); end 
                        
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
                NIRS.Dt.fir.pp(lst+1).pre = 'remove_chn_stdev';
                NIRS.Dt.fir.pp(lst+1).job = job;
            end
            NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
            NIRS.Dt.fir.pp(lst+1).kept{f,1} = k2; %kept channels
        end 
        %keep copy of original NIRS structure
        [dir1 fil1 ext1] =fileparts(job.NIRSmat{Idx,1});
        %the old file will be dNIRS.mat
        %out1 = fullfile(dir1,[prefix fil1 ext1]);
        if NewDirCopyNIRS
            save(fullfile(dir2,'NIRS.mat'),'NIRS');            
        else
            save(job.NIRSmat{Idx,1},'NIRS'); 
        end

        save(job.NIRSmat{Idx,1},'NIRS');    
    catch
        disp(['Could not remove channels based on standard deviation for subject' int2str(Idx)]);
    end
    
end
out.NIRSmat = job.NIRSmat;