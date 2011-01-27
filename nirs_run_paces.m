function out = nirs_run_paces(job)
prefix = 'r'; %heart "rate"
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry
% 2010-10

%Parameters for the Short Term Fourier Transform :
windo = @hann;                          % windo = job.STFT_param.win_type; %Hamming not coded up?% Window Hanning
windo_width = job.STFT_param.win_width; % width : 6 seconds
n = job.STFT_param.Nprobe;              % number of probes
fft_size = job.STFT_param.fft_size;     % size of fft

remove_no_heartbeat = job.remove_no_heartbeat; %Boolean to remove channels with no heartbeat from data files
%detect_wavelength = job.detect_wavelength; %Wavelength number(s) for detection of heart beat

% parameters for Philippe's method to select channels
MinHeartRate = job.choice_method.MinHeartRate;
MaxHeartRate = job.choice_method.MaxHeartRate;
InternalMinHeartRate = job.choice_method.InternalMinHeartRate;
InternalMaxHeartRate = job.choice_method.InternalMaxHeartRate;
MaxHeartStdev = job.choice_method.MaxHeartStdev;
%to store output information for all subjects on heartrate
heartpace_all = {};

%list of NIRS.mat locations, one per subject
for Idx=1:size(job.NIRSmat,1)
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});         %Load NIRS.mat information
        heartpace = {};
        
        fs = NIRS.Cf.dev.fs;        %NIRS sampling frequency
        NC = NIRS.Cf.H.C.N;        %NIRS total number of channels
        nc = NIRS.Cf.H.C.nc;        %NIRS number of channels (counting 1 for all wavelengths)
        
        lst = length(NIRS.Dt.fir.pp);        %use last step of preprocessing
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        
        slab_width = floor(windo_width*fs)-1;      %width of a slab of signal
        win = window(windo,floor(windo_width*fs)); %win is ?
        
        for f=1:size(rDtp,1)        %loop over data files
            d = fopen_NIR(rDtp{f},NC);            %load data
            ns = size(d,2);
            const1 = 1/(n-1)*(ns-windo_width*fs);             %const1 is the window size in data points appproximately ns/n.
            
            % Heart Rate
            heart.pace = zeros(nc,ns);            %????????????CAREFUL, pace and energie have been transposed
            heart.energie = zeros(nc,ns);
            
            ci =0;
            
            k1 = []; % liste des Channels to keep
            
            %Loop over channels
            for Ci=1:NC
                %Only check channels corresponding to HBO hence 830
                wl_nm = NIRS.Cf.dev.wl;
                if wl_nm(1,NIRS.Cf.H.C.wl(Ci))== 830
                    ci = ci+1;
                    dCi = d(Ci,:);
                    
                    %Selects pairs where physiology appears (from Script_ProcessDataAccelerometer)
                    for ispectre = 1:n % n spectra to be analysed...
                        Ni = floor(1+(ispectre-1)*const1);    %Ni is the beginning of the window in data points
                        slab = dCi(Ni:(Ni+slab_width)).*win'; %slab is the data in that window multiplied by win
                        fft_slab = abs(fft(slab,fft_size));   %the absolute value of the Fourier transform of slab
                        fft_freq_step = fs/fft_size;
                        % for plotting purpose: fft_freq_scale = (1:fft_size)*fft_freq_step;
                        
                        %analyse if spectrum contains peak showing there is a regularity in
                        %the signal
                        outbeattest = nirs_criugm_trackpace(fft_slab,fft_freq_step,InternalMinHeartRate,InternalMaxHeartRate);
                        
                        heart.energie(ci,Ni:Ni+slab_width) = outbeattest.heart.energie;
                        heart.pace(ci,Ni:Ni+slab_width) = outbeattest.heart.pace;
                        
                        %         mayer_energie(Ni:Ni+slab_width,Ci) = outbeattest.mayer.energie;
                        %         mayer_pace(Ni:Ni+slab_width,Ci) = outbeattest.mayer.pace;
                    end
                    
                    %% on garde ou pas les paires
                    if job.choice_method.heart_method==0
                        v = heart.pace(ci,:);
                        median1 = median(v);
                        std1 = std(v);
                        count = 0;
                        if MinHeartRate < median1 && median1 < MaxHeartRate && ...
                                std1 < 2*MaxHeartStdev
                            %keep channel
                            %count=count+1;
                            k1 = [k1 Ci];
                            %clean up aberrant values
                            %for i4=1:length(v)
                            %if v(i4) <= MinHeartRate || v(i4) >= MaxHeartRate
                            %quick fix, better would be to take average
                            %value at ends of bad intervals
                            %    v(i4) = median1;
                            %end
                            %end
                            heart.pace(ci,:) = v;
                        else
                            heart.pace(ci,:) = zeros(1,ns);
                        end
                        
                    elseif  job.choice_method.heart_method==1
                        test_C = sum(heart.pace(ci,:)); % on a le battement ou 0 si no battementm il faut juste voir si on a assez de battement d'ou la somme
                        %%%%%%la valeur reste a fixer...automatiquement
                        if(test_C<2*10^4)%median(Cok_temp))%2.3*10^4) %Conditions for a good heart beat
                            heart.pace(ci,:) = zeros(1,ns);
                            heart.energie(ci,:) = zeros(1,ns);
                            k1 = [k1 Ci];
                        end
                    end %end choice channels               
                end %end : boucle traitement des channels a 830
            end % end channels Ci
            
            %remove only channels that were not detected in the first
            %session - to harmonize all the sessions
            %§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             if harmonize==1
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 if f == 1
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                     first_k1 = k1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             end
            
            % % % % % % % %             if remove_no_heartbeat
% % % % % % % %                 %shrink the data to keep only desired channels
% % % % % % % %                 d = d(first_k2,:);
% % % % % % % %                 heart.energie = heart.energie(first_k2,:);
% % % % % % % %                 heart.pace = heart.pace(first_k2,:);
% % % % % % % %             end
            
            %output a heart rate per minute
            heart.pace = 60* heart.pace;
            %last entry is zero - replace by previous entry
           %????? heart.pace(:,end) = heart.pace(:,end-1);
            
            %needs to be generalized to more than one session
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            outheartfile = fullfile(dir1, [fil1 '_heart_pace.mat']);
            save(outheartfile,'heart');
            
            % on calcule les decours temporels de ces bonnes paires et on enleve les
            % artefacts de mouvements...(voir si on peut pas trier a plus
            % haur niveau sur le rapport de l'energie du battement par rapport
            % a l'energie totale dans le signal)
            
            %median1 = median(heart.pace);
            %std1 = std(heart.pace);
            
%             [dir1,fil1,ext1] = fileparts(rDtp{f});
            outfile = fullfile(dir1,[prefix fil1 ext1]);
            fwrite_NIR(outfile,d);
            %add outfile name to NIRS
            if f == 1
                NIRS.Dt.fir.pp(lst+1).pre = 'heart_rate';
                NIRS.Dt.fir.pp(lst+1).job = job;
            end
            NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
            bpi = [];
            bpd = [];
            NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
            NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
            NIRS.Dt.fir.ht{f,1} = outheartfile;
            %Take median heart rate on the better channels
            if remove_no_heartbeat
                cR = median(heart.pace(size(heart.pace,1)/2,:),1)'; %cR, fR are column vectors
            else
                cR = median(heart.pace(first_k1,:),1)';
            end
            NIRS.Dt.fir.Sess(f).cR = {cR};
            try
                fR = NIRS.Dt.fir.Sess(f).fR{1};
                %correlation of heart rates between ECG rate and NIRS rate
                NIRS.Dt.fir.Sess(f).cCor = corr(fR,cR);
                %correlation of changes in heart rates
                NIRS.Dt.fir.Sess(f).cdCor = corr(diff(fR),diff(cR));
            catch
            end
            heartpace  = [heartpace; outheartfile];
        end
%         if remove_no_heartbeat
%             %update the NIRS structure
%             NIRS.Cf.H.C.N = length(first_k2);
%             try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(first_k2); end
%             try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,first_k2); end
%             try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(first_k2); end
%             try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(first_k2); end
             try NIRS.Cf.H.C.ok = NIRS.Cf.H.C.ok(k1); end
%         end
        save(job.NIRSmat{Idx,1},'NIRS');
        heartpace_all = [heartpace_all; heartpace];
    catch
        disp(['Could not analyze heart rate for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;
out.heartpace = heartpace_all;