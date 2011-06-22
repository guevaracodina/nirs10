function out = nirs_run_criugm_paces(job)
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry
% 2010-10

prefix = 'r'; %heart "rate"
DelPreviousData  = job.DelPreviousData;
try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

%Boolean to remove channels with no heartbeat from data files
remove_no_heartbeat = job.remove_no_heartbeat;

%Parameters for the FFT
% Short Term Fourier Transform :
try
    %resting state
    % width : 6 seconds
    switch job.heart_rate_cfg.heart_resting.STFT_param.win_type
        case 0
            windo = @hann;
        case 1
            windo = @hamm;
        case 2
            windo = @rect;
    end
    windo_width = job.heart_rate_cfg.heart_resting.STFT_param.win_width;
    % number of probes
    n = job.heart_rate_cfg.heart_resting.STFT_param.Nprobe;
    fft_size = job.heart_rate_cfg.heart_resting.STFT_param.fft_size;
    MinHeartRate = job.heart_rate_cfg.heart_resting.MinHeartRate/60;
    MaxHeartRate = job.heart_rate_cfg.heart_resting.MaxHeartRate/60;
    InternalMinHeartRate = job.heart_rate_cfg.heart_resting.InternalMinHeartRate/60;
    InternalMaxHeartRate = job.heart_rate_cfg.heart_resting.InternalMaxHeartRate/60;
    MaxHeartStdev = job.heart_rate_cfg.heart_resting.MaxHeartStdev/60;
    ex =0;
catch
    %aerobic exercise
    % width : 6 seconds
    switch job.heart_rate_cfg.heart_exercise.STFT_param.win_type
        case 0
            windo = @hann;
        case 1
            windo = @hamm;
        case 2
            windo = @rect;
    end
    windo_width = job.heart_rate_cfg.heart_exercise.STFT_param.win_width;
    n = job.heart_rate_cfg.heart_exercise.STFT_param.Nprobe;
    fft_size = job.heart_rate_cfg.heart_exercise.STFT_param.fft_size;
    InternalMinHeartRate = 0.5;
    InternalMaxHeartRate = 4; %or 3?
    ex =1;
end

%to store output information for all subjects on heartrate
heartpace_all = {};

for Idx=1:size(job.NIRSmat,1)
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        % % % % % % %         heartpace = {};
        fs = NIRS.Cf.dev.fs;
        NC = NIRS.Cf.H.C.N;
        wl = NIRS.Cf.dev.wl;
        
        slab_width = floor(windo_width*fs)-1;   % width of a slab of signal
        win = window(windo,floor(windo_width*fs)); % window to attenuate high frequencies when working on time slabs
        
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
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        %loop over data files
        for f=1:size(rDtp,1)
            if strcmp(NIRS.Cf.dev.n,'CW6')
                d = fopen_NIR(rDtp{f},NC);
            elseif strcmp(NIRS.Cf.dev.n,'CW5')
                fd =load(rDtp{f});
                d =fd.d';
            end
            ns = size(d,2);
            % Heart Rate
            % Frequency of Mayer Waves not calculated: to be done
            
            %CAREFUL, pace and nrgy have been transposed
            heart.from = rDtp{f};
            heart.pace = zeros(NC,ns);
            heart.nrgy = zeros(NC,ns);
            
            const1 = 1/(n-1)*(ns-windo_width*fs);%const1 is the window size in data points : appproximately ns/n.
            
            % on repere la qualite de chaque paire en regardant si on a les battements
            % physiologiques :
            %Channels to keep:
            k1 = [];
            for Ci=1:NC            %Loop over channels
                %Only check channels corresponding to selected wavelength(s)
                if any(NIRS.Cf.H.C.wl(Ci)== HbO_like)
                    dCi = d(Ci,:);
                    
                    %Selects pairs where physiology appears (from Script_ProcessDataAccelerometer)
                    for ispectre = 1:n % n spectra to be analysed...
                        %Ni is the beginning of the window in data points
                        Ni = floor(1+(ispectre-1)*const1);
                        %slab is the data in that window multiplied by win
                        slab = dCi(Ni:(Ni+slab_width)).*win';
                        %the absolute value of the Fourier transform of
                        %slab
                        fft_slab = abs(fft(slab,fft_size));
                        fft_freq_step = fs/fft_size;
                        % for plotting purpose: fft_freq_scale = (1:fft_size)*fft_freq_step;figure;plot(fft_freq_scale,fft_slab)
                        
                        %analyse if spectrum contains peak showing there is a regularity in
                        %the signal
                        outbeattest = nirs_criugm_trackpace(fft_slab,fft_freq_step,InternalMinHeartRate,InternalMaxHeartRate);
                        
                        heart.nrgy(Ci,Ni:Ni+slab_width) = outbeattest.heart.nrgy;
                        heart.pace(Ci,Ni:Ni+slab_width) = outbeattest.heart.pace;
                    end
                    
                    if ex==0 % resting state
                        k1 = [k1 Ci];
                        %                         v = heart.pace(Ci,:);
                        %                         median1 = median(v);
                        %                         std1 = std(v);
                        %                         count = 0;
                        %                         if MinHeartRate < median1 && median1 < MaxHeartRate && ...
                        %                                 std1 < 2*MaxHeartStdev
                        %                             %keep channel
                        %                             %count=count+1;
                        %                             %channel list that we keep:
                        %                             k1 = [k1 Ci];
                        %                             %clean up aberrant values
                        %                             %for i4=1:length(v)
                        %                             %if v(i4) <= MinHeartRate || v(i4) >= MaxHeartRate
                        %                             %quick fix, better would be to take average
                        %                             %value at ends of bad intervals
                        %                             %    v(i4) = median1;
                        %                             %end
                        %
                        %                             %end
                        %                             heart.pace(Ci,:) = v;
                        %                         else
                        %                             heart.pace(Ci,:) = zeros(1,ns);
                        %                         end
                    else % during exercise
                        test_C = sum(heart.pace(Ci,:)); % on a le battement ou 0 si no battementm il faut juste voir si on a assez de battement d'ou la somme
                        %%%%%%la valeur reste a fixer...automatiquement
                        if(test_C<2*10^4)%median(Cok_temp))%2.3*10^4) %Conditions for a good heart beat
                            heart.pace(Ci,:) = zeros(1,ns);
                            heart.nrgy(Ci,:) = zeros(1,ns);
                        else
                            k1 = [k1 Ci];
                        end
                    end
                end %end if any
            end
            
            %%% a cause du fait qu'on recupere le numero de la ligne et
            %%% qu'ensuite on complete, on inverse l'ordre des longueurs
            %%% d'ondes... DANGEUREUX, non ???
            k2 = k1;
            %only valid if detection was done on first wavelength only
            %if detect_wavelength == 1
            %complete to all wavelengths
            wl = NIRS.Cf.dev.wl;
            nc = NC/length(wl);
            %channel indices for all wavelengths
            for i3=1:length(wl)-1
                %include HbR_like channels
                for i4=1:length(HbR_like)
                    k2 = [k2 k2+nc*(-HbO_like+HbR_like(i4))];
                end
            end
            %end
            %remove only channels that were not detected in the first
            %session - to harmonize all the sessions
            if f == 1
                first_k1 = k1;
                first_k2 = k2;
            end
            if remove_no_heartbeat
                %shrink the data to keep only desired channels
                d = d(first_k2,:);
                heart.nrgy = heart.nrgy(first_k2,:);
                heartfRpace = heart.pace(k2,:);
                heart.pace = heart.pace(first_k2,:);
            end
            
            %output a heart rate per minute
            heart.pace = 60* heart.pace;
            heartfRpace = 60*heartfRpace;
            %last entry is zero - replace by previous entry
            % cb: pk  ?? moi j'ai pas ce probleme. Peut etre a regler au
            % niveau anterieur...
            heart.pace(:,end) = heart.pace(:,end-1);
            heartfRpace(:,end) = heartfRpace(:,end-1);
            
            % % % % % %             est ce qu'on en a vraiment besoin puisqu'on
            % conserve tout dans heart ??????
            %needs to be generalized to more than one session
            %                         outheartfile = fullfile(NIRS.Dt.s.p,'heart_pace.mat');
            %                         save(outheartfile,'heart');
            
            % on calcule les decours temporels de ces bonnes paires et on enleve les
            % artefacts de mouvements...(voir si on peut pas trier a plus
            % haur niveau sur le rapport de l'nrgy du battement par rapport
            % a l'nrgy totale dans le signal)
            
            %median1 = median(heart.pace);
            %std1 = std(heart.pace);
            
            
            %Conditions for a good heart beat
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            if NewDirCopyNIRS
                dir2 = [dir1 filesep NewNIRSdir];
                if ~exist(dir2,'dir'), mkdir(dir2); end;
                outfile = fullfile(dir2,[prefix fil1 ext1]);
            else
                outfile = fullfile(dir1,[prefix fil1 ext1]);
            end
            if DelPreviousData % est ce qu;on devrait pas le supprimer dans la matrice aussi ????
                delete(rDtp{f,1});
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % tant qu'on y met pas autre chose que d, a mon avis, pas grand
            % interet
            fwrite_NIR(outfile,d);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            NIRS.Dt.fir.ht{f,1} = heart;          %heart pace and energy
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cb:utile pour la methode dans le cas de tests sans
            % exercices...
            % cb: est ce que c'est definitif ? si oui, il faut le
            % rentrer dans la matrice NIRS
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
            
            %%%%%%% ici, je rebidouille pour pouvoir revenir dans le bon
            %%%%%%% ordre....
            
            %%% ATTENTION : dna la partie regresseur, on ne tient plus
            %%% compte de savoir si les paires sont les memes pour toutes
            %%% les sessions car seul compte d'obtenir un bon regresseur
            heart_regressor=1;
            if heart_regressor
                % Attention deux pb
                %-- on a parfois une mauvaise detection des rythmes cardiaques et
                %ceux-ci sont consideres comme bons...
                %-- il faut raisonner en temporel pour voir quels rythmes se
                %correspondent d'une paire a l'autre et ainsi moyenner que ce qui
                %est bon
                %>> on cherche une paire qui est fiable (critere : elle n'a pas de grands sauts...)
                %--- on cherche le fichier du rythme qui correspond au fichier de
                %donnees
                try
                    hp = heartfRpace;
                    if sum(sum(hp))==0
                        % arbitraire
                        interpx  =[1 size(hp,2)];
                        interpY  = [60 65];
                        interpxi = (1:size(hp,2));
                        reg(1:size(hp,2)) = interp1(interpx,interpY,interpxi,'linear');
                    else
                        % on doit travailler sur les paires de HbO
                        C_HbO =[];
                        for Ci=1:NC            %Loop over channels
                            %Only check channels corresponding to selected wavelength(s)
                            if any(NIRS.Cf.H.C.wl(Ci)== HbO_like)
                                C_HbO = [C_HbO Ci];
                            end
                        end
                        
                        %%%%%
                        whp = hp/max(max(hp));
                        whp_b = zeros(size(whp));
                        
                        %%% attention : a partir de cette ligne, on
                        %%% travaille avec l'ordre des lignes correspondant
                        %%% a celui des canaux mais pas avant.....
                        % Affichage du resultat
                        if strcmp(NIRS.Cf.dev.n,'CW6')
                            whpR = zeros(NC,size(whp,2)); whpR([k1 k1-(NC/2)],:) = whp; figure;imagesc(whpR);title(['Heart pace: ' rDtp{f}]);
                        elseif strcmp(NIRS.Cf.dev.n,'CW5')
                            whpR = zeros(NC,size(whp,2)); whpR([k1 k1+(NC/2)],:) = whp; %figure;imagesc(whpR);title(['Heart pace: ' rDtp{f}]);
                        end
                        
                        level = graythresh(whpR(C_HbO,:));% Otsu
                        whp_b(C_HbO,:) = im2bw(whpR(C_HbO,:),level);
                        
                        % reconstruction
                        bouchetrou = whpR.*whp_b;
                        testq = sum(whp_b,2);
                        [val,canal] =max(testq);
                        reg = bouchetrou(canal,:);
                        
                        if val<length(reg)
                            holes =[];
                            for t=1:length(reg)
                                if reg(t)==0
                                    %cherche sur une autre paire
                                    try
                                        i=1;
                                        while bouchetrou(i,t)==0
                                            i = i+1;
                                        end
                                        reg(1,t) = bouchetrou(i,t);
                                    catch
                                        %pas de valeur disponible, on interpolera
                                        holes = [holes t];
                                    end
                                end
                            end
                            if ~isempty(holes)
                                debuts=holes(1);
                                fins=[];
                                for ih=2:length(holes)-1;
                                    if holes(ih) ~= holes(ih-1)+1
                                        debuts = [debuts holes(ih)];
                                        fins = [fins holes(ih-1)];
                                    end
                                end
                                fins=[fins holes(end)];
                                
                                if length(debuts)==1
                                    if debuts(1)==1
                                    else
                                    end
                                else
                                    % moment de l'interpolation sur les holes
                                    if(debuts(1)==1)
                                        reg(1:fins(1))=reg(fins(1)+1);
                                    end
                                    for i=1:length(debuts)
                                        if fins(i)==size(reg,2)
                                            reg(debuts(i):end)=reg(debuts(i)-1);
                                        else
                                            interpx  =[debuts(i)-1 fins(i)+1];
                                            interpY  = [reg(debuts(i)-1) reg(fins(i)+1)];
                                            interpxi = (debuts(i):fins(i));
                                            reg(debuts(i):fins(i)) = interp1(interpx,interpY,interpxi,'linear');
                                        end
                                    end
                                end
                            end
                        end
                    end
                catch
                end
                if sum(abs(diff(reg)))==0 %cst alors il faut eviter //////
                    interpx  =[1 size(hp,2)];
                    interpY  = [60 65];
                    interpxi = (1:size(hp,2));
                    reg(1:size(hp,2)) = interp1(interpx,interpY,interpxi,'linear');
                end
                NIRS.Dt.fir.Sess(f).fR{1} = reg';
                save(job.NIRSmat{Idx,1},'NIRS');
            end
        end
        
        if remove_no_heartbeat
            %update the NIRS structure
            NIRS.Cf.H.C.N = length(first_k2);
            try NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(first_k2); end
            try NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,first_k2); end
            try NIRS.Cf.H.C.wl = NIRS.Cf.H.C.wl(first_k2); end
            try NIRS.Cf.H.C.gp = NIRS.Cf.H.C.gp(first_k2); end
        end
        
        try NIRS.Cf.H.C.ok = first_k2; end %gives good channels whether they
        
        %were removed or not
        if NewDirCopyNIRS
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{Idx,1} = newNIRSlocation;
        else
            save(job.NIRSmat{Idx,1},'NIRS');
        end
        %heartpace_all = [heartpace_all; heartpace]; %cb: vide, est ce qu'on supprime ???
    catch
        disp(['Could not analyze heart rate for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;
out.heartpace = heartpace_all;