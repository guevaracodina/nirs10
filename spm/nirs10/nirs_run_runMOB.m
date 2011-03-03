function out = nirs_run_runMOB(job)

% Script_ProcessDataAccelerometer
% Determines variability during walk
% Clément Bonnéry - clement.bonnery@polymtl.ca
% Written during 04/10

% First variable is vertical axis
% Second variable is horizontal axis (direction perpendiculary to the axis
% of walking)
% Third variable is horizontal axis (direction of walking)

%Big loop over all subjects
sN = size(job.subj,2);
for is=1:sN
    
    %%  From Accelerometer to Matlab
    
    
    % (uncomment(ctrl+T) following lines (from 20 to 24))to import datas to Matlab :
    % - double clic on the datas (.csv) you want to import (import wizard) in
    % current directory window (next then finish : the datas are stored in the workspace)
    % - enter the right path and the right name then launch the script :
    filemat = csv2cell(job.subj(1,is).acc_file{:});
    data = cell2mat(filemat(11:end,:));
    textdata = filemat(1:10,1);
    Te = 33.3*10^(-3);
    [dir, name, ~] = fileparts(job.subj(1,is).acc_file{:});
    % acc_file = job.acc_file{:};
    % load(acc_file,'csv');
    save([dir '\' name '.mat'], '-mat', 'data', 'Te', 'textdata');
    
    
    Fe = 1/Te;
    t =1:length(data(:,1)); % time scale
    t = Te*t;
    
    %% Extract period while the subject is walking
    % criteron : each step correponds to a period on the low passed signal
    %-- positive acceleration (superior to mean) correponds to the falling of
    %the foot [at each maximum, vibrations are seen]
    %-- Low Pass Filter to filter all the noise : callibrated around 2 steps a
    %second
    
    datV = data(:,2);
    
    % Low Pass Filter
    Wn = 2*2*Te;                % Cutting frequency : 2Hz
    [b,a] = butter(4,Wn,'low');
    datV_F = filtfilt(b,a,datV);
    
    % figure;plot(datV_F);title('Vertical axis datas low pass filtered');
    
    %% Dicrimine periode while subjects are walking (sinusoidal curve) or
    % resting (noise filtered to a certain extend thanks to Low Pass Filter)
    
    % Short Term Fourier Transform :
    windo = @hann;                          % Window Hanning
    windo_width = 6;                        % width : 6 seconds
    n = 4000;                               % number of probes
    fft_size = 1024;                        % size of fft
    slab_width = floor(windo_width*Fe)-1;   % width of a slab of signal
    
    walkBool = zeros(floor(slab_width/2)+length(datV),1);
    win = window(windo,floor(windo_width*Fe));
    
    offset = 70;                            % time offset (to make walk and signal match)
    count = 1;                              % to remember moments when bool changes
    bool_mem = 0;
    
    for i = 1:n % n spectra to be analysed...
        Ni = floor(1+(i-1)/(n-1)*(length(datV)-windo_width*Fe));
        slab = datV(Ni:(Ni+slab_width)).*win;
        fft_slab = abs(fft(slab,fft_size));
        fft_freq_step = Fe/fft_size;
        fft_freq_scale = (1:fft_size)*fft_freq_step;
        
        %analyse if spectrum contains pick showing there is a regularity in
        %the signal
        [bool, E] = walktest(fft_slab,fft_size);
        if bool_mem~=bool
            time2remember(count) = Ni+offset;
            count = count+1;
        end
        
        E_evo(i) = E;
        walkBool(Ni+offset:Ni+offset+slab_width,1) = bool;
        bool_mem = bool;
    end
    
    walk = walkBool;
    
    % % display in red walking period of time
    figure;
    subplot(2,1,1)
    hold on
    plot(datV_F);
    plot(datV,'g');
    plot(mean(datV)+300*walk(1:length(datV)),'r');
    hold off
    
    subplot(2,1,2)
    plot(E_evo)
    
    
    %% Determine period while variability can be accurately measured
    
    % Security margins of 2s inside the detected time periods
    time_margin = 2;%2s
    %walk_stationary = zeros(length(datV_F));
    
    for i_t2r = 0:floor(length(time2remember)/2)-1
        N_min = time2remember(2*i_t2r+1)+time_margin/Te;
        N_max = time2remember(2*(i_t2r+1))-time_margin/Te;
        if N_min < N_max % Too small walks are not consistent
            walk_stationary(N_min:N_max) = 1;%datV_F(N_min:N_max);
        end
    end
    
    
    %% Measure variability and give results
    
    % Maxima detection
    compteurBloc= 0;
    nBlocs = 0;
    begin_step = zeros(0);
    
    for i=2:length(walk_stationary)-1
        if walk_stationary(i)==1
            if walk_stationary(i-1)==0
                nBlocs = nBlocs+1;
            end
        end
    end
    
    compteur = zeros(nBlocs,1);
    
    for i=2:length(walk_stationary)-1
        if walk_stationary(i)==1
            if walk_stationary(i-1)==0
                compteurBloc = compteurBloc+1;
            end
            if (datV_F(i)>datV_F(i+1) && datV_F(i)>datV_F(i-1))% then step
                compteur(compteurBloc,1) = compteur(compteurBloc,1)+1;
                begin_step(compteurBloc,compteur(compteurBloc)) = i;%remember beginning of step
            end
        end
    end
    
    % figure;
    % subplot(2,1,1)
    % plot(walk);
    % hold on
    % plot(walk_stationary,'r');
    % hold off
    % subplot(2,1,2)
    % %plot(begin_step)
    
    % Traitement des temps entre les pas
    length_step = zeros(0);
    moy_length_step = zeros(size(begin_step,1),1);
    std_length_step = zeros(size(begin_step,1),1);
    
    for i=1:size(begin_step,1)
        for j=2:size(begin_step(i,:),2)
            length_step(i,j-1) = (begin_step(i,j) - begin_step(i,j-1))*Te;
            if length_step(i,j-1)<0
                length_step(i,j-1) =0;
            end
        end
        moy_length_step(i) = mean(length_step(i,:));
        std_length_step(i) =  std(length_step(i,:));
    end
    
    %% Visualisation des résultats
    indice =1;
    visu = zeros(0);
    
    for i =1:length(compteur)
        if compteur(i,1)>20
            visu(indice,1) = compteur(i,1);
            visu(indice,2) = moy_length_step(i,1);
            visu(indice,3) = std_length_step(i,1);
            visu(indice,4) = begin_step(i,1)*Te;
            indice = indice+1;
        end
    end
    
    visu = visu';
    
    figure;
    plot(datV(1:length(walk_stationary)));
    hold on
    plot(mean(datV)+300*walk_stationary,'r');
    hold off
    
    save([dir '\visu.mat'],'visu');
end
out = 1;
end