function histogram_for_simulated_spikes
%select a job 
n = 100; %number of protocols
plain_rnd = 0; %Boolean to generate a plain random protocol, 
%otherwise load a real protocol from job file
dt = [];
if plain_rnd
    ns = 400; %number of spikes
    for ROCiter=1:n
        ons = sort(unifrnd(0,900,1,ns));
        d = diff(ons);
        dt = [dt d];
    end
else
    t = spm_select([1 1],'mat');

    dir_dataSPM = 'dataSPM';
    clear LoopJob
    %Load/read job
    LoopJob = load(t(:));
    NIRSmat = LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets.addTestStimuli.NIRSmat;
    [dir1,dummy,dummy2] = fileparts(NIRSmat{1});
    dir1 = [dir1 filesep dir_dataSPM];
    %dir1 = ['W:\epiNIRSj\epi127SD' filesep dir_dataSPM]; 
    for ROCiter=1:n
        %build folder name
        testName = LoopJob.matlabbatch{1}.spm.tools.nirs10.readOnsets. ...
                addTestStimuli.testStimulusName;
        testFullName = [testName int2str(ROCiter)];
        dir_stat = LoopJob.matlabbatch{2}.spm.tools.nirs10.model_specify. ...
            wls_bglm_specify.dir1;
        dir_spm = [dir1 filesep testFullName filesep dir_stat];
        NIRS = [];
        load(fullfile(dir_spm,'NIRS.mat'));
        ons = NIRS.Dt.fir.Sess.U.ons;


        d = diff(ons);
        dt = [dt d];
        %figure;
        %hist(d,80);
    end
end
    
    
figure;
hist(dt,80)
xs = length(find(dt<1))/length(dt); %45%
xm = length(find(dt<5))/length(dt)-xs; %32%
s = std(dt); %10s
md = median(dt); %1.2s
m = mean(dt); %5.1s
xl = length(find(dt>10))/length(dt); % 
%find number of frequent and infrequent spikes per bunch
fr = []; %array of number of frequent spikes per bunch
nfr = 0; %number of frequent spikes in this bunch
for i=1:length(dt)
    if dt(i) < 3 %cutoff for frequent spikes
        nfr = nfr + 1;
    else
        if nfr
            fr = [fr nfr];
        end
        nfr = 0;
    end
end
mfr = mean(fr);
mdfr = median(fr);
sfr = std(fr);

ifr = []; %array of number of infrequent spikes per bunch
nfr = 0; %number of infrequent spikes in this bunch
for i=1:length(dt)
    if dt(i) > 5 %cutoff for infrequent spikes
        nfr = nfr + 1;
    else
        if nfr 
            ifr = [ifr nfr];
        end
        nfr = 0;
    end
end
mifr = mean(ifr);
mdifr = median(ifr);
sifr = std(ifr);
end