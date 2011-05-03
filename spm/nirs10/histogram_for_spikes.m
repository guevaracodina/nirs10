function histogram_for_spikes

t = spm_select([1 Inf],'mat');
dt = [];

for f=1:size(t)
    clear names onsets durations
    load(t(f,:));
    ons = onsets{1};
    d = diff(ons);
    dt = [dt d];
    figure;
    hist(d,80);
            
end
figure;
hist(dt,80)
%values for SD 5 sessions
xs = length(find(dt<1))/length(dt); %45%
xm = length(find(dt<5))/length(dt)-xs; %32%
s = std(dt); %10s
md = median(dt); %1.2s
m = mean(dt); %5.1s
xl = length(find(dt>10))/length(dt); % 
%values for EA

%values for ASC

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