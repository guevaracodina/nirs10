function Y = ButterLPF(fs,cutoff,FilterOrder,Y)
Wn=cutoff*2/fs;  % normalised cutoff frequency
[fb,fa]=butter(FilterOrder,Wn);            % butterworth filter
Y=filtfilt(fb,fa,Y); 
end