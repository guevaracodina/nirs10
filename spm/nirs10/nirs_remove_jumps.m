function b = nirs_remove_jumps(d,OP)
n = length(d);
% Sb = 4; %bound on standard deviation
% Nr = 5;
% Mp = 50; %fill gaps of less than 10 seconds (50 bad points)
Sb = OP.Sb;
Nr = OP.Nr;
Mp = OP.Mp;
%Use high pass filtered data
if OP.ubf
    f = ButterHPF(OP.sf,OP.bf,OP.bo,d);
else
    f = d;
end
s = std(f);
m = median(f);
i = find(abs(f-m)>Sb*s); %bad points
%k = find(abs(f-m)<=Sb*s); %good points
%fill gaps in bad points
ni = [];
for i0 = 2:length(i)
    if i(i0) - i(i0-1) <= Mp
        ni = [ni max(1,i(i0-1)-Nr):min(i(i0)+Nr,n)];
    else
        ni = [ni max(1,i(i0)-Nr):min(i(i0)+Nr,n)];
    end
end
ni = unique(ni);
niI = int16(ni);
a = 1:n;
a(niI) = [];
g = d;
g(niI) = [];
j = interp1(a,g,ni,'cubic'); %interpolate bad points
b = d;
b(niI) = j; %replace bad points by interpolated points


%remove jumps options:
%OP.Sb = 4; %number of standard deviations
%OP.Nr = 1/TR; %number of points removed before and after
%OP.Mp = 10/TR; %size of gaps to be filled
%OP.sf = 1/TR; %sampling frequency  fs = NIRS.Cf.dev.fs
%OP.ubf = 1; %use Butterworth filter
%OP.bf = 0.01; %Butterworth HPF cutoff
%OP.bo = 2; %Butterworth order