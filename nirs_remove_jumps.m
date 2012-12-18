function e = nirs_remove_jumps(d,OP)
%n = length(d);
% Sb = 4; %bound on standard deviation
% Nr = 5;
% Mp = 50; %fill gaps of less than 10 seconds (50 bad points)
Sb = OP.Sb;
Nr = OP.Nr;
Mp = OP.Mp;
%Use high pass filtered data
if OP.ubf
    d = ButterHPF(OP.sf,OP.bf,OP.bo,d);
else
    d = d;
end
N = length(d);
e = zeros(1,N);
Wi = round(OP.Wi*OP.sf); %window size in data points
Ns = length(d);
Nwi = floor(Ns/Wi);
Nh = round(Wi/2);
for w0=1:(Nwi+1)
    for h0=1:2 %to do overlap
        if w0 < Nwi+1 || (w0 == Nwi+1 && h0 == 1)
            if w0 == Nwi + 1
                si = N-Wi; %start index
                ei = N; %end index
            else
                si = 1+(w0-1)*Wi+(h0-1)*Nh; %start index
                ei = w0*Wi+(h0-1)*Nh; %end index
            end
            lp = si:ei;
            n = length(lp);
            f = d(lp);
            %treatment
            s = std(f);
            m = median(f);
            i = find(abs(f-m)>Sb*s); %bad points
            if ~isempty(i)
                w0
                h0
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
                g = f;
                g(niI) = [];
                j = interp1(a,g,ni,'cubic'); %interpolate bad points
                b = f;
                b(niI) = j; %replace bad points by interpolated points
                %return value into the output
            else
                b = f;
            end
            e(lp) = b;
        end       
    end
end



%remove jumps options:
%OP.Sb = 4; %number of standard deviations
%OP.Nr = 1/TR; %number of points removed before and after
%OP.Mp = 10/TR; %size of gaps to be filled
%OP.sf = 1/TR; %sampling frequency  fs = NIRS.Cf.dev.fs
%OP.ubf = 1; %use Butterworth filter
%OP.bf = 0.01; %Butterworth HPF cutoff
%OP.bo = 2; %Butterworth order