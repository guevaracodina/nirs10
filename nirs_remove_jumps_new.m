function d = nirs_remove_jumps_new(d,OP)
%d0 = d;
Sb = OP.Sb; %bound on standard deviation
Nr = 10; %OP.Nr;
Mp = 50; %OP.Mp; %size of window
N = length(d);
%obtain a lower bound on the standard deviation
sL = std(d);

%e = zeros(1,N);
for i0=Mp:N-1
    si = i0-Mp+1;
    ei = i0;
    lp = si:ei;
    f = d(lp);
    p = d(ei+1);
    m = median(f);
    s = std(f);
    if abs(p-m)>Sb*max(s,sL)  
        %replace by median of previous Nr points
        q = median(d(ei-Nr+1:ei));
        d(ei+1) = q;
    end
end
