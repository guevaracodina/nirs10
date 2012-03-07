function [threshold] = calc_EC(LKC, p_value, stat, df)
% calculating a threshold for LKC-based expected Euler characteristics 
% inputs: 
% LKC: Lipschitz-Killing curvatures 
% p_value: specified p-value
% stat: T or F field 
% df: degrees of freedom 
% output: 
% threshold: calculated threshold according to T or F statistics 
L0 = LKC(1);
L1 = LKC(2);
L2 = LKC(3);
v = df(2);
cstT = 1/ power(2*pi,3/2) / sqrt(v/2) / gamma(v/2);
dv2 = (1-v)/2;
gam2 = gamma((v+1)/2);
a = 1/(2*pi);
switch stat
    case 'T'
        t = 0:0.001:7;
        p = zeros(1,length(t));
        for i=1:length(t)
            ti = t(i);
            pti2 = 1+ti^2/v;
            p0 = 1-spm_Tcdf(ti,v);
            p1 = power(pti2, dv2)*a;
            p2 = gam2*ti*power(pti2, dv2) * cstT ;
            p(i) = p0*L0 + p1*L1 + p2*L2;
        end
    case 'F'
        k = df(1);        
        b = gammaln(v/2) + gammaln(k/2);
        sqa = a^(1/2)*exp(gammaln((v+k-1)/2)-b)*2^(1/2);
        acst = a*exp(gammaln((v+k-2)/2)-b);
        k1 = 1/2*(k-1);
        k2 = 1/2*(k-2);
        kv = -1/2*(v+k-2);
        kdv = k/v;
        t = 0:0.01:70;
        p = zeros(1,length(t));
        for i=1:length(t)
            ti = t(i);
            kdti = kdv*ti;
            ktiv = (1+kdti).^kv; 
            p0 = 1 - spm_Fcdf(ti,[k,v]);
            p1 = sqa*kdti.^k1.*ktiv;
            p2 = acst*kdti.^k2.*ktiv.*((v-1)*kdti-(k-1));
            p(i) = p0*L0 + p1*L1 + p2*L2;
        end        
end
index_th=[];
dp = 10^(-8);
while isempty(index_th) == 1;
    index_th=find(p>p_value - dp & p<p_value + dp);
    dp = dp*10;
end
threshold = t(index_th(end));