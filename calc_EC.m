function [threshold] = calc_EC(LKC, p_value, stat, df)
% calculating a threshold for LKC-based expected Euler characteristics
% inputs:
% LKC: Lipschitz-Killing curvatures
% p_value: specified p-value
% stat: T or F field
% df: degrees of freedom
% output:
% threshold: calculated threshold according to T or F statistics
try
    pb = 0; %if problem
    L0 = LKC(1);
    L1 = LKC(2);
    L2 = LKC(3);
    if isnan(L1)
        L1 = 0;
        disp('Careful, L1 is NaN, setting it to zero')
    end
    if isnan(L2)
        L2 = 0;
        disp('Careful, L2 is NaN, setting it to zero')
    end
    pOK = 0;
    lc = 0;
    fp = 1;
    dof_cutoff = 340;
    while ~pOK && lc < 2
        lc = lc+1;
        v = df(2);
        if fp
            if v > dof_cutoff
                vOld = v;
                v = dof_cutoff;
                disp(['First pass, # d.o.f. too large (' num2str(vOld) '), reset to ' int2str(v)]);
                fp = 0;
            end
        else
            if v > dof_cutoff
                vOld = v;
                v = min(vOld/2,dof_cutoff);
            end
        end
        cstT = 1/ power(2*pi,3/2) / sqrt(v/2) / gamma(v/2);
        dv2 = (1-v)/2;
        gam2 = gamma((v+1)/2);
        a = 1/(2*pi);
        switch stat
            case 'T'
                t = 0:0.001:11; %Could perhaps increase max value? initially: 7
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
                t = 0:0.01:150; %initially: 70; %Could perhaps increase max value?
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
        if isfinite(sum(p(2:end)))
            pOK = 1;
        else
            if lc == 1
                df(2) = 340; %quick fix to get out of infinite loop -- but only for 1st level of GLM, not for 2nd level
                disp('problem with calc_EC'); % sum(p)= ' sum(p(2:end))]);
                pb = 1;
            else
                disp('Still problem with calc_EC')
            end
        end
    end
    index_th=[];
    dp = 10^(-8);
    while isempty(index_th) == 1;
        index_th=find(p>p_value - dp & p<p_value + dp);
        dp = dp*10;
    end
    if dp>0.1 %not converged, take minimum value for
        disp('EC threshold calculation did not work -- imposing threshold of 3');
        threshold = 3.9;
    else
        threshold = t(index_th(end));
        if pb
            disp(['Problematic threshold: ' stat ' = ' num2str(threshold)]);
        else
            disp(['Calc_EC OK; threshold: ' stat ' = ' num2str(threshold)]);
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1));
    disp('Could not evaluate EC threshold.');
end