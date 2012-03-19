function A = calc_hfgg(cbeta,A,ns,B,nfac)
%sz = size(cbeta);
%cbeta = reshape(cbeta,[sz(1,end-1),s1,s2]);
[min0 imin0] = min(A.F);
[max0 imax0] = max(A.F);
if abs(min0) > abs(max0)
    id0 = imin0;
else
    id0 = imax0;
end

switch nfac
    case 1
        Y = squeeze(cbeta(:,id0));
        S = [];
        WInFacs = [];
        BTFacs = zeros(length(B),ns);
        for i1 =1:length(B)
            BTFacs(i1,B{i1}) = 1;
        end
    case 2
        [ns0 nS0 nC0 np] = size(cbeta);
        Y = squeeze(cbeta(:,:,:,id0));
        Y = Y(:);
        S = zeros(length(Y),1);
        %determine subjects -- same as WInFacs?
        for i1=1:nS0
            for j1=1:nC0
                for k1=1:ns0
                    S(k1+((i1-1)+(j1-1)*nS0)*ns0) = k1;
                end
            end
        end
        %within subject effects -- if requested
        if A.includeSubjectEffects
            WInFacs = B.Xs;
        else
            WInFacs = [];
        end
        %between subject effects
        BTFacs = B.Xbw;
end
        
%Calculate Huynh-Feldt and Greenhouse-Gasser corrections
%Where: at the site of maximal activation?

[A.EpsHF A.EpsList A.EpsGG] = GenCalcHFEps(Y,BTFacs,WInFacs,S);

%For post-hoc contrasts only -- not done yet
switch Z.CorrectionMethod
    case 1 %Huynh-Feldt
        A.erdf = A.erdf;
        A.eidf = A.eidf;
    case 2 %Bonferroni
        p_value = Z.p_value/length(B-1);
    case 3 %Greenhouse-Gasser
        A.erdf = A.erdf;
        A.eidf = A.eidf;
end
%Recommendation: If epsilon is >0.75, use the Huynh-Feldt correction. 
%If epsilon is <0.75, or nothing is known about sphericity at all, 
%use the Greenhouse-Geisser correction

