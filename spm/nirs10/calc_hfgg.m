function A = calc_hfgg(cbeta,A,ns,B,nfac)
%Inputs - for 2-anova
%cbeta: data, with shape [ns0 nS0 nC0 np] = size(cbeta);
%where ns0: number of subjects, nS0: number of sessions, nC0: number of
%contrasts, np: number of space points, equal to s1*s2;
%A: structure from previously estimated 2-anova
%ns = ns0: number of subjects
%B: information on between factors
%nfac: number of factors
%s1,s2: size of image

%Output: A, to which the following results are added:

%sz = size(cbeta);
%cbeta = reshape(cbeta,[sz(1,end-1),s1,s2]);
if nfac == 2
    A.F = A.Tmap(:);
end

%Find bounds on F map
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
        %S = zeros(length(Y),1);
        %determine subjects -- same as WInFacs?
        %         for i1=1:nS0
        %             for j1=1:nC0
        %                 for k1=1:ns0
        %                     S(k1+((i1-1)+(j1-1)*nS0)*ns0) = k1;
        %                 end
        %             end
        %         end
        S = B.S; % B.WInFacs;
        %within subject effects -- if requested
        %if A.includeSubjectEffects
        %    WInFacs = B.WInFacs;
        %else
        WInFacs = B.WInFacs;
        %end
        %between subject effects
        BTFacs = B.BTFacs;
end

%Calculate Huynh-Feldt and Greenhouse-Gasser corrections
%Where: at the site of maximal activation?

[EpsHF EpsList EpsGG] = GenCalcHFEps(Y,BTFacs,WInFacs,S);
A.Eps.EpsHF = EpsHF;
A.Eps.EpsList = EpsList;
A.Eps.EpsGG = EpsGG;
