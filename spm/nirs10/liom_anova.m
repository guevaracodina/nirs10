function A = liom_anova(cbeta,ccov_beta,s1,s2,nsubj,min_subj,level_subj)
try %leave ccov_beta dependency for testing purposes...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation #1: standard one-Anova with mask, not using GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Careful -- subjects must be in correct order, none missing
    %number of levels
    nl = length(level_subj);
    for l1=1:nl
        msk{l1} = zeros(1,s1*s2); %a mask for each level
    end
    msk_subj = zeros(nsubj,nl);
    X = zeros(nsubj,nl+1);
    for i1 = 1:nsubj
        for l1 = 1:nl
            if any(i1 == level_subj{l1})
                X(i1,l1) = 1;
                index = find(cbeta(i1,:) ~= 0);
                msk{l1}(index) = msk{l1}(index) + 1;
                %clear index
                msk_subj(i1,l1) = 1;
            end
        end
        X(i1,nl+1) = 1; %constant
    end
    for l1 = 1:nl
        index_msk{l1} = find(msk{l1} > min_subj-1);
    end
    
    ind = [];
    for kk = 1:length(index_msk{1})
        %only look at voxels activated for all groups
        voxel_good = 1;
        for l1 = 2:nl
            if ~any(index_msk{1}(kk)==index_msk{l1})
                voxel_good = 0;
            end
        end
        if voxel_good
            ind = [ind index_msk{1}(kk)];
        end
    end
    lm = length(ind);
    indiv_beta = cbeta(:,ind);
    mean_beta = mean(indiv_beta,1);
    %average for each group
    for l1 = 1:nl
        msk_subj2{l1} = repmat(msk_subj(:,l1),[1 lm]);
        avg_beta{l1} = mean(indiv_beta.*msk_subj2{l1},1);
        rep_avg_beta = repmat(avg_beta{l1},[nsubj 1]);
        SSwithin{l1} = sum(((indiv_beta - rep_avg_beta).*msk_subj2{l1}).^2); %sum of squares within group
    end
    
%     %Residuals
%     res = zeros(nsubj,lm);
%     rep_mean =repmat(mean_beta,[nsubj 1]);
%     for l1=1:nl
%         res = res + (avg_beta{l1}-rep_mean(level_subj{l1},:)); %not quite
%     end
    MSwithin = zeros(1,lm);
    dfwithin = zeros(1,lm);

    for l1 = 1:nl
        MSwithin = MSwithin + SSwithin{l1};
        dfwithin = dfwithin + msk{l1}(ind)-1;
        %mean_beta = mean_beta+ avg_beta{l1};
    end
    MSwithin = MSwithin./dfwithin;
    %mean_beta = mean_beta/nl;
    SSbetween = zeros(1,lm);
    for l1 = 1:nl
        SSbetween = SSbetween + msk{l1}(ind).*(avg_beta{l1}-mean_beta).^2;
    end
    dfbetween = nl-1;
    MSbetween = SSbetween/dfbetween;
    %F statistic, with F(dfbetween,dfwithin) distribution
    F = MSbetween./MSwithin;
    %effect size:
    eta_squared = SSbetween./(SSbetween+MSwithin.*dfwithin); %biased -- it overestimates 
    %how much variance is explained
    %partial eta squared: not recommended
    %partial_eta_squared = SSbetween./(SSbetween+SSwithin+SSerror);
    
    %omega_squared is preferred (but not for repeated measures designs)
    %w2 = (SSeffect - (dfeffect)(MSerror)) / MSerror + SStotal 
    
    %reshape
    F2 = zeros(1,s1*s2);
    F2(ind) = F(:);
    F2 = reshape(F2,s1,s2);
    
    mean_beta2 = zeros(1,s1*s2);
    mean_beta2(ind) = mean_beta(:);
    mean_beta2 = reshape(mean_beta2,s1,s2);
    
    dfwithin2 = zeros(1,s1*s2);
    dfwithin2(ind) = dfwithin(:);
    dfwithin2 = reshape(dfwithin2,s1,s2);
    
    eta_squared2 = zeros(1,s1*s2);
    eta_squared2(ind) = eta_squared(:);
    eta_squared2 = reshape(eta_squared2,s1,s2);
    
%Just for display    
%     MSbetween2 = zeros(1,s1*s2);
%     MSbetween2(ind) = MSbetween(:);
%     MSbetween2 = reshape(MSbetween2,s1,s2);
%     
%     SSbetween2 = zeros(1,s1*s2);
%     SSbetween2(ind) = SSbetween(:);
%     SSbetween2 = reshape(SSbetween2,s1,s2);
%     
%     MSwithin2 = zeros(1,s1*s2);
%     MSwithin2(ind) = MSwithin(:);
%     MSwithin2 = reshape(MSwithin2,s1,s2);
%     
%     SSwithin2 = zeros(1,s1*s2);
%     SSwithin2(ind) = SSwithin{1}(:);
%     SSwithin2 = reshape(SSwithin2,s1,s2);
    
    %store
    A.F = F2;
    A.df = dfwithin2;
    A.mbeta = mean_beta2;
    A.dfbetween = nl-1;
    A.eta_squared = eta_squared2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation #2: using a GLM, over the whole data (no mask)
%needed just to calculate the residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    pX = pinv(X);
    b = pX*cbeta;
    
%just for display
%    b2 = reshape(b(1,:)-b(2,:),s1,s2);
    %residuals
    res = cbeta - X*b;
%     rs = sum(res.^2);
%     
%     Z.includeSubjectEffects = 1;
%     if Z.includeSubjectEffects
%         X0 = ones(nsubj,1);
%         pX0 = pinv(X0);
%         b0 = pX0*cbeta;
%         r0 = cbeta - X0*b0;
%         rs0 = sum(r0.^2);
%     else
%         M = mean(cbeta,1);
%         rs0 = sum((cbeta-repmat(M,size(cbeta,1),1)).^2);
%     end
%     G.erdf = nsubj-rank(X0);
%     G.eidf = rank(X) - rank(X0);
%     F = ((rs0-rs)./rs)*(G.erdf-G.eidf)/G.eidf ;
%     F(isnan(F)) = 0;
%     F = reshape(F,s1,s2);
%     G.Tmap = F; %OK
    
    
%     %Just for display of residuals
%     rs = reshape(rs,s1,s2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of LKC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp_mask = zeros(size(cbeta));
tmp_mask(cbeta == 0) = 1;
mask = sum(tmp_mask,1);
index_group = find(mask == 0);
L2 = calc_LKC(index_group,[s1 s2], reshape(res,[nsubj s1 s2]), 'group');
r = sqrt(L2./pi);
L1 = pi * r;
L0 = 1;
A.LKC = [L0 L1 L2];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for l1 = 1:nl
%         msk2{l1} = zeros(1,s1*s2);
%         msk2{l1}(ind) = msk{l1}(ind);
%         msk2{l1} = reshape(msk2{l1},s1,s2);
%     end
%     
%     for l1 = 1:nl
%         msk3{l1} = zeros(1,s1*s2);
%         msk3{l1}(index_msk{l1}) = msk{l1}(index_msk{l1});
%         msk3{l1} = reshape(msk3{l1},s1,s2);
%     end
%     
%     for l1 = 1:nl
%         msk4{l1} = zeros(1,s1*s2);
%         msk4{l1} = msk{l1};
%         msk4{l1} = reshape(msk4{l1},s1,s2);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation #3: 2-sample t-test, done on whole data -- no mask --
%gives slightly different results from calculation #4, due to difference
%between t-stat and F-stat when the samples have unequal size -- results
%should be the same (i.e. t^2=F, when the samples have equal size)
%SOUNDS WRONG: t^2=F... thus there should be some mistake in here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %if nl==2, try to do a 2 sample t-test
%     %formula: t = (X1-X2)/s, with s^2 = s1^2/n1+s2^2/n2
%     if nl == 2
%         x1b = mean(cbeta(level_subj{1},:),1);
%         x2b = mean(cbeta(level_subj{2},:),1);
%         s1b = std(cbeta(level_subj{1},:),0,1);
%         s2b = std(cbeta(level_subj{2},:),0,1);
%         n1 = length(level_subj{1});
%         n2 = length(level_subj{2});
%         s = (s1b.^2/n1+s2b.^2/n2).^0.5;
%         t = (x1b-x2b)./s;
%         s0 = reshape(s,s1,s2);
%         s1b0 = reshape(s1b,s1,s2);
%         s2b0 = reshape(s2b,s1,s2);
%         t0 = reshape(t,s1,s2);
%         G.Tmap = t0; 
%         A.F = t0.^2;
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation #4: 1-Anova, done on whole data -- no mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %One way anova without masks -- only when using a template for coregistration
%     M = mean(cbeta,1);
%     for l1=1:nl
%         Mi{l1} = mean(cbeta(level_subj{l1},:),1);
%     end
%     BSS = zeros(1,length(M));
%     for l1=1:nl
%         BSS = BSS + length(level_subj{l1})*(Mi{l1}-M).^2;
%         sd{l1} = std(cbeta(level_subj{l1},:),0,1);
%     end
%     dfb = nl-1;
%     BMS = BSS/dfb;
%     WSS = zeros(1,length(M));
%     for l1=1:nl
%         WSS = WSS + (length(level_subj{l1})-1)*sd{l1}.^2;
%     end
%     dfw = size(cbeta,1)-nl;
%     WMS = WSS/dfw;
%     F = BMS./WMS;
%     F0 = reshape(F,s1,s2);
%     TSS = BSS + WSS;
%     TSS0 = reshape(TSS,s1,s2);
%     TSS2 = TSS-sum((cbeta-repmat(M,size(cbeta,1),1)).^2,1);
%     TSS20 = reshape(TSS2,s1,s2); %OK, this is 0 to machine-precision
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end