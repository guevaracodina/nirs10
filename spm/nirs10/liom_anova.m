function A = liom_anova(cbeta,ccov_beta,s1,s2,nsubj,min_subj,level_subj)
try
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
    indiv_beta = zeros(nsubj,lm);
    indiv_cov  = zeros(nsubj,lm);
    
    for i1 = 1:nsubj
        indiv_beta(i1,:) = cbeta(i1,ind);
        indiv_cov(i1,:) = ccov_beta(i1,ind);
    end
    
    %average for each group
    for l1 = 1:nl
        msk_subj2{l1} = repmat(msk_subj(:,l1),[1 lm]);
        avg_beta{l1} = sum(indiv_beta.*msk_subj2{l1}); %./msk(l1,ind); %xbar_i
        rep_avg_beta = repmat(avg_beta{l1},[nsubj 1]);
        SSwithin{l1} = sum(((indiv_beta - rep_avg_beta).*msk_subj2{l1}).^2); %sum of squares within group
    end
    MSwithin = zeros(1,lm);
    dfwithin = zeros(1,lm);
    mean_beta = zeros(1,lm);
    for l1 = 1:nl
        MSwithin = MSwithin + SSwithin{l1};
        dfwithin = dfwithin + msk{l1}(ind)-1;
        mean_beta = mean_beta+ avg_beta{l1};
    end
    MSwithin = MSwithin./dfwithin;
    mean_beta = mean_beta/nl;
    SSbetween = zeros(1,lm);
    for l1 = 1:nl
        SSbetween = SSbetween + msk{l1}(ind).*(avg_beta{l1}-mean_beta).^2;
    end
    dfbetween = nl-1;
    MSbetween = SSbetween/dfbetween;
    %F statistic, with F(dfbetween,dfwithin) distribution
    F = MSbetween./MSwithin;
    
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
    
    MSbetween2 = zeros(1,s1*s2);
    MSbetween2(ind) = MSbetween(:);
    MSbetween2 = reshape(MSbetween2,s1,s2);
    
    SSbetween2 = zeros(1,s1*s2);
    SSbetween2(ind) = SSbetween(:);
    SSbetween2 = reshape(SSbetween2,s1,s2);
    
    MSwithin2 = zeros(1,s1*s2);
    MSwithin2(ind) = MSwithin(:);
    MSwithin2 = reshape(MSwithin2,s1,s2);
    
    SSwithin2 = zeros(1,s1*s2);
    SSwithin2(ind) = SSwithin{1}(:);
    SSwithin2 = reshape(SSwithin2,s1,s2);
    
    %store
    A.F = F2;
    A.df = dfwithin2;
    A.mbeta = mean_beta2;
    A.dfbetween = nl-1;
    
    % LKC calculation
    pX = pinv(X);
    b = pX*cbeta;
    
    %b2 = zeros(1,s1*s2);
    %b2(ind) = b(1,:);
    b2 = reshape(b(1,:)-b(2,:),s1,s2);
    %residuals
    res = cbeta - X*b;
    rs = sum(res.^2);
    
    Z.includeSubjectEffects = 1;
    if Z.includeSubjectEffects
        X0 = ones(nsubj,1);
        pX0 = pinv(X0);
        b0 = pX0*cbeta;
        r0 = cbeta - X0*b0;
        rs0 = sum(r0.^2);
    else
        rs0 = sum(cbeta.^2);
    end
    G.erdf = nsubj-rank(X0);
    G.eidf = rank(X) - rank(X0);
    F = ((rs0-rs)./rs)*(G.erdf-G.eidf)/G.eidf ;
    F(isnan(F)) = 0;
    F = reshape(F,s1,s2);
    G.Tmap = F;
    
    rs = reshape(rs,s1,s2);
    
    tmp_mask = zeros(size(cbeta));
    tmp_mask(cbeta == 0) = 1;
    mask = sum(tmp_mask,1);
    index_group = find(mask == 0);
    %r = cbeta - repmat(mean_beta2,[nsubj 1 1]) - ;
    [L2] = calc_LKC(index_group,[s1 s2], reshape(res,[nsubj s1 s2]), 'group');
    r = sqrt(L2./pi);
    L1 = pi * r;
    L0 = 1;
    A.LKC = [L0 L1 L2];
    
    for l1 = 1:nl
    msk2{l1} = zeros(1,s1*s2);
    msk2{l1}(ind) = msk{l1}(ind);
    msk2{l1} = reshape(msk2{l1},s1,s2);
    end
    
    for l1 = 1:nl
    msk3{l1} = zeros(1,s1*s2);
    msk3{l1}(index_msk{l1}) = msk{l1}(index_msk{l1});
    msk3{l1} = reshape(msk3{l1},s1,s2);
    end
    
    for l1 = 1:nl
    msk4{l1} = zeros(1,s1*s2);
    msk4{l1} = msk{l1};
    msk4{l1} = reshape(msk4{l1},s1,s2);
    end
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end