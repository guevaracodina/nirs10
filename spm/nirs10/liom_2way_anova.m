function A = liom_2way_anova(cbeta,ccov_beta,s1,s2, nsubj, min_subj)
try
    %number of levels
    nl1 = size(cbeta,2);
    nl2 = size(cbeta,3);
    msk = zeros(nl1,nl2,s1*s2); %a mask for each level
    msk_subj = zeros(nsubj,nl); 
    for i1 = 1:nsubj
        for l1 = 1:nl
            if any(i1 == level_subj{l1})
                index = find(cbeta(i1,:) ~= 0);
                msk(l1,index) = msk(l1,index) + 1;
                clear index
                msk_subj(i1,l1) = 1;
            end
        end
    end
    for l1 = 1:nl
        index_msk{l1} = find(msk(l1,:) > min_subj-1);
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
        dfwithin = dfwithin + msk(l1,ind)-1;
        mean_beta = mean_beta+ avg_beta{l1};
    end
    MSwithin = MSwithin./dfwithin;
    mean_beta = mean_beta/nl;
    SSbetween = zeros(1,lm);
    for l1 = 1:nl
        SSbetween = SSbetween + msk(l1,ind).*(avg_beta{l1}-mean_beta).^2;
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
    
    %store
    A.F = F2;
    A.df = dfwithin2;
    A.mbeta = mean_beta2;
    A.dfbetween = nl-1;
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

end