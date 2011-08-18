function G = liom_group(cbeta,ccov_beta,s1,s2, nsubj, min_subj,FFX)
try 
msk = zeros(1,s1*s2);
for i1 = 1:nsubj
    index = find(cbeta(i1,:) ~= 0);
    msk(index) = msk(index) + 1;
    clear index
end

index_msk = find(msk > min_subj-1);
indiv_beta = zeros(nsubj, length(index_msk));
indiv_cov  = zeros(nsubj, length(index_msk));

for i1 = 1:nsubj
    indiv_beta(i1,:) = cbeta(i1,index_msk);
    indiv_cov(i1,:) = ccov_beta(i1,index_msk);
end

avg_beta = sum(indiv_beta)./msk(index_msk);
tmap = zeros(1,length(index_msk));
beta_group_tmp = zeros(1,length(index_msk));
std_group_tmp = zeros(1,length(index_msk));
erdf = zeros(1,length(index_msk));
for kk = 1:length(index_msk)
    X = ones(msk(index_msk(kk)),1);
    index = find(indiv_beta(:, kk) ~= 0);
    beta_tmp = indiv_beta(index, kk);
    res = beta_tmp - X *avg_beta(kk);
    if length(X) == 1
        var_bs = 0;
        disp('err');
        return
    else
        erdf(1, kk) = msk(index_msk(kk)) - rank(X);
        var_bs = sum(res.^2)./erdf(1,kk); % variance between subjects
    end
    if FFX
        tmp_denum = indiv_cov(index, kk); %no intersubject variance
    else
        tmp_denum = indiv_cov(index, kk) + var_bs; %%% covariance individual + inter- subject variance
    end
    numer = indiv_beta(index, kk)./tmp_denum;
    numer = sum(numer);
    denum = sqrt(sum(1./tmp_denum));
    tmap(1,kk) = numer./denum;
    beta_group_tmp(1,kk) = numer;
    std_group_tmp(1,kk) = denum;
end


tmap_group = zeros(1, s1*s2);
tmap_group(index_msk) = tmap(:);
tmap_group = reshape(tmap_group, s1, s2);

erdf_group = zeros(1, s1*s2);
erdf_group(index_msk) = erdf(:);
erdf_group = reshape(erdf_group, s1, s2);

beta_group = zeros(1, s1*s2);
beta_group(index_msk) = beta_group_tmp(:);
beta_group = reshape(beta_group, s1, s2);

std_group = zeros(1,s1*s2);
std_group(index_msk) = std_group_tmp(:);
std_group = reshape(std_group,s1,s2);

G.tmap_group = tmap_group;
G.erdf_group = erdf_group;
G.var_bs = var_bs;
G.beta_group = beta_group;
G.std_group = std_group;
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
    
end