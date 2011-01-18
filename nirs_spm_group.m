function [tmap_group, erdf_group] = nirs_spm_group(fname_cinterp_SPM_nirs, nsubj, min_subj)

for kk = 1:nsubj
    load(fname_cinterp_SPM_nirs{kk});
    if kk == 1
        msk = zeros(1, cinterp_SPM_nirs.s1 * cinterp_SPM_nirs.s2);
    end
    index = find(cinterp_SPM_nirs.cbeta ~= 0);
    msk(index) = msk(index) + 1;
    clear index
end

s1 = cinterp_SPM_nirs.s1;
s2 = cinterp_SPM_nirs.s2;

index_msk = find(msk > min_subj-1);
indiv_beta = zeros(nsubj, length(index_msk));
indiv_cov  = zeros(nsubj, length(index_msk));

for kk = 1:nsubj
    load(fname_cinterp_SPM_nirs{kk});
    indiv_beta(kk,:) = cinterp_SPM_nirs.cbeta(index_msk);
    indiv_cov(kk,:) = cinterp_SPM_nirs.ccov_beta(index_msk);
end

avg_beta = sum(indiv_beta)./msk(index_msk);

h_wait = waitbar(0, 'Please wait... ');
for kk = 1:length(index_msk)
    waitbar(kk/length(index_msk))
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
    tmp_denum = indiv_cov(index, kk) + var_bs; %%% covariance individual + inter- subject variance
    numer = indiv_beta(index, kk)./tmp_denum;
    numer = sum(numer);
    denum = sqrt(sum(1./tmp_denum));
    tmap(1,kk) = numer./denum;
end
close(h_wait);

clear indiv_beta
clear indiv_cov

tmap_group = zeros(1, s1*s2);
tmap_group(index_msk) = tmap(:);
tmap_group = reshape(tmap_group, s1, s2);

erdf_group = zeros(1, s1*s2);
erdf_group(index_msk) = erdf(:);
erdf_group = reshape(erdf_group, s1, s2);

