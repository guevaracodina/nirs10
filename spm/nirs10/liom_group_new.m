function G = liom_group_new(cbeta,s1,s2)
try
    ns = size(cbeta,1);
    tmp_mask = zeros(size(cbeta));
    tmp_mask(cbeta == 0) = 1;
    mask = sum(tmp_mask,1);
    index_group = find(mask == 0);
    %design matrix
    X = ones(1,ns);
    pX = pinv(X);
    b = pX'*cbeta; %simple average
    s = std(cbeta,0,1);
    
    beta_group = reshape(b,s1,s2);
    std_group = reshape(s,s1,s2);
    t = b./s;
    t(isnan(t)) = 0;
    res = cbeta - X' * b;
    res = reshape(res,ns,s1,s2);
    % LKC calculation
    [L2] = calc_LKC(index_group,[s1 s2], res, 'group');
    r = sqrt(L2./pi);
    L1 = pi * r;
    L0 = 1;
    
    G.erdf_group = ns-1;
    G.var_bs = [];
    G.beta_group = beta_group;
    G.std_group = std_group;
    G.tmap_group = reshape(t,s1,s2);
    G.LKC = [L0 L1 L2];
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end