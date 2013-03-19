function G = liom_2samples_t_test(Z,cbeta,nx,ny)
try
    [ns np] = size(cbeta);
    tmp_mask = zeros(size(cbeta));
    tmp_mask(cbeta == 0) = 1;
    %mask counts the number of missing subjects at each pixel
    mask = sum(tmp_mask,1);   
    index_group = find(mask == 0); %intersection of all the subjects
    
    n1 = length(Z.group1);
    n2 = length(Z.group2);
    %design matrix
    X = ones(1,ns);
    X(Z.group2) = -1;
    pX = pinv(X);
    b = pX'*cbeta; %b = m1-m2
    res = cbeta - X' * b;
    res = reshape(res,ns,nx,ny);
    % LKC calculation
    L2 = calc_LKC(index_group,[nx ny], res, 'group');
    r = sqrt(L2./pi);
    L1 = pi * r;
    L0 = 1;
    %means and standard deviations
    %s = std(cbeta,0,1);
    s1 = std(cbeta(Z.group1,:),0,1);
    s2 = std(cbeta(Z.group2,:),0,1);
    m1 = mean(cbeta(Z.group1,:),1);
    m2 = mean(cbeta(Z.group2,:),1);
    
    if Z.eq_variance
        s12 = (((n1-1)*s1.^2+(n2-1)*s2.^2)/(n1+n2-2))^(1/2); %should be identical to s = std_group
        t = (m1-m2)./(s12*(1/n1+1/n2)^(1/2));
        MM = max(mask(:));
        df = mask*(MM-1)/MM;
        erdf_group = ns-1;
    else
        s12 = (s1.^2/n1+s2.^2/n2).^(1/2);
        t = (m1-m2)./s12;
        df = (s1.^2/n1 + s2.^2/n2).^2./( (s1.^2/n1).^2/(n1-1) + (s2.^2/n2).^2/(n2-1));
        df(isnan(df)) = 0;
        erdf_group = max(df(:));
    end
    t(isnan(t)) = 0;
   
%     s1_2stt = reshape(s1,nx,ny);
%     s2_2stt = reshape(s2,nx,ny);
%     m1_2stt = reshape(m1,nx,ny);
%     m2_2stt = reshape(m2,nx,ny);
    s12 = reshape(s12,nx,ny);
    beta_group = reshape(m1-m2,nx,ny);
    G.ns = ns;
    G.erdf_group = erdf_group; %ns-1;
    G.var_bs = [];
    G.beta_group = beta_group;
    G.std_group = s12; %std_group;
    G.tmap_group = reshape(t,nx,ny);
    G.LKC = [L0 L1 L2];
    G.df = reshape(df,nx,ny);
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end