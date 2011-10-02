function G = liom_group_F(cbeta,s1,s2)
try
    ns = size(cbeta,1);
    %design matrix
    X = ones(1,ns);
    pX = pinv(X);
    b = pX'*cbeta; %simple average
    s = std(cbeta,0,1);   
    
    beta_group = reshape(b,s1,s2);
    std_group = reshape(s,s1,s2);
    t = b./s;
    t(isnan(t)) = 0;
    G.erdf_group = ns-1;
    G.var_bs = [];
    G.beta_group = beta_group;
    G.std_group = std_group;
    G.tmap_group = reshape(t,s1,s2);
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

end