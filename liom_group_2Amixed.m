function G = liom_group_2Amixed(cbeta,X,X0,s1,s2,Z)
try
    %reshape cbeta
    [ns0 nB np] = size(cbeta);
    %cbeta = reshape(cbeta,[ns0 nS0 nC0 s1 s2]);
    sX = ns0*nB;
    cbeta = reshape(cbeta,[sX,np]);
    pX = pinv(X);
    b = pX*cbeta;
    %residuals
    res = cbeta - X*b;
    rs = sum(res.^2);
    if Z.includeSubjectEffects
        pX0 = pinv(X0);
        b0 = pX0*cbeta;
        r0 = cbeta - X0*b0;
        rs0 = sum(r0.^2);
        %G.erdf = sX-rank(X0); %sX-ns0;
    else
        rs0 = sum(cbeta.^2);
        %G.erdf = sX-rank(X0); %sX;
    end
    G.erdf = sX-rank(X0);
    %Only valid if ...
    %G.eidf = size(X,2)-size(X0,2); %(nS0-1)*(nC0-1);
    G.eidf = rank(X) - rank(X0); 
    F = ((rs0-rs)./rs)*(G.erdf-G.eidf)/G.eidf ;
    F(isnan(F)) = 0;
    F = reshape(F,s1,s2);
    G.Tmap = F;
       
    % LKC calculation
    tmp_mask = zeros(size(cbeta));
    tmp_mask(cbeta == 0) = 1;
    mask = sum(tmp_mask,1);
    index_group = find(mask == 0);
    [L2] = calc_LKC(index_group,[s1 s2], reshape(res,[size(res,1) s1 s2]), 'group');
    r = sqrt(L2./pi);
    L1 = pi * r;
    L0 = 1;
    G.LKC = [L0 L1 L2];
    
    %effect size:
    %eta_squared = SSbetween./(SSbetween+SSwithin); %biased -- it overestimates 
    %how much variance is explained
    %Omega_squared: do not use for repeated measures! (i.e. not for this
    %2-anova as coded up so far
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end