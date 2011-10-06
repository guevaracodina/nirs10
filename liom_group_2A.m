function G = liom_group_2A(cbeta,X,X0,s1,s2,Z)
try
    %reshape cbeta
    [ns0 nS0 nC0 np] = size(cbeta);
    nR = ns0*nS0*nC0;
    cbeta = reshape(cbeta,[nR,np]);
    pX = pinv(X);
    b = pX*cbeta;
    %residuals
    r = cbeta - X*b;
    rs = sum(r.^2);
    if Z.includeSubjectEffects
        pX0 = pinv(X0);
        b0 = pX0*cbeta;
        r0 = cbeta - X0*b0;
        rs0 = sum(r0.^2);
        G.erdf = nR-ns0;
    else
        rs0 = sum(cbeta.^2);
        G.erdf = nR;
    end
    G.eidf = nS0*nC0;
    
    F = ((rs0-rs)./rs)*(G.erdf-G.eidf)/G.eidf ;
    F(isnan(F)) = 0;
    F = reshape(F,s1,s2);
    G.tmap_group = F;
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
end

end