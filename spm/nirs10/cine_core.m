function CINE = cine_core(Z,W,CINE,f1,u1,k0)
Q = W.Q;
try    
    for h1=Z.select_chromophore  
        switch h1
            case 1
                ch = W.ch_HbO;
            case 2
                ch = W.ch_HbR;
            case 3
                ch = W.ch_HbT;
        end
        hb = get_chromophore(h1);
        %Interpolation of beta and covbeta
        Q.ibeta = W.beta(:,ch) * Q.B;
        for i0 = 1:size(W.covbeta,1)
            Q.ivar(i0,:) = sum(Q.B.* (diag(W.covbeta(i0,ch)) * Q.B),1);
        end
        C = loop_cine(Q,W,Z);
        interpolated_cine(Z,W,C,Q,f1,Z.erdf_default,hb,u1,k0);
        %fill CINE
        %CINE = fill_TOPO(CINE,C,W.side_hemi,f1,hb,thz);        
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp(['Problem with ' hb ' CINE -- setting up figures or saving them']);
end
