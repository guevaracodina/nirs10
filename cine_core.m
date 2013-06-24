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
                chO = W.ch_HbO;
                chR = W.ch_HbR;
        end
        hb = get_chromophore(h1);
        %Interpolation of beta and covbeta
        if h1 == 3
            Q.ibeta = (W.beta(:,chO)+W.beta(:,chR)) * Q.B;
            ivar_temp = sum(Q.B.* (diag(W.covbeta(1,chO)+W.covbeta(1,chR)) * Q.B),1);
            for i0 = 1:size(W.covbeta,1)
                Q.ivar(i0,:) = ivar_temp;
            end
        else
            Q.ibeta = W.beta(:,ch) * Q.B;
            ivar_temp = sum(Q.B.* (diag(W.covbeta(1,ch)) * Q.B),1);
            for i0 = 1:size(W.covbeta,1)
                Q.ivar(i0,:) = ivar_temp;
            end
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
