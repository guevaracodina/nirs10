function TOPO = pre_calculate_F(SPM,Z,TOPO)
xCon = TOPO.xCon;
try
    if Z.GroupMultiSessions
        wX = TOPO.xX.xKXs.X;
    else
        wX = SPM.xXn{1}.wX;
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem with wX -- whitened design matrix');
    wX = xX.xKXs; %quick fix for GroupMultiSessions - more properly, full xKXs
    %structure needs to be generated at the subject level from the
    %multisessions, in order to treat F contrasts
end

for c1=1:length(xCon)
    Fdone = 0;
    if xCon(c1).STAT == 'F'
        for s1=1:length(SPM.xXn)
            if Z.sessions == 0 || any(s1 == Z.sessions)
                if ~Fdone
                    xCon(c1).h = spm_FcUtil('Hsqr',xCon(c1), wX); %Need filtered design matrix
                    switch xX.K.LParam.type
                        case {'hrf', 'Gaussian'}
                            S = xX.K.KL;
                        case 'none'
                            S = speye(nScan);
                    end
                    switch xX.K.HParam.type
                        case 'DCT'
                            S = S - xX.K.X0 * (xX.K.X0' * S);
                            %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
                    end
                    %Calculate modified number of degrees of freedom due to filtering
                    %for the sum of squares difference between the full and reduced models
                    trRV2 = approx_trRV(xX.xKXs.X,xX.pKX,S,xCon(c1).c);
                    xCon(c1).trRV = trRV2;                    
                    Fdone = 1;
                end
            end
        end
    end
end

TOPO.xCon = xCon;