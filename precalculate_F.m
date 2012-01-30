function TOPO = precalculate_F(SPM,Z,TOPO)
xCon = TOPO.xCon;
TOPO.cSigma = [];
try
    %loop over 2 options
    if Z.GroupMultiSession
        %Sessions are processed as one
        Nsess = 1;
    else
        %Sessions are processed individually
        Nsess = length(SPM.xXn);
    end
    NC = SPM.xY.Cf;
    %load the filtered data
    
    KY = [];
%     dir_done = 0;
    for s1=1:length(SPM.xXn)
        if Z.sessions == 0 || any(s1 == Z.sessions)
            KYfile = SPM.xY.Pf{s1};
            KYtmp = fopen_NIR(KYfile,NC)';
%             if ~dir_done
%                 [dir fil ext] = fileparts(KYfile);
%                 dir_done = 1;
%             end
            if Z.GroupMultiSession
                KY = [KY;KYtmp];
            else
                KY0{s1} = KYtmp;
            end
        end
    end
    
    for c1=1:length(xCon)
        if xCon(c1).STAT == 'F'
            for s1=1:Nsess
                if Z.GroupMultiSession
                    try
                        wX = TOPO.xX.xKXs.X;
                    catch exception
                        disp(exception.identifier);
                        disp(exception.stack(1));
                        disp('Problem with wX -- whitened design matrix');
                        wX = xX.xKXs; %quick fix for GroupMultiSessions - more properly, full xKXs
                        %structure needs to be generated at the subject level from the
                        %multisessions, in order to treat F contrasts
                    end
                    xX = TOPO.xX;
                else
                    wX = SPM.xXn{s1}.xKXs;
                    xX = SPM.xXn{s1};
                end
                if Z.sessions == 0 || any(s1 == Z.sessions)
                    xCon(c1).h = spm_FcUtil('Hsqr',xCon(c1), wX); %Need filtered design matrix
                    switch xX.K.LParam.type
                        case {'hrf', 'Gaussian'}
                            S = xX.K.KL;
                        case 'none'
                            S = speye(nScan);
                    end
                    switch xX.K.HParam.type %not used
                        case 'DCT'
                            S = S - xX.K.X0 * (xX.K.X0' * S);
                            %note NIRS_SPM has a catch if out of memory occurs (- deleted here)
                    end
                    %Calculate modified number of degrees of freedom due to filtering
                    %for the sum of squares difference between the full and reduced models
                    [trMV trMVMV] = approx_trRV(xX.xKXs.X,xX.pKX,S,xCon(c1).c);
                    xCon(c1).trRV = trMV;
                    xCon(c1).trRVRV = trMVMV;
                    xCon(c1).eidf = trMVMV ./ (trMV.^2);
                    xCon(c1).erdf = xX.erdf;
                                       
                    nCont = size(wX, 2);
                    c = xCon(c1).c;
                    c0 = eye(nCont) - c * pinv(c);
                    X0 = wX.X * c0;
                    if ~Z.GroupMultiSession
                        KY = KY0{s1};
                    end
                    %Calculate and save cSigma, for each session and each F-contrast;
                    cSigma = (KY' * KY) - (KY' * X0) * pinv(X0) * KY;
                    %fname = fullfile(dir,['cSigma' gen_num_str(s1,2) 'x' gen_num_str(c1,2) '.nii']);
                    TOPO.cSigma{s1,c1} = cSigma;
                    %fwrite_NIR(fname,cSigma);
                end
            end
        end
    end
    TOPO.xCon = xCon;
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem with precalculation of F statistics (cSigma)');
end