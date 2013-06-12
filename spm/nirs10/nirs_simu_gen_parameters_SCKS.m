function M = nirs_simu_gen_parameters_SCKS(M)
if M.S.simuOn
    simuR = M.S.simuR;
    simuP = M.S.simuP;
    simuIt = M.S.simuIt;
    simuPrior = M.S.simuPrior;
    pA = repmat(M.pE',[simuIt 1]); %zeros(simuIt,length(pE));
    ct = 0;
    for pE1=1:length(M.pE)
        if simuP == 0 || any(pE1 == simuP)
            ct = ct+1;
            %initialize the stream
            mtstream = RandStream('mt19937ar','Seed',pE1);
            %RandStream.setDefaultStream(mtstream);
            RandStream.setGlobalStream(mtstream);
            for it1=1:simuIt
                %generate the random numbers
                if ~isempty(simuPrior)
                    pA(it1,pE1) = simuPrior(ct);
                end
                tpA = pA(it1,pE1);
                if length(simuR) == 1
                    pA(it1,pE1) = unifrnd(tpA*(1-simuR/100),tpA*(1+simuR/100));
                else
                    pA(it1,pE1) = unifrnd(tpA*(1-simuR(ct)/100),tpA*(1+simuR(ct)/100));
                end
            end
        end
    end
    M.S.pA = pA;
    %stimuli types to include, if including all
    if M.S.simuS == 0
        M.S.simuS = 1:size(M.U.u,2);
    end
end
end