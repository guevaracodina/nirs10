function pA = nirs_generate_random_param(S,pE)

pA = repmat(pE',[S.simuIt 1]); %zeros(S.simuIt,length(pE));
ct = 0;
for pE1=1:length(pE)
    if any(S.simuP == 0) || any(pE1 == S.simuP)
        ct = ct+1;
        % Initialize the stream
        mtstream = RandStream('mt19937ar','Seed',pE1);
        RandStream.setDefaultStream(mtstream);
        % Generate the random numbers
        
        switch S.simuPriorDistr
            case 1 % Uniform distribution
                for it1=1:S.simuIt % for each simulation
                    % User-specified priors
                    if ~isempty(S.simuPrior)
                        pA(it1,pE1) = S.simuPrior(ct);
                    end
                    tpA = pA(it1,pE1);
                    % User-specified range of simulated
                    % values for parameters
                    if length(S.simuR) == 1
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR/100),tpA*(1+S.simuR/100));
                    else
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR(ct)/100),tpA*(1+S.simuR(ct)/100));
                    end
                end
                
            case 2 % 2 Gaussians
                for it1=1:round(S.simuIt/2) % for 1st half of simulations
                    % User-specified priors
                    if ~isempty(S.simuPrior1)
                        pA(it1,pE1) = S.simuPrior1(ct);
                    end
                    tpA = pA(it1,pE1);
                    % User-specified range of simulated
                    % values for parameters
                    if length(S.simuR1) == 1
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR1/100),tpA*(1+S.simuR1/100));
                        %normrnd(tpA,S.simuR1/100);
                    else
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR1(ct)/100),tpA*(1+S.simuR1(ct)/100));
                        %pA(it1,pE1) = normrnd(tpA,S.simuR1(ct)/100);
                    end
                end
                
                for it1=round(S.simuIt/2)+1:S.simuIt % for 2nd half of simulations
                    % User-specified priors
                    if ~isempty(S.simuDiffPrior)
                        pA(it1,pE1) = pA(it1-round(S.simuIt/2),pE1)*(1+S.simuDiffPrior(ct)/100);
                    end
                    tpA = pA(it1,pE1);
                    % User-specified range of simulated
                    % values for parameters
                    if length(S.simuR2) == 1
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR2/100),tpA*(1+S.simuR2/100));
                        %pA(it1,pE1) = normrnd(tpA,S.simuR2);
                    else
                        pA(it1,pE1) = unifrnd(tpA*(1-S.simuR2(ct)/100),tpA*(1+S.simuR2(ct)/100));
                        %pA(it1,pE1) = normrnd(tpA,S.simuR2(ct));
                    end
                end
                
        end
    end
end