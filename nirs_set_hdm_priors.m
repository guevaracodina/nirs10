function M=nirs_set_hdm_priors(M)
%Name of parameters
switch M.O.PhysioModel_Choice
    case 0 %Buxton-Friston
        M.name={'Signal_decay', 'Feedback', 'Transit_time',...
            'Exponent', 'Extraction', 'eff'};
    case 1 %Zheng-Mayhew
        M.name={'Signal_decay', 'Feedback', 'Transit_time',...
            'Exponent', 'Extraction','Vascular_tone' 'Gain_parameter', 'eff'};
    case 2 %Huppert1
        BH.name={'ksr', 'kr', 'ksm', 'km', 'P0r', 'Vw0', 'beta', ...
            'delta', 'HGB', 'Ra0r', 'M0', 'effCMRO', 'effFlow'};
end
%load priors
switch M.O.prior_choice
    case 0 %default
        switch M.O.PhysioModel_Choice
            case 0 %Buxton-Friston
                [pE,pC] =nirs_SCKS_priors_base_case(M.m,5); %keeping all modes
                pE(end)=0.025; % put a better approximation around effFlow to extract good estimate of covariance
            case 1 %Zheng-Mayhew
                [pE,pC] = nirs_priors_ZM(M.m);
                pE(pE==0)=min(pE(pE~=0));%obligatoire que cela vaille pas zéro
            case 2 %Huppert1
                BH.ksr=0.69; % 0.69;
                BH.kr=0.27;% .17;
                BH.ksm=0.96;%; %gros retour lent    kx et km gros(montée rapide)
                BH.km=0.27;%gros retour rapide
                BH.P0r=1; %same as PwO20
                BH.Vw0=0.0388 ;%Vvasc*75%
                BH.beta=2.33;
                BH.delta=0.5;
                BH.HGB=0.16;%fixed, but depends on subject
                BH.Ra0r=1;
                BH.M0=0.0387;
                BH.effCMRO=0.051; %[.051 .006];
                BH.effFlow= 0.4; %[.4 0.2];
                
                BH.Fin0=.0148; %Initial flow, normally derived from model
                
                %%% fixed physical constants
                BH.a=23400;%fixed - from dissociation curve of hemoglobin
                BH.b=150;%fixed - from dissociation curve
                BH.Hn=1.39;%fixed
                BH.rho=1.04; %brain density in g/mL
                BH.omega=44.6429;
                
                BH.alphap=3.9e-5;
                BH.alphat=1.18e-4;
                
                BH.P00=47;
                BH.P0 = BH.P0r*BH.P00;
                BH.Ra00=4.23e3;
                BH.Ra0 = BH.Ra0r*BH.Ra00;
                BH.PaO2=75; % 75 ou 100 à voir
                BH.PtO2=16;
                BH.Rw0=2.53e3;
                BH.Vt=0.9483;
                BH.Va0=0.0388/3; %Vvasc*25%
                BH.zeta1=1;
                BH.d0=300;
                BH.d1=0.5;
                
                BH.HbT0=0.107;
                BH.phiw0=0.75;
                BH.scFlow=1 ;
                BH.scHbR=1;
                
                BH.noCompliance=0;
                if BH.noCompliance==1
                    BH.effCMRO=[0.00000 0];
                    BH.ksr=069;
                    BH.kr=170;
                    BH.Vw0=BH.Fin0;
                end
                for i1=1:length(BH.name)
                    var=BH.(BH.name{i1});
                    pE(i1,1)=var(1);
                    if length(var)>=2
                        pC(i1,i1)=var(2);
                    else
                        pC(i1,i1)=.2^2*var(1)^2; %?variance approximated by a fraction of the square of the mean???
                    end
                end
                M.BH = BH;
            otherwise
        end
    case 1 %finger tapping priors
        switch M.O.PhysioModel_Choice
            case 0 %Buxton-Friston
                pE = [0.6479,0.4097,0.9818,0.3220,0.3410,0.0238]';
                pC = 1e-4*[ 4.6636    0.0059   -0.0190   -0.0006   -0.0018    0.0083
                            0.0059    0.6240    0.0114   -0.0020   -0.0026    0.0353
                           -0.0190    0.0114   17.7676   -0.0129   -0.0300    0.0032
                           -0.0006   -0.0020   -0.0129    0.4051    0.0018   -0.0063
                           -0.0018   -0.0026   -0.0300    0.0018    0.7670    0.0069
                            0.0083    0.0353    0.0032   -0.0063    0.0069    0.0644];
            case 1 %Zheng-Mayhew
                pE = [0.6197    0.3864    1.0364    0.3677    0.3603    0.6769   -0.1340    0.0226]';
                pC = 1e-4*[ 120.0322    5.3362  -13.6570   -0.3334   -0.1305   -1.4070  -14.2291    0.5041
                            5.3362   16.2744    8.6576    0.1753   -0.0469    0.2510    2.0955    0.8045
                          -13.6570    8.6576  492.8406   -0.1529    0.1921   -1.2540   35.0596    0.7939
                           -0.3334    0.1753   -0.1529    8.9465   -2.1060   -0.0226   10.8486   -0.1577
                           -0.1305   -0.0469    0.1921   -2.1060   22.9810   -0.0946    1.9984    0.2372
                           -1.4070    0.2510   -1.2540   -0.0226   -0.0946   98.2691   -7.5497   -0.0007
                          -14.2291    2.0955   35.0596   10.8486    1.9984   -7.5497  473.2026    0.4263
                            0.5041    0.8045    0.7939   -0.1577    0.2372   -0.0007    0.4263    0.1136];
            case 2 %Boas-Huppert
                pE =  [0.4689    0.2910    0.9662    0.2132    0.9684    0.0402    2.8129    0.4941    0.1556    1.0000    0.0386    0.0522    0.2830]';
                pC = 1e-3*[ 11.7197    0.4145    0.0046   -0.0040    0.3165   -0.0145    0.4317   -0.0789   -0.0460    0.0000    0.0146    0.0000    0.5875
                            0.4145    1.4774   -0.0005    0.0042   -0.5180   -0.0005   -1.0353    0.0519    0.0126    0.0001   -0.0008   -0.0000    0.9606
                            0.0046   -0.0005   36.3675    0.0941    0.1200   -0.0018   -0.5311    0.0023    0.0005   -0.0000    0.0039    0.0072   -0.0400
                           -0.0040    0.0042    0.0941    2.2318    0.0233   -0.0024   -0.3204   -0.0002   -0.0026   -0.0000   -0.0000    0.0363   -0.0379
                            0.3165   -0.5180    0.1200    0.0233   12.6975    0.1446    3.7417    0.2071    0.0792    0.0003   -0.0091   -0.0739    0.1436
                           -0.0145   -0.0005   -0.0018   -0.0024    0.1446    0.0577   -0.1018   -0.0060   -0.0017   -0.0000    0.0007    0.0015    0.0431
                            0.4317   -1.0353   -0.5311   -0.3204    3.7417   -0.1018  194.4688    0.2157    0.4062   -0.0001   -0.0579    0.1631    2.7162
                           -0.0789    0.0519    0.0023   -0.0002    0.2071   -0.0060    0.2157    9.9402   -0.0361    0.0000    0.0042   -0.0002    0.0650
                           -0.0460    0.0126    0.0005   -0.0026    0.0792   -0.0017    0.4062   -0.0361    0.9755    0.0000    0.0075   -0.0000    0.0456
                            0.0000    0.0001   -0.0000   -0.0000    0.0003   -0.0000   -0.0001    0.0000    0.0000   40.0000   -0.0000    0.0000   -0.0004
                            0.0146   -0.0008    0.0039   -0.0000   -0.0091    0.0007   -0.0579    0.0042    0.0075   -0.0000    0.0547   -0.0003   -0.0025
                            0.0000   -0.0000    0.0072    0.0363   -0.0739    0.0015    0.1631   -0.0002   -0.0000    0.0000   -0.0003    0.0928    0.0178
                            0.5875    0.9606   -0.0400   -0.0379    0.1436    0.0431    2.7162    0.0650    0.0456   -0.0004   -0.0025    0.0178    1.6575];
             M.BH = BH;
        end
end
M.pE = pE;
M.pC = pC;