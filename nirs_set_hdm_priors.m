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
                %pE(pE==0)=min(pE(pE~=0));%obligatoire que cela vaille pas zéro
            case 2 %Huppert1
                %Prior values for the 13 parameters
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
                BH.effCMRO=0.0051; %[.051 .006];
                BH.effFlow= 0.04; %[.4 0.2];
                %All the quantities below should be either
                %physical/physiological constants, or derived
                
                BH.Fin0=0.0148; %Initial flow, normally derived from model
                
                %%% fixed physical constants
                BH.a=23400;%fixed - from dissociation curve of hemoglobin (in (mmHg)^3)
                BH.b=150;%fixed - from dissociation curve (in (mmHg)^2)
                BH.Hn=1.39;%fixed, Huefner's constant
                BH.rho=1.04; %Considered fixed: brain density in g/mL
                BH.omega=44.6429; %Fixed: Conversion factor from ideal gas law, in umol O2/mL O2
                
                BH.alphap=3.9e-5; %Solubility of O2 in plasma in mL O2/mmHg . mL blood
                BH.alphat=1.18e-4; %Solubility of O2 in tissue -- different from Simon's article? 
                
                BH.P00=47;%?
                BH.P0 = BH.P0r*BH.P00;
                BH.Ra00=4.23e3;
                BH.Ra0 = BH.Ra0r*BH.Ra00;
                BH.PaO2=75; % 75 ou 100 à voir
                BH.SaO2=(BH.a/(BH.PaO2^3+BH.b*BH.PaO2)+1)^-1;
                BH.T0=16; %PtO2
                BH.Rw0=2.53e3;
                BH.Vt=0.9483;
                BH.Va0=0.0388/3; %Vvasc*25%
%                 BH.zeta1=1;
%                 BH.d0=300;
%                 BH.d1=0.5;
                
                BH.HbT0=0.107;
                BH.phiw0=0.75;
%                 BH.scFlow=1 ;
%                 BH.scHbR=1;
%                 
%                 BH.noCompliance=0;
%                 if BH.noCompliance==1
%                     BH.effCMRO=[0.00000 0];
%                     BH.ksr=069;
%                     BH.kr=170;
%                     BH.Vw0=BH.Fin0;
%                 end
                for i1=1:length(BH.name)
                    var=BH.(BH.name{i1});
                    pE(i1,1)=var(1);
                    if length(var)>=2
                        pC(i1,i1)=var(2);
                    else
                        pC(i1,i1)=.2^2*var(1)^2; %?variance approximated by a 20%^2 fraction of the square of the mean???
                    end
                end
                M.BH = BH;
            otherwise
        end
    case 1 %finger tapping priors
        switch M.O.PhysioModel_Choice
            case 0 %Buxton-Friston
                pE = [0.6479,0.4097,0.9818,0.3220,0.3410, 0.0238]';
                pC = 1e-4*[ 4.6636    0.0059   -0.0190   -0.0006   -0.0018
                            0.0059    0.6240    0.0114   -0.0020   -0.0026
                            -0.0190    0.0114   17.7676   -0.0129   -0.0300
                            -0.0006   -0.0020   -0.0129    0.4051    0.0018
                            -0.0018   -0.0026   -0.0300    0.0018    0.7670  ];
                pC    = blkdiag(pC,0.2^2*pE(end)^2);
            case 1 %Zheng-Mayhew
                %pE = [0.6197    0.3864    1.0364    0.3677    0.3603    0.6769   -0.1340 0.0226]';
                pE = [0.5612    0.3750    1.0302    0.3668    0.3605    0.6779   -0.1471    0.0358]'; %Don't forget the prime!
                pC = 1e-4*[ 109.4961    6.0781  -18.0980   -0.5301   -0.2014   -1.8868  -18.8710    
                    6.0781   15.4633   11.5976    0.2812   -0.0498    0.2296    3.1991   
                    -18.0980   11.5976  452.3288   -0.2350    0.4392   -2.4493   45.1437    
                    -0.5301    0.2812   -0.2350    7.9794   -2.6269   -0.1063   11.3792   
                    -0.2014   -0.0498    0.4392   -2.6269   22.6507   -0.1495    2.9548   
                    -1.8868    0.2296   -2.4493   -0.1063   -0.1495   97.1279   -8.3872   
                    -18.8710    3.1991   45.1437   11.3792    2.9548   -8.3872  415.4385 ];
                
%                 pC = 1e-4*[ 120.0322    5.3362  -13.6570   -0.3334   -0.1305   -1.4070  -14.2291  
%                             5.3362   16.2744    8.6576    0.1753   -0.0469    0.2510    2.0955   
%                           -13.6570    8.6576  492.8406   -0.1529    0.1921   -1.2540   35.0596   
%                            -0.3334    0.1753   -0.1529    8.9465   -2.1060   -0.0226   10.8486   
%                            -0.1305   -0.0469    0.1921   -2.1060   22.9810   -0.0946    1.9984    
%                            -1.4070    0.2510   -1.2540   -0.0226   -0.0946   98.2691   -7.5497  
%                           -14.2291    2.0955   35.0596   10.8486    1.9984   -7.5497  473.2026  ];
                %Append for vascular tone and gain parameters
                %pC = blkdiag(pC,diag([0.01 0.1]));
                 pC    = blkdiag(pC,10); %0.2^2*pE(end)^2);
            case 2 %Boas-Huppert
                pE =  [0.4689    0.2910    0.9662    0.2132    0.9684    0.0402    2.8129    0.4941    0.1556    1.0000    0.0386    0.00522    0.02830]';
%                 pC = 1e-3*[ 11.7197    0.4145    0.0046   -0.0040    0.3165   -0.0145    0.4317   -0.0789   -0.0460    0.0000    0.0146   
%                             0.4145    1.4774   -0.0005    0.0042   -0.5180   -0.0005   -1.0353    0.0519    0.0126    0.0001   -0.0008   
%                             0.0046   -0.0005   36.3675    0.0941    0.1200   -0.0018   -0.5311    0.0023    0.0005   -0.0000    0.0039    
%                            -0.0040    0.0042    0.0941    2.2318    0.0233   -0.0024   -0.3204   -0.0002   -0.0026   -0.0000   -0.0000    
%                             0.3165   -0.5180    0.1200    0.0233   1e-1*12.6975    0.1446    3.7417    0.2071    0.0792    0.0003   -0.0091   
%                            -0.0145   -0.0005   -0.0018   -0.0024    0.1446    0.0577   -0.1018   -0.0060   -0.0017   -0.0000    0.0007   
%                             0.4317   -1.0353   -0.5311   -0.3204    3.7417   -0.1018  1e-2*194.4688    0.2157    0.4062   -0.0001   -0.0579   
%                            -0.0789    0.0519    0.0023   -0.0002    0.2071   -0.0060    0.2157    9.9402   -0.0361    0.0000    0.0042   
%                            -0.0460    0.0126    0.0005   -0.0026    0.0792   -0.0017    0.4062   -0.0361    0.9755    0.0000    0.0075  
%                             0.0000    0.0001   -0.0000   -0.0000    0.0003   -0.0000   -0.0001    0.0000    0.0000   40.0000   -0.0000   
%                             0.0146   -0.0008    0.0039   -0.0000   -0.0091    0.0007   -0.0579    0.0042    0.0075   -0.0000    0.0547  ];
%                 pC    = blkdiag(pC,0.2^2*pE(end-1)^2,0.2^2*pE(end)^2);


                %Put parameters in BH structure
                pC = zeros(length(BH.name));
                for i1=1:length(BH.name)
                     BH.(BH.name{i1}) = pE(i1); 
                     pC(i1,i1)=0.5^2*pE(i1)^2;
                end
                %Physical constants and derived quantities: should be
                %recalculated
                %BH.Fin0=0.014; %Initial flow, normally derived from model
                
                %%% fixed physical constants
                BH.a=23400;%fixed - from dissociation curve of hemoglobin (in (mmHg)^3) %A in Simon's article
                BH.b=150;%fixed - from dissociation curve (in (mmHg)^2) %B, in Simon's article
                BH.Hn=1.39;%fixed, Huefner's constant
                BH.rho=1.04; %Considered fixed: brain density in g/mL -- rhob in article
                BH.omega=44.6429; %Fixed: Conversion factor from ideal gas law, in umol O2/mL O2
                
                BH.alphap=3.9e-5; %Solubility of O2 in plasma in mL O2/mmHg . mL blood; instead of gammap?
                BH.alphat=1.18e-4; %Solubility of O2 in tissue -- different from Simon's article? ; instead of gammat?
                
                BH.P00=47; %Peu important, puisque valeur effective redéfinie par le prior 
                BH.Ra00=4.23e3; %Not consistent with Simon's article, with Rw00/(Ra00+Rw00) = 0.3;
                BH.PaO2=75; % 75 ou 100 à voir
                BH.SaO2=(BH.a/(BH.PaO2^3+BH.b*BH.PaO2)+1)^-1;
                BH.T0=16; %PtO2
                BH.Rw0=2.53e3; %Peu important, puisque valeur effective redéfinie par le prior 
                BH.Vt=0.9483;
                BH.Va0=0.0388/3; %Vvasc*25% -- why so small?
                %BH.zeta1=1; %some adjustment factor, set to 1 (no effect, on Fin)                        
                M.BH = BH;
        end
end
M.pE = pE;
M.pC = pC;