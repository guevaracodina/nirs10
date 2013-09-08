function BH = nirs_fillBH(x,P,BH)
%Recover parameters:
for i1=1:length(BH.name)  
   BH.(BH.name{i1}) = P(i1);
end
BH.P0 = BH.P0r*BH.P00; %update #1 from variable P0r
BH.Ra0 = BH.Ra0r*BH.Ra00; %update #2 from variable Ra0r
%Additional parameters in IOI extrait var
%calculate variables for Huppert
BH.lambda=BH.omega*BH.Hn*BH.HGB; %update #3 from variable HGB
%BH.T0=BH.PtO2; %This is a constant

BH.A=BH.a/BH.P0^3; %Normalisation of a by P0, %In Simon's article, A = k_A 
BH.B=BH.b/BH.P0^2; %Normalisation of b by P0  %In Simon's article, B = k_B
%from dissociation curve -- 
BH.S=1./(BH.A./(x(7,:).^3+BH.B*x(7,:))+1); %OK
%BH.SwO20=(BH.a/( BH.P0^3+BH.b*BH.P0 )+1)^-1; %this used in nirs_gx for HbR
BH.SwO20=(BH.A/(1+BH.B)+1)^-1; %this is used in nirs_gx for HbR

BH.CaO2=BH.lambda*BH.SaO2+BH.omega*BH.alphap*BH.PaO2; %OK
% if BH.noCompliance==1
%     BH.CaO2=BH.CaO2/x(2,:);
% end
BH.mu=BH.Rw0/BH.Ra0; %Corrected, was BH.mu=BH.Rw0/(BH.Ra0*BH.Ra00); mu is k6 in Simon's article

BH.fv = x(3,:).^(2+BH.beta); %OK
% fin = inflow
BH.fin = (1+BH.mu*(1-x(3,:).^BH.beta))./x(2,:); %update #4 from variable beta -- OK

%Fin0 is a constrained variable: it needs to be recalculated as needed
% global Fin0 dataOLD
% data= [BH.alphap BH.P0*BH.P00 BH.M0 BH.delta BH.SaO2 BH.Vt BH.omega BH.alphap BH.lambda BH.rho]; %11 parameters
% if  isempty(Fin0) || any(abs((data-dataOLD)./data)>0.000001)
%     Fin0=1/(1+BH.delta/(1-BH.delta))*(BH.SaO2-(BH.a/(BH.P0.^3+BH.b*BH.P0)+1).^-1 + ...
%         BH.omega*BH.alphap*BH.PaO2/BH.lambda-BH.omega*BH.alphap*BH.P0./BH.lambda)^(-1)*BH.rho*BH.M0*BH.Vt/BH.lambda; %changé 24 fev 11 (pas d'effet car on calcule cette valeur lorsque fin=fv=1, mais plus logique)
%     dataOLD=data;
% end
% BH.Fin0=Fin0;
BH.Upsilon = BH.omega*BH.alphap*BH.P0/BH.lambda; %In Simon's article,Upsilon = k4; alphap ou gammap???
BH.Cw0 = BH.lambda*(BH.Upsilon + BH.SwO20);
BH.Fin0=(1-BH.delta)*BH.rho*BH.M0*BH.Vt/(BH.CaO2-BH.Cw0); %Wrong sign in Simon's article and missing lambda in definition of Cw0

BH.Gamma0=BH.P0/BH.T0; %In Simon's article, Gamm0=k7 
BH.K=BH.Vt*BH.rho*BH.M0/(BH.T0*(BH.Gamma0-1));%eq36 rho=rho_b,
BH.gamma=BH.K/(BH.omega*BH.alphat*BH.Vt);%used in f(6,:) %In Simon's article, gamma = k1.
BH.eta=BH.Fin0/BH.Vw0; %used in f(3,:) %In Simon's article, eta = k2
BH.Chi=BH.CaO2/BH.lambda; %In Simon's article,Chi = k5
BH.Psi=BH.K*BH.T0/(BH.Fin0*BH.lambda); %In Simon's article,Psi = k3 -- OK
%BH.transit=1/BH.eta=BH.Vw0/BH.Fin0; %Definition
%BH.Grubb=1/(2+BH.beta); %Definition



% 
% if BH.noCompliance==1
%     BH.fin      = 1/x(2,:) ; 
%     BH.fv       = BH.fin;
%     BH.fin      = 1 ; 
%     BH.fv       = 1;
%     BH.CaO2=1/x(2,:);    
% end