function BH = nirs_fillBH(x,P,BH)
%Recover parameters:
for i1=1:length(BH.name)  
   BH.(BH.name{i1}) = P(i1);
end
BH.P0 = BH.P0r*BH.P00;
BH.Ra0 = BH.Ra0r*BH.Ra00;
%Additional parameters in IOI extrait var
%calculate variables for Huppert
BH.lambda=BH.omega*BH.Hn*BH.HGB; %
BH.T0=BH.PtO2;
BH.SaO2=(BH.a/(BH.PaO2^3+BH.b*BH.PaO2)+1)^-1;
BH.CaO2=BH.lambda*BH.SaO2+BH.zeta1*BH.omega*BH.alphap*BH.PaO2;
if BH.noCompliance==1
    BH.CaO2=BH.CaO2/x(2,:);
end

BH.mu=BH.Rw0/(BH.Ra0*BH.Ra00);
BH.fin      = 1/x(2,:) *(1+BH.mu*(1-x(3,:)^BH.beta)); %%% revoir le -

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
BH.fv       = x(3,:)^(2+BH.beta);

% fin = inflow
BH.fin      = 1/x(2,:) *(1+BH.mu*(1-x(3,:)^BH.beta)); %%% revoir le -

global Fin0 dataOLD
data= [BH.alphap  BH.P0*BH.P00 BH.M0 BH.zeta1 BH.delta BH.SaO2 BH.Vt BH.omega BH.alphap BH.lambda BH.rho]  ;
if  isempty(Fin0) || any(abs((data-dataOLD)./data)>0.000001)
    Fin0=1/(1+BH.delta/(1-BH.delta))*(BH.SaO2-(BH.a/(BH.P0.^3+BH.b*BH.P0)+1).^-1 + ...
        BH.zeta1*BH.omega*BH.alphap*BH.PaO2/BH.lambda-BH.zeta1*BH.omega*BH.alphap*BH.P0./BH.lambda)^(-1)*BH.rho*BH.M0*BH.Vt/BH.lambda; %changé 24 fev 11 (pas d'effet car on calcule cette valeur lorsque fin=fv=1, mais plus logique)
    dataOLD=data;
end
BH.Fin0=Fin0;
BH.SwO20=(BH.a/( BH.P0^3+BH.b*BH.P0 )+1)^-1;
BH.Gamma0=BH.P0/BH.T0;
BH.K=BH.Vt*BH.rho*BH.M0/BH.T0/(BH.Gamma0-1);%eq36
BH.gamma=BH.K/BH.omega/BH.alphat/BH.Vt;%not used
BH.eta=BH.Fin0/BH.Vw0;%not used
BH.A=BH.a/BH.P0^3;
BH.B=BH.b/BH.P0^2;
BH.Upsilon=BH.zeta1*BH.omega*BH.alphap*BH.P0/BH.lambda;
BH.Chi=BH.CaO2/BH.lambda;
BH.Psi=BH.K*BH.T0/BH.Fin0/BH.lambda;
BH.transit=BH.Vw0/BH.Fin0;
BH.Grubb=1/(2+BH.beta);
%from dissosciation curve
BH.S=1/( BH.A/(x(7,:)^ 3+BH.B*x(7,:) )+1);

if BH.noCompliance==1
    BH.fin      = 1/x(2,:) ; 
    BH.fv       = BH.fin;
    BH.fin      = 1 ; 
    BH.fv       = 1;
    BH.CaO2=1/x(2,:);    
end