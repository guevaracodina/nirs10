function y = nirs_gx(x,u,P,M)
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% exponentiation of hemodynamic state variables
x     = exp(x);
switch M.O.PhysioModel_Choice
    case 0 %Buxton, 
        HbT = (x(3,:)-1); %Assume v = HbT = HbO+ HbR
        HbR = (x(4,:)-1); %Assume q = HbR 
    case 1 %, Zheng Mayhew and Zheng Decay
        HbT = (x(3,:)-1); %Assume v = HbT = HbO+ HbR
        HbR = (x(4,:)-1); %Assume q = HbR 
    case 2 %Huppert1        
        BH = nirs_fillBH(x,P,M.BH); %Does BH need to be recalculated here?
        HbT=(BH.Vw0*(x(3,:)-1 )-0.5*BH.Va0 *(x(2,:)-1))/(BH.Vw0+BH.Va0); %HbTratio %Why additional term?
        HbR=(BH.Vw0*(  (1-BH.S).*x(3,:) - (1-BH.SwO20)) ...
            -0.5*BH.Va0 *(1-BH.SaO2) *(x(2,:)-1))./(BH.Vw0*(1-BH.S)+ BH.Va0*(1-BH.SaO2));
end
y = [];
if M.IC.include_HbR, y = [y; HbR]; end
if M.IC.include_HbT, y = [y; HbT]; end  