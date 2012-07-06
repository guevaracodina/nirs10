function y = nirs_gx(x,u,P,M)
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% exponentiation of hemodynamic state variables
x     = exp(x);
switch M.O.PhysioModel_Choice
    case 0 %Buxton, 
        out.HbT = (x(3,:)-1); %Assume v = HbT = HbO+ HbR
        out.HbR = (x(4,:)-1); %Assume q = HbR 
    case 1 %, Zheng Mayhew and Zheng Decay
        out.HbT = (x(3,:)-1); %Assume v = HbT = HbO+ HbR
        out.HbR = (x(4,:)-1); %Assume q = HbR 
    case 2 %Huppert1        
        BH = M.BH;
        BH = nirs_fillBH(x,P,BH);
        out.HbT=(BH.Vw0*(x(3,:)-1 )-.5*BH.Va0 ...
            *(x(2,:)-1))/(BH.Vw0+BH.Va0); %hbTratio
        out.HbR=(BH.Vw0*(  (1-BH.S)*x(3,:) - (1-BH.SwO20)*1  ) ...
            -.5*BH.Va0 *(1-BH.SaO2) *(x(2,:)-1))/(BH.Vw0*(1-BH.S)+ BH.Va0*(1-BH.SaO2));
end
y = [];
if M.O.include_HbR, y = [y; out.HbR]; end
if M.O.include_HbT, y = [y; out.HbT]; end   
