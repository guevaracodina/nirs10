function [f] = nirs_fx(x,u,P,M)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1,:) - vascular signal                                    s
%   x(2,:) - rCBF                                           log(f)
%   x(3,:) - venous volume                                  log(v)
%   x(4,:) - dHb                                            log(q)
% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(3) - transit time                                   (t0)
%   P(4) - exponent for Fout(v)                           (alpha)
%   P(5) - resting oxygen extraction                      (E0)
%
%   P(6:5+m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hdm.m 2495 2008-11-27 12:18:33Z karl $
switch M.O.PhysioModel_Choice
    case 0 %B-F
        x(2:end,:) = exp(x(2:end,:));
        % Fout = f(v) - outflow
        %--------------------------------------------------------------------------
        fv       = x(3,:).^(1/P(4));
        % e = f(f) - oxygen extraction
        %--------------------------------------------------------------------------
        ff       = (1 - (1 - P(5)).^(1./x(2,:)))./P(5);
        % implement differential state equations
        %--------------------------------------------------------------------------
        f(1,:)     = P(6).*u(:,:) - P(1)*x(1,:) - P(2)*(x(2,:) - 1);
        f(2,:)     = x(1,:)./x(2,:);
        f(3,:)     = (x(2,:) - fv)./(P(3).*x(3,:));
        f(4,:)     = (ff.*x(2,:) - fv.*x(4,:)./x(3,:))./(P(3).*x(4,:));
    case 1 %Z-M
        x(2:end,:) = exp(x(2:end,:));    
        % Fout = f(v) - outflow
        %--------------------------------------------------------------------------
        fv       = x(3,:).^(1/P(4))./x(5,:);        
        % e = f(f) - oxygen extraction
        %--------------------------------------------------------------------------
        ff       = (1 - (1 - P(5)).^(1./x(2,:)))/P(5);        
        % implement differential state equations
        %--------------------------------------------------------------------------
        f(1,:)     = P(8)*u(:,:) - P(1)*x(1,:) - P(2)*(x(2,:) - 1);
        f(2,:)     = x(1,:)./x(2,:);
        f(3,:)     = (x(2,:) - fv)./(P(3)*x(3,:));
        f(4,:)     = (ff.*x(2,:) - fv.*x(4,:)./x(3,:))./(P(3)*x(4,:));
        f(5,:)     = (-x(5,:)+exp(-10*P(7).*f(3)))./(10*P(6)*x(5,:)); %%Why 10???
    case 2 %B-H
        x([2 3 5 6 7],:) = exp(x([2 3 5 6 7],:)); 
        
        BH = M.BH;
        BH = nirs_fillBH(x,P,BH);
        %--------------------------------------------------------------------------
        f(1,:)     = -BH.effFlow'*u(1,:)/10 - BH.ksr*x(1,:) - BH.kr*(x(2,:) -1);
        f(2,:)     = x(1,:)/x(2,:);
        f(3,:)     = BH.eta*(BH.fin - BH.fv)/x(3,:);
        f(4,:)     =  BH.effCMRO'*u(1,:)/10 - BH.ksm*x(4,:) - BH.km*(x(5,:) - 1);
        f(5,:)     = x(4,:)/x(5,:);
        f(6,:)     = BH.gamma*   (BH.Gamma0*x(7,:)-x(6,:)-x(5,:)*(BH.Gamma0-1))   /x(6,:);
        f(7,:)     = BH.eta*(  (BH.fin+BH.fv*BH.delta/(1-BH.delta)) *...
            (BH.Chi -BH.S -BH.Upsilon*x(7,:))  - BH.Psi*( BH.Gamma0*x(7,:) -  x(6,:) )  )/...
            (x(3,:)*((1-BH.S)^2*(3*x(7,:)^2+BH.B)/BH.A  + BH.Upsilon )*x(7,:));%        
end