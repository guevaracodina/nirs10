% =========================================================================
% estime_frh_param.m
% 13/12/2005
% Guillaume Marrelec 
% modified by Machado Alexis 03/05/2010 Ecole polytechnique de
% Montreal
%
% Estimation de l'intensite de la reponse hemodynamique par maximum de
% vraisemblance sur un modele parametrique (reponse canonique SPM)
%
%--INPUT
% t             nTimepts by 1
% X             Design Matrix (nTimepts by nCond)
% y             Data to estimate ( Concentration) nTimePts by 1
% D             Design Matrix for nuisancess(nTimePts by non correlated nDrifts) 

%--OUTPUT
% degh          degree of freedom
% mh            estimator of betas dim ( 1 by nCond)
% vh            covariance of betas    ( 1 by nCond)
% degs2         degree of freedom for the noise  dim scalar
% vs2           scale matrix for the noise
% degl          derive degree of freedom
% ml            derive mean
% Vl            derive variance
% yHrf          estimated signals of interest nTimePts by nCond
% yDrift        signal of drifts  nTimepts by 1
% 
% =========================================================================






function [degh,mh,vh,vs2,degl,ml,Vl,yHrf,yDrift]= ...
    estime_frh_param(t,y,X,D,prevent_value,perform_variance)


% PARAMETERS
nTimePts=length(t);
if (nargin<=4), prevent_value = 5; end                                           % prevent bad conditionning for inversion of D'*D
if (nargin<=5), perform_variance = true; end                                     % Compute variance of drifts 




%--COMPUTATION OF PROJECTOR J
if (size(D,2) > size(D,1)-prevent_value)
    D = D(:,1:size(D,1)-prevent_value);                                         % on enleve les derniers drifts
end
J=-D/(D'*D)*D';                                                              % J : nTimePts by nTimePts
for i=1:nTimePts
    J(i,i)=J(i,i)+1;
end


% ESTIMATION -----------------------------
M=X'*J*X;
%invM=inv(M);                                                                    % Attention on peut rencontrer probleme de singularite
yJy=y'*J*y;


% Estimation of betas (mean of posterior distribution)
mh=M\(X'*J*y);
if ~isreal(mh)
    display( 'attention on a des valeurs complexes')
end
clear J;

% Degre of freedom
nDrift=size(D,2);                                                               % nb de drift
nCond=size(X,2);                                                                % nb de conditions expérimentales
degh=nTimePts-(nDrift+1)-(nCond+1);
degs2=degh;


% Variance of noise
vs2=(yJy-mh'*M*mh)/degs2;
% Variance of betas
vh=vs2/M;

% Estimated signal
for iCond = 1:nCond
    yHrf(:,iCond)=X(:,iCond)*mh(iCond);                                         % yHrf:nTimePts by nCond
end


% --- NUISANCES --
%--COMPUTATION OF PROJECTOR K
K=-X/(X'*X)*X'; 
for i=1:nTimePts
    K(i,i)=K(i,i)+1;
end
yKy=y'*K*y;

% Degree of freedom
degl=degh;

% estimation of thetas
ml=(D'*K*D)\(D'*K*y); 

% Variance of thetas
if perform_variance
    Vl=(yKy-ml'*(D'*K*D)*ml)*inv(D'*K*D)/degl;  
else
    Vl=0;
end

% Nuisance signal
yDrift=(D*ml);                                                                %  yDrifts:nTimePts by 1


