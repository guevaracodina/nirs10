%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machado Alexis 03/05/2010 Ecole polytechnique de Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--INPUTS
% t: nTimePts  by 1
% conc: nTimePts  by nPairs
% X: nTimePts by nCond
% D: nTimePts by nDrifts


%--OUTPUTS
% Structures Betas,Thetas,Sigma2,Modelisation


function [Betas,Thetas,Sigma2,Modelisation] = glm(t,conc,X,D)

nPairs=size(conc,2);
for iPairs=1:nPairs
    display(sprintf('pair number %g',iPairs))
    
    [degh,mh,vh,vs2,degl,ml,vl,yHrf,yDrift]=estime_frh_param(t,...
        conc(:,iPairs),...
        X,...
        D);
    
    %---Betas 
    Betas.betas(iPairs,:)=mh;
    Betas.betasVariances(iPairs,:)=diag(vh);
    Betas.betasFreedom(iPairs,:)=degh;
    
    
    %---Drifts
    Thetas.thetas(iPairs,:)=ml;
    Thetas.thetasVariances(iPairs,:)=diag(vl);
    Thetas.thetasFreedom(iPairs,:)=degl;
   
    
    %---Sigma (noise variance)
    Sigma2.variance(iPairs,:)=vs2;
    
    %---Reconstructed Reponses
    Modelisation.yHrf{iPairs}=yHrf;                                                 
    Modelisation.yDrift{iPairs}=yDrift;                                                 
    Modelisation.yModel(:,iPairs)=yDrift+sum(yHrf,2);                             
    Modelisation.residual(:,iPairs)=conc(:,iPairs)-Modelisation.yModel(:,iPairs);
    
end

Modelisation.t=t;
Modelisation.conc=conc;
