function SCKS = nirs_SCKS_priors_Buxton(SCKS,h)
% returns priors for a hemodynamic dynamic causal model
[pE,pC] = ioi_SCKS_priors_base_case(SCKS.m,5);
pE(end)=0.025; % put a better approximation around effFlow to extract good estimate of covariance

SCKS.name={'Signal_decay', 'Feedback', 'Transit_time',...
        'Exponent', 'Extraction',  'log_signal_ratio' 'eff'};
% %add scaleGx
% pC=blkdiag(pC,zeros(2,2));
%  pE=[pE ; 1; 1];
%  if isfield(SCKS,'scaleGx') && isstruct(SCKS.scaleGx)
%     pE(end-1:end)=pE(end-1:end).*[SCKS.scaleGx.Flow; SCKS.scaleGx.HbR];
%  end

% pEreal = pE;
% %normalise
%  pE=ones(length(pE),1 );
%  %pC=pC(index,index);
%  pC2= pC./ (pEreal*pEreal');
%  pC(isfinite(pC2))=pC2(isfinite(pC2));
%  pC=min(pC,100);
 
%  if isfield(SCKS,'covarianceIni') && SCKS.covarianceIni>0
%     pC=eye(length(pE))*SCKS.covarianceIni;    
%   end
%  %remove null variance
%  dd=diag(pC);
%  pC=pC-diag(dd);
%  dd(dd==0)=0.04;
%  pC=pC+diag(dd);
%  
% SCKS.PS.pEreal = pEreal;
SCKS.PS.pE=pE;
SCKS.PS.pC=pC;
