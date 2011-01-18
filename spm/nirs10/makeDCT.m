
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machado Alexis 03/05/2010 Ecole polytechnique de Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% makeDesign.m

% INPUT
% t:             nTimePtsData by 1
% tStimuli:      nTimePtsStim by 1
% stimuli        nTimePtsStim by nCond
% degre          Polynom degre
% nCos           number of cosinus
% treshold_corr

% OUTPUT
% X: nTimePtsdata by nCond+(nDeriv*nCond)
% D: nTimePtsData by nDrifts non correlated
% correlationInfos: matrices des indices des drifts correlees


function [D correlated_infos] = makeDCT(t,fmax,degre,threshold_corr,X)

%---PARAMETERS OF DATA AND STIMULI
nTimePtsData=length(t);
dtData=t(2)-t(1);

Dphysio = [];
start_correlation = 2;

%------CREATION OF DRIFT MATRIX ---
nCos = fix(2*fmax*nTimePtsData*dtData+1); % le +1 est la car il y a creation un regresseur cst ds spm_dctmtx
fRealMax = (nCos-1)/(2*nTimePtsData*dtData);
display(sprintf('fmax for cosinusoides is %f: ',fRealMax));
Dini = create_base(nTimePtsData,t,degre,nCos,Dphysio);                          % Dini: nTimePtsData by nDrifts total   



%-----CORRELATION ANALYSIS
corr=zeros(size(Dini,2),size(X,2));
correlated_drifts = [];
corr_values=[];

nReg = size(X,2);

for iReg=1:nReg 
    coeff = sum(X(:,iReg).^2)/length(X(:,iReg));                                % E[X^2]
    X_norm = X(:,iReg)/sqrt(coeff);                                             % Normalization for correlation analysis
    %---Perform correlation between paradigm and drifts
    nDrifts=size(Dini,2);
    for iDrifts=start_correlation:nDrifts                                    
        tmp2 = corrcoef(X_norm,Dini(:,iDrifts));
        corr(iDrifts,iReg) = tmp2(1,2);
        %---Check correlation coefficients upon an arbitrary threshold
        if (abs(corr(iDrifts,iReg)) > threshold_corr)
            correlated_drifts = cat(1,correlated_drifts,iDrifts);                 % On enregistre les index des drifts correles
            corr_values= cat(1,corr_values,corr(iDrifts,iReg));                   % Ainsi que les valeurs
        end
    end
    clear X_norm;
end

%----SUMMARYZE INFORMATION ABOUT CORRELATION
if isempty(correlated_drifts)
    correlated_infos=[];    
else
    correlated_infos(:,1)=correlated_drifts;                                    % indices des drifts correles
    correlated_infos(:,2)=corr_values;                                          % valeur de la correlation
    correlated_infos(:,3)=(correlated_drifts-1)/(2*nTimePtsData*dtData);        % frequence correspondante
end
correlated_drifts=unique(correlated_drifts);
display(['Drifts number ' num2str(correlated_infos(:,1)') ' are correlated'])
display(['It corresponds to frequencies: ' num2str(correlated_infos(:,3)')])
display(['Corr coeff are: ' num2str(correlated_infos(:,2)')])

Dcorr=Dini(:,correlated_drifts);                                                % Dcorr=Corelated drifts

%---REMOVE  DRIFTS THAT CORRELATES
if ~isempty(correlated_drifts)
    j=1;
    for iDrifts=1:nDrifts 
        if iDrifts~=correlated_drifts(:)
            D(:,j)=Dini(:,iDrifts);
            j=j+1;
        end
    end
    
 
end


display(sprintf('%g drifts are not correlated on initialy %g drifts \n',[size(D,2) size(Dini,2)]));








% % quelques test pr les inversions
% DD=D'*D; % nb drift by nb drift
% display(sprintf('determinant de Dt*D %g',det(DD)))
% display(sprintf('rank de Dt*D %g',rank(DD)))
% display(sprintf('cond number de Dt*D %g \n',cond(DD)))
% clear DD
% 
% J=-D*inv(D'*D)*D'; % J est nb_sample by nb_sample
% for i=1:nb_ech
%     J(i,i)=J(i,i)+1;
% end
% M=X'*J*X; %nb cond by nb cond
% display(sprintf('determinant de XhJX %g',det(M)))
% display(sprintf('rank de XhJX %g',rank(M)))
% display(sprintf('cond number de XhJX %g \n\n',cond(M)))
% clear M J





