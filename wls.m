%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machado Alexis 03/05/2010 Ecole polytechnique de Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--INPUTS
% fs:           frequency of data
% fu:           frequency of stimuli data
% conc:         nTimePts  by nPairs
% X: design matrix, nTimePts  by nPairs
% tStimuli:     nTimePtsStim by 1
% stimuli       nTimePtsStim by nCond
% HRF



function [Betas,spectralExponents,Modelisation,Design] = wls(fs,conc,X,Opt)



%-----PARAMETERS
%-Data
nTimePts=size(X,1);

nPairs=size(conc,2);

%-Stimuli parameters
%nTimePtsStim=size(stimuli,1);
%dtStim=1/fu;
nCond=size(X,2); 

%factor=round(1/(dtStim*fs));


display ('initialisation')
%-------GET WAVELET FILTER
f = MakeONFilter('Daubechies',10);                                         % 5 vanishing moment:degree 4
display('we use a daubechie filter with 5 vanishing moment \n')

%------For decomposition
J = floor(log2(nTimePts));


L0=Opt.L0;
display(sprintf('scale max of decomposition J-L0= %g',J-L0))                    % Maximun scale of decomposition
N=2^J; 
%PP - get size larger (or smaller, is enough) than the initial data size
%N=2*N; 
% Total number of wavelet coefficients ( after decompition)
display(sprintf('nb total de coeff= %g',N))
% J-LO= effective scale of reconstruction ( low frequency are not considered)

%%%%%%%%INTERPOLATE TO MAKE DATA LENGTH A POWER OF 2 %%%%%%%%%%%%%%%%%%
conc_resized = zeros(N,nPairs);
for Idx=1:nPairs
    conc_resized(:,Idx) = resample(conc(:,Idx),N,nTimePts);
end
%need to resample X too
X_resized = zeros(N,nCond);
for Jdx=1:nCond
    X_resized(:,Jdx) = resample(X(:,Jdx),N,nTimePts);
end
X = X_resized;
%%%%%%%%%%%%

%------For nuisance model
J0=Opt.J0;
frequence_pseudo_max=(2^-J0)*fs;
n0 = (2^(-J0+1))*N ;                                                            % nb of drift coeff to estimate (J...J0)
degre_of_freedom=N-n0-nCond;
display(sprintf('minimum scale  pour les drifts J0=%g',J0))
display(sprintf('Which corresponds to a pseudo frequency %g',frequence_pseudo_max))
display(sprintf('To model nuisances there are %g coefficients to estimate',n0))
display(sprintf('Degre of freedom= %g',degre_of_freedom))





% % 
% % %--------------RESIZE DATA
% % display (sprintf('\nData are resized'))
% % conc_resized(:,:) = conc(1:N,:); %OOOOOOOOOOOOOOOPPPPPPPPPPPPPSSSSSSSSSS
% % display(sprintf('%g time samples has been removed\n',nTimePts-N))
% % clear conc
% % 
% % 



% % %-----CREATE DESIGN MATRIX
% % display ('Create design matrix')
% % for iCond=1:nCond
% %     temp_up= conv(stimuli(:,iCond),HRF.hrf);
% %     temp_up= temp_up(1:nTimePtsStim,1);
% %     % Downsample at data rate
% %     temp=downsample(temp_up,factor);
% %     clear temp_up
% %     if size(temp,1)~=nTimePts
% %         error('PROBLEM IN THE DOWNSAMPLING OF THE DESIGN MATRIX');
% %     end
% %     % resize design matrix
% %     X(:,iCond)=temp(1:N,1);
% %     clear temp
% % end

% % %--If  first and second derivatives of canonical HRF are included
% % if ~isempty(HRF.dhrf)
% %     nDer=size(HRF.dhrf,2);
% %     for ider=1:nDer
% %         for iCond= 1:nCond
% %             temp_up=conv(stimuli(:,iCond),HRF.dhrf(:,ider));
% %             temp_up= temp_up(1:nTimePtsStim,1);
% %             % dyadic..
% %             temp=downsample(temp_up,factor); %on repasse au sampling originel.
% %             clear temp_up
% %             if size(temp,1)~=nTimePts
% %                 error('PROBLEM IN THE DOWNSAMPLING OF THE DESIGN MATRIX');
% %             end
% %             X(:,ider*nCond+iCond)=temp(1:N,1);
% %             clear temp
% %         end
% %     end
% % end

nReg=size(X,2);
display(sprintf('Number of experimental conditions =%1.f',nCond))
display(sprintf('Total number of regressors (derivatives...) =%1.f',nReg))



%------CORRELATION ANALYSIS:Perform correlation between protocol and wavelet
%atoms
performCorr=0;
if performCorr
    for iReg=1:nReg
        display(sprintf('\nCorrelation analysis started'))
        tic;
        corr_stim(iReg,:) = atom_correlation(f,X(:,iReg)',L0,N);           % corr_stim: dim nCond by N
        toc
    end
else
  corr_stim=zeros(size(X,2),N);  
  display('No correlation analysis done')                                  %code pour les tests pour eviter de faire corr stim lors des tests
end


%

temp=[];
for iReg=1:nReg
    corr= find(Opt.threshold_drift < corr_stim(iReg,:));                   %stim_corr=correlation paradigm atom nb_condition by N
    temp=[temp corr];
end
correlated_drifts=unique(temp)';                                           %1 by nb_correlated drift

clear temp corr
% pour eviter d'enlever des coeffs qui existent pas
correlated_drifts=correlated_drifts(correlated_drifts <= n0);






%-----ANALYSIS
display (sprintf('\nAnalysis starts'))
for iPairs=1:nPairs
    display(sprintf('Pair number %g',iPairs))
    
    [betas,varBetas,sExp,yDrift]=weighted_estimate(X',conc_resized(:,iPairs)...
                                                    ,L0,J0,J,N,...
                                                    correlated_drifts,f);
    for iReg=1:nReg
        yHrf(:,iReg)=betas(iReg).*X(:,iReg);
    end
    yModel=yDrift+sum(yHrf,2);
    
    %--Betas
    Betas.betas(iPairs,:)=betas;
    Betas.betasVariances(iPairs,:)=varBetas;
    %--Spectral exponent
    spectralExponents(iPairs,:)=sExp;
    
    %--Reconstructed reponses
    Modelisation.yHrf{iPairs}=yHrf;
    Modelisation.yDrift{iPairs}=yDrift;
    Modelisation.yModel(:,iPairs)=yModel;
    Modelisation.residual(iPairs,:)=conc_resized(:,iPairs)-...
                                    Modelisation.yModel(:,iPairs);
    
end

Modelisation.conc=conc_resized;


%--Design
Design.J0=J0;
Design.fmax=frequence_pseudo_max;
Design.n0=n0;
Design.df=degre_of_freedom;
Design.X=X;












