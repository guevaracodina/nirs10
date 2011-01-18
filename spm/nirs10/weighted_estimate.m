%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machado Alexis 03/05/2010 Ecole polytechnique de Montreal from code of
% Carl Matteau Pelletier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weighted_estimate() - Estimate the HRF strength using the WLS algorithm

% Inputs:
%   response : Design matrix Attention:nCond by nTimePts
%   y        : Data nTimePts by 1
%   L0       : J-L0 effective degree of decompostion
%   J0       : Scale minimum for nuisance representation
%   J        :  Possible Scale maximum of decomposition
%   N        : Number of wavelet coefficients
%   correlated_drifts: coefficients to set to zero 
%
% Outputs:
%   beta               : Estimated beta  1 by nCond
%   spectral_exponent  : Estimated spectral exponent scalar
%   varBetas           : Estimated  variance of beta 1 by nCond
%   yDrift             : Signal of nuisance nTimePts by 1



function [betas, varBetas, spectral_exponent, yDrift] = ...
    weighted_estimate(response,y,L0,J0,J,N,correlated_drifts,f)


%--Parameters
nCond=size(response,1);
n0 = 2^(-J0+1)*N;

%--------------------OLS pour la 1re iteration--------------------------
display('iteration 1: OLS')
%--Wavelet Transform of paradigm
for iCond=1:nCond
Wx(iCond,:)=FWT_PO(response(iCond,:),L0,f);                                     % Wx: nCond by N ( nb wavelet coeff )
end
Wx=Wx';                                                                         % N by nCond
%-- Wavelet Transform of raw data
Wy=FWT_PO(y,L0,f);                                                              % Wy: N by 1 
   
%--Make A with multi conditions
%--First make identity matrix part of A in sparse representation
index_i = 1:n0;                                                                 % [1 2...n0]
index_j = 1:n0;                                                                 % [1 2...n0]
value_ij = ones(1,n0);                                                          % [1...1] n0 times
			
%--Last column of A
for iCond=1:nCond
temp = 1:N;                                                                     %[1 2...N]
index_j = [index_j temp];                                                       % 1 by N+n0 [1 2..n0 1 2..N 1 2...N......1 2...N]
end

for iCond=1:nCond
temp = (n0+iCond)*ones(1,N);
index_i=[index_i temp];                                                         % 1 by N+n0 [1 2..n0 no+1...no+1......no+5...no+5]
value_ij = [value_ij Wx(:,iCond)'];                                             % 1 by N+n0 [1...1 w1...wn............w1...wn]			
end

A=sparse(index_j,index_i,value_ij,N,n0+nCond);                                  % A: N by n0+nCond

%--We remove colunms of A correpsonding to correlated drift coeffcients
    if ~isempty(correlated_drifts)
        A(:,correlated_drifts(:)) = [];                                         % Matrix B in Matteau Article correlated coeff are already removed
    end			


%------ESTIMATION OLS
B=A'*A;
C=A'*Wy;
XI=B\C;                                                                         % estimation de XI tild (15) equivalent to inv(B)*C inverse par elimination gausienne
                                                                                % XI=[n0 coeff for nuisances... betasCond1 betasCond2 ....]



% zeros padding for XI ( correspond to coefficient that have been not considered or removed by correlation analysis)
XItemp=zeros(n0+nCond,1); 				 
if ~isempty(correlated_drifts)
    j=1;
    for k=1:n0+nCond
        temp = find(k==correlated_drifts);
        iselement = nnz(temp);                                                  % nb of non zero elements
            if ~iselement
                XItemp(k,1)=XI(j,1);
                j=j+1;
            else
				XItemp(k,1)=0;
            end
    end
    clear XI;
    XI=XItemp;
    clear XItemp;
end
				
% Coefficients of Nuisances
temp=XI(1:(end-nCond));
temp2=zeros(1,N);
temp2(1:n0) = temp;                                                             %[1..n0coeff 0 0 0 0 ... 0] 1 by N				
betas(1,:) = XI(end-(nCond-1):end);                                             % betas ( 1 by nCond) [betas1 betas2 betas 3.....]

residual = Wy - temp2' - Wx*betas';                                             % matrice des coeff DWT du bruit









%--------------------Iteration 2:Weigthed Least square
display('iteration 2 WLS')
%Generate L matrix
L = adapt_wavedec(residual,J);                                                  % 1 by scale_max+2 ici 1 by J+2 utilisé in loglike.m et detcoef.m (voir matlab help)

%Initial value of cGamma
cGamma = 0.5;
spectral_exponent = 0;
Sigma2= 1.0;
cutoff = 1e-7;                                                                  % epsilon is little and limits tends to zero

nbr_level = J; 
var_approx = SaJ(cGamma,spectral_exponent,J,cutoff);                            % RaJ(gamma) in article
for index_j = 1:nbr_level
    var_detail(index_j)=SdJ(cGamma,spectral_exponent, index_j);                 % Rdj(gamma)
end

residual_cell = num2cell(residual);
c=detcoef(residual_cell,L,'cells');                                              % retourne dans la cell c les detail coeff de chaque level (dkj du residu)

%ML (total) residual variance
Sigma2 = sigma_ml(residual(1), c, var_approx, var_detail,nbr_level, J);         % equation 21: estimation de sigma carre

approx_coeff = residual(1); % coeff d'approx scale J
% on cherche le spectral exponent qui minime la fonction de vraisemblance
[spectral_exponent,fval,exitflag] = fminsearch(@(spectral_exponent) ...
                                    loglike(cGamma,spectral_exponent,N,J,...
                                    nbr_level,L,approx_coeff,c,cutoff),0.5);


% reconstruction de la matrice sigma majuscule
value_ij =  SaJ(cGamma,spectral_exponent, J,cutoff);


for index_j = 1:nbr_level
    new_SdJ(index_j)=SdJ(cGamma, spectral_exponent, index_j);
end
new_SdJ = new_SdJ(end:-1:1);

for index_j = 1:J
    value_ij = [value_ij new_SdJ(index_j)*ones(1,L(index_j+1)) ];
end

index_line = 1:N;
inverse_sigma_matrix =sparse(index_line,index_line,1./value_ij,N,N);            % matrice inverse SIGMA


% Estimation de xi
%display('iteration 1 estimation')
B=A'*inverse_sigma_matrix*A;
C=A'*inverse_sigma_matrix*Wy;
XI=B\C;

XItemp=zeros(n0+nCond,1);					 
if ~isempty(correlated_drifts)
    j=1;
    for k=1:n0+nCond
        temp = find(k==correlated_drifts);
        iselement = nnz(temp);
            if ~iselement
                XItemp(k,1)=XI(j,1);
                j=j+1;
            else
				XItemp(k,1)=0;
            end
    end
    clear XI;
    XI=XItemp;
    clear XItemp;
end
	

% on garde les theta
temp=XI(1:end-(nCond));
temp2=zeros(1,N);
temp2(1:n0) = temp;		
% yDrift = IWT_PO(temp2,L0,f); % reconstruction des drifts
betas(1,:) = XI(end-(nCond-1):end);                                          % betas ( 1 by nCond) [betas1 betas2 betas 3.....]









%%%%%%%%%%%%%%%%%%%%   Iterations supplementaires  %%%%%%%%%%%%%%%%%%%%%%%%


n_iteration = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iteration = 1:n_iteration
display(sprintf('iteration %1.f',2+iteration))

%Residual
residual = Wy - temp2' -Wx*betas';

%Generate L matrix
L = adapt_wavedec(residual,J);

nbr_level = J; % number of level for detail coefficients
var_approx = SaJ(cGamma,spectral_exponent,J,cutoff);

for index_j = 1:nbr_level
    var_detail(index_j)=SdJ(cGamma,spectral_exponent, index_j);
end

residual_cell = num2cell(residual);
c=detcoef(residual_cell,L,'cells');

%ML (total) residual variance
Sigma2 = sigma_ml(residual(1), c, var_approx, var_detail,nbr_level, J);

approx_coeff = residual(1);

[spectral_exponent,fval,exitflag] = fminsearch(@(spectral_exponent) ...
                                    loglike(cGamma,spectral_exponent,N,J,...
                                    nbr_level,L,approx_coeff,c,cutoff),[0]);

% spectral_exponent

% reconstruction de la matrice sigma majuscule
value_ij =  SaJ(cGamma,spectral_exponent, J,cutoff);

for index_j = 1:nbr_level
    new_SdJ(index_j)=SdJ(cGamma, spectral_exponent, index_j);
end

new_SdJ = new_SdJ(end:-1:1);

for index_j = 1:J
    value_ij = [value_ij new_SdJ(index_j)*ones(1,L(index_j+1)) ];
end

index_line = 1:N;
inverse_sigma_matrix =sparse(index_line,index_line,1./value_ij,N,N);

%Estimatin xi
%display(sprintf('estimation iteration %1.f',iteration+1))
B=A'*inverse_sigma_matrix*A;
C=A'*inverse_sigma_matrix*Wy;
XI=B\C;

					 
XItemp=zeros(n0+nCond,1);					 
if ~isempty(correlated_drifts)
    j=1;
    for k=1:n0+nCond
        temp = find(k==correlated_drifts);
        iselement = nnz(temp);
            if ~iselement
                XItemp(k,1)=XI(j,1);
                j=j+1;
            else
				XItemp(k,1)=0;
            end
    end
    clear XI;
    XI=XItemp;
    clear XItemp;
end

%reconstruction des drift
temp=XI(1:end-(nCond));%n0 coeff
temp2=zeros(1,N);
temp2(1:n0) = temp;	%[1..n0coeff 0 0 0 0 0] 1 by N			
%xi = IWT_PO(temp2,L0,f);

betas(1,:) = XI(end-(nCond-1):end);

end
display ('Estimation finished')





% Signal of nuisance
yDrift=IWT_PO(temp2,L0,f); 
yDrift=yDrift';                                                             % nTimePts by 1


% Calcul des variances
residual = Wy - temp2' -Wx*betas';
sigma2_WLS = Sigma2;                                                        % scalar
var_xi  = inv(B)*sigma2_WLS;                                                % B=A'*inverse_sigma_matrix*A;
temp=diag(var_xi);
varBetas(1,:)=temp(end-(nCond-1):end);
clear temp






% % calcul des tstats OLD version which is good to see how use contrasts
% size_design_matrix = size(var_xi,1);
% for iCond=1:nCond
% contrast = zeros(size_design_matrix,1);
% contrast(end-(nCond-1)+(iCond-1)) = 1;
% varBetas(1,iCond)=(contrast'* var_xi*contrast);
% t_wls(1,iCond) = (betas(iCond)-0)/sqrt(varBetas(1,iCond));
% end













