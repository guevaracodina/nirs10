
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machado Alexis 03/05/2010 Ecole polytechnique de Montreal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% makePca.m
% Description: Fonction qui fait le fitrage des donnees NIRS avec une PCA

%--INPUT
% d: data ; dim: timepts by nb channels
% channels: channels to consider for pca; dim vector
% nsv: number of singular values to remove; dim scalar

%--OUTPUT
% d_PCA: filtered data; dim: timepts by nb channels







function [d_PCA] = makePca(d,channels,nSV)


display(sprintf('les donnees sont compose de %g channels',size(d,2)));
display(sprintf('la pca est faite sur les channels %g', channels));



%-----------------PREPARATION DES DONNEES---------------------------------
%--On normalise toute les channels pour quelles soit comparables en terme d'intensite, 
d_Norm=d./(ones(timepts,1)*mean(d));                                            % ici moyenne a 1
d_Norm_zM=d_Norm -1;                                                            % substract the mean to have zero mean

%------------------SINGULAR VALUE DECOMPOSITION--------------------------------
%  covarience matrix:c=cov(Intens_Norm_zM(:,channels)); %cov(X1,X2)=(1/n-1)*sum_i((xi-mux)(yi-muy)) avec mux et muy=0
c = d_Norm_zM(:,channels)'* d_Norm_zM(:,channels);

%single value decomposition  X = U*S*V' (reels V'=V*(adjoint) #voir conjugate transpose sur wiki)
[dummy,S,V] = svd(c);
                                                                                % U: matrice carre nb_channels by nb_channels  
                                                                                % S:nb_channels by nb_channels 
                                                                                % V:nb_channels by nb_channels 
                                                                                % On a U=V car c est une matrice normale donc SVS equivalente a X = U*S*U'  avec U la matrice des vecteurs propres et S contenant les valeurs propres
 
                                                                                
                                                                                
%-----------------------------DISPLAY SV-SPECTRUM ----------                                                                         
% svs = diag(S);                                                                  
% nb_sv=size(svs,1);
% relativ_var=svs/sum(svs);
% cumsvs=cumsum(relativ_var);

% figure
% subplot(1,2,1)
% plot([1:nb_sv],relativ_var,'b.')
% title('relative variance')
% subplot(1,2,2)
% plot([1:nb_sv],cumsvs,'g.')
% title('cumulative variance')





%--------------------------FILTERING---------------------------------------
W = d_Norm_zM(:,channels)*V/S;                                             % W: timepts by nb_channels
lstSV = 1:nSV;                                                                  % liste des singular value a enlever
d_PCA = d_Norm;                                                                 % A verifier
%on ne soustrait que les channels choisies
d_PCA(:,channels) = d_Norm_zM(:,channels) - W(:,lstSV)*S(lstSV,lstSV)*V(:,lstSV)';
%add one to undo the effect of zero-mean
d_PCA=d_PCA+1;  



