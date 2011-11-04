function out = nirs_run_reconstruct2(job)
% Clément Bonnéry October 2011
% CONCLUSION : 
% --la PCA n explique rien, les coefficients montre que chaque
% vecteur a son importance... il faudrait peut etre passer en referentiel
% polaire!
% --la SVD crashe par manque de memoire...

% Reconstruction TEMPORELLE
% / dOD(channel,t0) = K dmua(vx) + e1
%| svd ou PCA(T1seg) soit dmua(vx) = W dmua_MG + e2
%| dmua_MG(vx) = [ext]*[DHbO(vx);DHbR(vx)] + e3
% \ [DHbO(vx);DHbR(vx)] = 0 + omega
%
% Reconstruction SPATIO-TEMPORELLE
% / dOD(channel,t) = K dmua(vx,t) + e1
%| svd ou PCA(T1seg) soit dmua(vx,t) = W dmua_MG(vx,t) + e2
%| dmua_MG(vx) = D(vx)B(vx) + e3 (voir comment reduire la taille temporellement)
% \ B = B0 + omega
%
%
% e1 : erreur entre pb direct et reqlite
% e2 : erreur commise en reduisant la taille des donnes et nuisances


% svd ou pca sur dOD, finalement, c est ce qui nous interesse
load('D:\Users\Clément\Etude2_Said_GLM\S004\fir\b04_stroop85.nirs','-mat');
[COEFF,SCORE,latent,tsquare] = princomp(d);

    % svd ou PCA sur les voxels :
    cs.segR = 'D:\Users\Clément\DPF_testDuncan3\S007\roi_all-channels\roi_00022_segmented_s201009151630-0003-00001-000192-01.nii';
    V = spm_vol(cs.segR);
    Y = spm_read_vols(V);
    % on passe en 2D pour pouvoir effectuer les calculs... est ce qu on ne
    % perd pas de l info la ?
    Y2 = Y(:);
    Y3 = zeros(numel(Y),4);
    count=1;
    % pour pouvoir faire les choses proprement, il faut un vecteur
    % Y3 =(x; y; z; layer)
    for i=1:size(Y,1)
        for j=1:size(Y,2)
            for k=1:size(Y,3)
                Y3(count,:) = [i,j,k,Y(i,j,k)];
                count = count+1;
            end
        end
    end
    %
    % SVD
    [U,S,V] = svd(Y3);
    % PCA
    [COEFF,SCORE,latent,tsquare] = princomp(Y3);

% sur le temps :


%
% ce que l on peut ecrire sous la forme
% /Y(t0)\   /Xsens*[ext] Xsens*[ext]\  /omega\   /epsilon_channelnoise\
% |  0  | = |      1           0    |* \beta0/ + |       -omega        |
% \  0  /   \     0           1     /            \       -beta0        /
%


end