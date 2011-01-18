% =========================================================================
% FUNCTION
% create_base.m
% Create Schmidt-orthonormalized base of polynomial and sinusoidal functions
% for derive modeling.
% 
% INPUT
% (-)
% 
% OUTPUT
% base      matrix      base of othonormalized functions
% 
% COMMENT
% julien cohen-adad 04/04/2006
% =========================================================================
function D = create_base(nb_samples,t,degre,nb_freq,Dphysio)



% =========================================================================
% Polynomial drifts

% polynomial interpolation
Dp=construit_base(degre,t);



% =========================================================================
% Sinusoidal drifts

Df = spm_dctmtx(nb_samples,nb_freq)*sqrt(nb_samples);



% =========================================================================
% Physiological drifts

% frequency decomposition of data
% for i=1:size(Dphysio,2)
%     load j_filter_LP-03;
%     Dphysio_lf(:,i) = filtfilt(Hd.Numerator,1,Dphysio(:,i));
%     load j_filter_HP-09;
%     Dphysio_hf(:,i) = filtfilt(Hd.Numerator,1,Dphysio(:,i));
% end

% normalize data
for i=1:size(Dphysio,2)
    coeff = sum(Dphysio(:,i).^2)/length(Dphysio(:,i));
    Dphysio(:,i) = Dphysio(:,i)/sqrt(coeff);
end
% for i=1:size(Dphysio_lf,2)
%     coeff = sum(Dphysio_lf(:,i).^2)/length(Dphysio_lf(:,i));
%     Dphysio_lf(:,i) = Dphysio_lf(:,i)/sqrt(coeff);
% end
% for i=1:size(Dphysio_hf,2)
%     coeff = sum(Dphysio_hf(:,i).^2)/length(Dphysio_hf(:,i));
%     Dphysio_hf(:,i) = Dphysio_hf(:,i)/sqrt(coeff);
% end


% =========================================================================
% Composition of all derives

D = cat(2,Dp,Dphysio,Df(:,2:end)); % on enleve le regresseur cst de Df qui existe deja ds Dp
% Dortho = orthonormalise(D);
% clear D;
% D = Dortho;
% D = cat(2,Dphysio_lf,Dphysio_hf,Dp,Df(:,2:end));
end

% construit_base.m
% 05/01/2001
%
% Construit une base de polynomes orthonormales pour un produit scalaire
% donnee sur une grille de discretisation par orthonormalisation de Schmidt.
%
% base=construit_base(degre,grille)
%
% ENTREE
% ------
% degre  : (entier)  degre maximal des polynomes
% grille : (vecteur) points d'echantillonage
%
% SORTIE
% ------
% base   : (matrice) coordonnees des vecteurs de la base (une colonne 
%                    correspond a un polynome)

function base=construit_base(degre,grille)

% OLD
base_canonique=(grille*ones(1,degre+1)).^(ones(size(grille))*(0:degre));
base=orthonormalise(base_canonique);



% JULIEN 04/04/2006
% base_canonique=(grille*ones(1,degre+1)).^(ones(size(grille))*[0:degre]);
% base=base_canonique;

end

% orthonormalise.m
% 18/06/2002
%
% Construit une base orthonormale (B.O.N) par orthonormalisation de Schmidt.
% bon=orthonormalise(grille,base)
%
% ENTREE
% ------
%
% SORTIE
% ------

function bon=orthonormalise(base)

[l,degre]=size(base);
bon=zeros(size(base));
for j=1:degre
    norme=sum(base(:,j).^2)/l;
    prodscal=zeros(degre,1);
    for k=1:j-1
        prodscal(k)=base(:,j)'*bon(:,k)/l;
    end
    coeff=1/sqrt(norme-sum(prodscal.^2));
    bon(:,j)=base(:,j);
	for k=1:j-1
%          if (j == 2)
%              coeff = 1;
%          elseif (j == 3)
%              coeff = 1/sqrt(norme);
%          else
            bon(:,j)=bon(:,j)-prodscal(k)*bon(:,k);
%          end
    end
    bon(:,j)=bon(:,j)*coeff;
end

% [l,degre]=size(base);
% bon=zeros(size(base));
% for j=1:degre
%     norme=sum(base(:,j).^2)/l;
%     prodscal=zeros(degre,1);
%     for k=1:j-1
%         prodscal(k)=base(:,j)'*bon(:,k)/l;
%     end
%     coeff=1/sqrt(norme-sum(prodscal.^2));
%     bon(:,j)=base(:,j);
%     for k=1:j-1
%         bon(:,j)=bon(:,j)-prodscal(k)*bon(:,k);
%     end
%     bon(:,j)=bon(:,j)*coeff;
% end
end