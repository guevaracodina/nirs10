function OD = nirs_I2OD(I, I_0, EPF2)
% Optical Density computation form measured intensity and mean intensity.
% SYNTAX
% OD = nirs_I2OD(I, I_0, EPF2)
% INPUTS
% I         NIRS measured intensity
% I_0       NIRS mean signal intensity
% EPF2      Effective PathLength Factor
% OUTPUT
% OD        Optical Density
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Optical density changes are given by:
OD = real(log10(I ./ I_0));
% OD(:, k) = real(log10(xhat(:, k) ./ xhat_DC(:, k)));
% Multiply by 1e6 to get micromolar units negative sign so that an increase
% in chromophore concentration corresponds to a decrease in intensity due to
% light absorption
OD = -1e6 * OD ./ EPF2;
% OD(:, k) = -1e6 * OD(:, k) ./ EPF2(k, :)'; %PP used to be d = d ./ EPF2;

% EOF
