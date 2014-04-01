function Hb = nirs_OD2Hb(inv_exs2, OD)
% Hemoglobin concentrations from optical density values
% SYNTAX
% Hb = nirs_OD2Hb(inv_exs2, OD)
% INPUTS
% inv_exs2      Inverse matrix of extinction coefficents 
%               Size #pairs*2 x #pairs*#wl, 
%               nChannels = number of pairs x number of wavelengths
% OD            Optical density column vector (nChannels x 1)
% OUTPUT
% Hb            Hemoglobin concentrations column vector (nChannels x 1)
%               First half is HbO (1:end/2), second half is HbR (end/2+1:end)
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Modified Beer-Lambert Law
% MBLL - Hb consists now of HbO and HbR, even if we had more than two
% wavelengths to begin
Hb = inv_exs2 * OD;
% Hb(:, k) = inv_exs2 * OD(:, k);

% EOF
