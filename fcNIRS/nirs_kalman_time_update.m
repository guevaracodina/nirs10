function [xhatminus, Pminus] = nirs_kalman_time_update(xhat, P, Q)
% Performs Kalman filter time update (predict)
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
% Time update
xhatminus   = xhat;
Pminus      = P + Q;
% xhatminus(:, k)    = xhat(:, k-1);
% Pminus(:, k)       = P(:, k-1) + Q;

% EOF
