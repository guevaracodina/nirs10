function [K, xhat, P] = nirs_kalman_measurement_update(Pminus, R, xhatminus, z)
% Performs Kalman filter measurement update (correct)
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
% Measurement update
K       = Pminus ./ (Pminus + R);
xhat    = xhatminus + K .* (z - xhatminus);
P       = (1 - K) .* Pminus;
% K(:, k)            = Pminus(:, k) ./ (Pminus(:, k) + R);
% xhat(:, k)         = xhatminus(:, k) + K(:, k) .* (z(:, k) - xhatminus(:, k));
% P(:, k)            = (1 - K(:, k)) .* Pminus(:, k);

% EOF
