function interpMap = nirs_interpolation_render_compute(Q, Dat, interpMap)
% This function renders interpolates NIRS channel 
%brain view (white matter layer).
%Ke Peng, LIOM, 2014-02
%Input: Dat - Channel data, Size: Channel-Number by 1, HbR or HbO or HbT
%Input: NIRS - NIRS structure after coregistration
%Input: config - configuration structure
%Output: interpMap - DF figure structure
%%*************************************************************************

% Start rendering
map_it = Dat' * Q.B;
% Interpolated topographical map
interpMap(Q.index_mask) = map_it(1,:);

% EOF
