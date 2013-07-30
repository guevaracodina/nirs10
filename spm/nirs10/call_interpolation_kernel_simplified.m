NC = 32; %Number of HbO channels
load('TopoData.mat'); %this loads the structure rendered_MNI

% VIEWS
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'
brain_view = 2; 
rchn = rendered_MNI{brain_view}.rchn;
cchn = rendered_MNI{brain_view}.cchn;
%Two options
W.AllowExtrapolation = 0; %Option: 0: do not extrapolate
W.no_interpolation = 0; %Option: 1: do not interpolate

%find channels which are visible from this projection view
W.index_ch = find(rchn ~= -1);
%rendering surface
brain = rendered_MNI{brain_view}.ren;
W.s1 = size(brain, 1);
W.s2 = size(brain, 2);
%split into HbO and HbR interpolations
W.ch_HbO = W.index_ch;
W.rchn = rchn(W.index_ch);
W.cchn = cchn(W.index_ch);

Q = interpolation_kernel_cine_simplified(W);
