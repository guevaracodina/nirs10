function [W, Q, interpMap, Dat] = nirs_interpolation_render_precompute(Dat, NIRS, config)
%%*************************************************************************
%This function helps to inerpolate channel data and render to a selected 
%brain view (white matter layer).
%Ke Peng, LIOM, 2014-02
%Input: Dat - Channel data, Size: Channel-Number by 1, HbR or HbO or HbT
%Input: NIRS - NIRS structure after coregistration
%Input: config - configuration structure
%Output: out - DF figure structure
%%*************************************************************************
eTime = tic;

%Config
Topodatafile = NIRS.Dt.ana.rend;
brain_view = config.brain_view; %Brain view, or a loop can be used here to generate all views
AllowExtrapolation = config.AllowExtrapolation; %Option: 0: do not extrapolate
no_interpolation = config.no_interpolation; %Option: 0: interpolate
new_path = config.new_path; %Path where interpolated images are saved
figure_name = config.figure_name; %Name of interpolated image
thz = config.thz; %Threshold value of the interpolated image. Attention: Cannot be zero
%**************************************************************************
%Interpolation kernel
load(Topodatafile);
[side_hemi spec_hemi] = nirs_get_brain_view(brain_view);
rchn = rendered_MNI{brain_view}.rchn;
cchn = rendered_MNI{brain_view}.cchn;
%Two options
W.AllowExtrapolation = AllowExtrapolation;
W.no_interpolation = no_interpolation;
%find channels which are visible from this projection view
W.index_ch = find(rchn ~= -1);
%rendering surface
brain = rendered_MNI{brain_view}.ren;
W.s1 = size(brain, 1);
W.s2 = size(brain, 2);
%split into HbO and HbR interpolations
W.ch = W.index_ch;
W.rchn = rchn(W.index_ch);
W.cchn = cchn(W.index_ch);
%Generate interpolation matrix
Q = interpolation_kernel_cine_simplified(W);
Dat = Dat(W.ch,:);
%******************************************************************
%Rendering config
[W, Z, F] = render_config(Dat, W, rendered_MNI, side_hemi, spec_hemi, new_path, brain_view, thz, figure_name);
if ~exist(F.pathn, 'dir')
    mkdir(F.pathn);
end
newsavepath_fig = [F.pathn '\fig'];
if ~exist(newsavepath_fig, 'dir')
    mkdir(F.pathn,'fig');
end
F.pathnsfig = newsavepath_fig;

% Interpolated topographical map
interpMap = zeros(W.s1, W.s2);

disp(['Time to pre-compute interpolation kernel: ' datestr(datenum(0,0,0,0,0,toc(eTime)),'HH:MM:SS')]);

% Start rendering
% map_it = Dat' * Q.B;
% map_it_reshaped = zeros(W.s1,W.s2);
% map_it_reshaped(Q.index_mask) = map_it(1,:);
% F.contrast_info_both = F.pathnsfig;
% F.contrast_info_for_fig = F.contrast_info_both;
% F.s_map = map_it_reshaped;
% DF = nirs_render_map(F,W,Z);
% %Save figure
% if isfield(DF, 'fh2')
%     set(DF.fh2, 'visible', 'on');
%     DF.fh2 = save_figure(DF.fh2, Z, F);
%     pause(0.5);
%     set(DF.fh2, 'visible', 'off');
% elseif isfield(DF, 'fh1')
%     set(DF.fh1, 'visible', 'on');
%     DF.fh1 = save_figure(DF.fh1, Z, F);
%     pause(0.5);
%     set(DF.fh1, 'visible', 'off');
% end
% out = DF;

function [W, Z, F] = render_config(estiHRF, W, rendered_MNI0, side_hemi, spec_hemi, new_path, brain_view, thz, figure_name)
W.side_hemi = side_hemi;
W.spec_hemi = spec_hemi;
%View dependent info for figures
brain = rendered_MNI0{brain_view}.ren;
msk = brain>1;brain(msk)=1;
msk = brain<0;brain(msk)=0;
brain = brain * 0.5;
W.brain = brain;
%For single subject group of sessions analysis
if isfield(rendered_MNI0{W.side_hemi},'view_mask_2d')
    W.brain_view_mask_2d = rendered_MNI0{W.side_hemi}.view_mask_2d;
end
W.thz = thz;
F.pathn = new_path;
load Split
F.split = split;
F.tstr = 'T';
F.contrast_info_both = 'test4';
F.contrast_info_both_for_fig = 'test5';
Z.GroupColorbars = 0;
F.str_cor = figure_name;
F.contrast_info = figure_name;
F.contrast_info_for_fig = 'test3';
Z.gen_fig = 0;
Z.gen_tiff = 0;
map_p_idx = find(estiHRF > 0);
map_n_idx = find(estiHRF < 0);
map_c_min = min(estiHRF(map_p_idx));
map_c_max = max(estiHRF(map_p_idx));
map_c_min2 = min(estiHRF(map_n_idx));
map_c_max2 = max(estiHRF(map_n_idx));
c_min = max(map_c_min, thz);
c_max = map_c_max;
c_min2 = map_c_min2;
c_max2 = min(map_c_max2, -thz);
Z.cbar.c_min = floor(c_min*10)/10;
Z.cbar.c_max = ceil(c_max*10)/10;
Z.cbar.c_min2 = floor(c_min2*10)/10;
Z.cbar.c_max2 = ceil(c_max2*10)/10;
Z.cbar.colorbar_override = 1;
Z.cbar.visible = 'off';
Z.write_neg_pos = 0;

function fh = save_figure(fh, Z, F)
if ~Z.gen_fig %If not saved before
    filen1 = fullfile(F.pathnsfig,[F.contrast_info '.fig']);
    saveas(fh,filen1,'fig');
end
if ~Z.gen_tiff
    filen2 = fullfile(F.pathn,[F.contrast_info '.png']);
    print(fh, '-dpng', filen2,'-r300');
end
