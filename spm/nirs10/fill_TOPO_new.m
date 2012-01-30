function TOPO = fill_TOPO_new(TOPO,C,side_hemi,f1,hb)
%Fills TOPO structure with interpolated beta, cov_beta, and F-contrast
hbi = inv_get_chromophore(hb);
%CONVENTION:
%'g': group of sessions
%'group': group analysis over subjects (see liom_group)
%'s': individual sessions
if f1 == 0
    fg = 'g';
else
    fg = 's';
end
