function Z = get_contrast_group_common_options(job)
%To generate inverted hemodynamic responses. Be careful when using this option
Z.GInv = job.display_options.GenerateInverted;
%Display feature - to group figures into subplots - careful, this may not
%group subplots in the correct order depending on which contrasts and other
%options you have selected
Z.GFIS = job.display_options.GroupFiguresIntoSubplots;
%Other options
Z.GroupColorbars = job.display_options.GroupColorbars;
Z.SmallFigures = job.display_options.SmallFigures;
Z.write_neg_pos = job.display_options.write_neg_pos;
Z.save_nifti_contrasts = job.display_options.save_nifti_contrasts;
Z.output_unc = job.display_options.output_unc;
Z.p_value = job.contrast_p_value;
Z.figures_visible = job.display_options.figures_visible;
%Colorbar - to allow specifying common min and max on all charts
if isfield(job.display_options.override_colorbar,'colorbar_override')
    Z.cbar.c_min = job.display_options.override_colorbar.colorbar_override.colorbar_min;
    Z.cbar.c_max = job.display_options.override_colorbar.colorbar_override.colorbar_max;
    Z.cbar.c_min2 = job.display_options.override_colorbar.colorbar_override.colorbar_min2;
    Z.cbar.c_max2 = job.display_options.override_colorbar.colorbar_override.colorbar_max2;
    Z.cbar.colorbar_override = 1;
else
    Z.cbar.colorbar_override = 0;
end
switch Z.figures_visible
    case 1
        Z.cbar.visible = 'on';
    case 0
        Z.cbar.visible = 'off';
end
%Booleans to choose which figures to write to disk, if any
switch job.display_options.contrast_figures
    case 0
        Z.gen_fig = 0;
        Z.gen_tiff = 0;
    case 1
        Z.gen_fig = 1;
        Z.gen_tiff = 1;
    case 2
        Z.gen_fig = 1;
        Z.gen_tiff = 0;
    case 3
        Z.gen_fig = 0;
        Z.gen_tiff = 1;
end