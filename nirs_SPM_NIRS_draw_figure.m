function [fh2] = nirs_SPM_NIRS_draw_figure(th_z,brain,contrast_info,contrast_info_for_fig,T_map,flag_figure,str_cor,split)
index_over = find(T_map > th_z);
if isempty(index_over) == 1
    disp('There is no significant voxel.!!');
    T_brain_over = brain;
    fh2 = figure('Name',[contrast_info 'There is no significant voxel.!!'],'NumberTitle','off');
    title(contrast_info_for_fig);
    imagesc(T_brain_over);
    colormap(split(1:64,:));
    axis off;
    axis image;
else
    min_T = min(T_map(index_over));
    max_T = max(T_map(index_over));
    smin_T = max_T - ((max_T - min_T)./63) * 127;

    sbar = linspace(smin_T, max_T, 128);
    T_brain_over = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
    T_brain_over(index_over) = T_map(index_over);

    if flag_figure == 1                   
        fh2 = figure('Name',[contrast_info '_' str_cor],'NumberTitle','off');
        title(contrast_info_for_fig);
        imagesc(T_brain_over);
        colormap(split)
        axis off
        axis image
        hc = colorbar; %PP How to save the colorbar with the figure???
        try
        set(hc, 'YLim', [sbar(65) sbar(128)]); %PP Bug when running 'right' view contrast for Volterra - NaN seen here -- problem due to index_over not empty
        %yet there are no relevant significant voxels -- how
        %could that be? (in the problematic case, there was
        %only one such voxel)
        y_tick = linspace(sbar(65), sbar(128), 6)';
        set(hc,'yticklabel',sprintf('%.2f |',get(hc,'ytick')'));
        set(hc, 'YTick', y_tick);
        set(hc, 'FontSize', 8);
        catch
        end
    end
end
%fh2 = imcapture(gcf, 'all', 150);
end