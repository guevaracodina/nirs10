function [fh1 ax1 hc] = nirs_draw_figure(fign,brain,T_map,contrast_info,...
    contrast_info_for_fig,split,pathn,erdf,sum_kappa,p_value,gen_fig,gen_tiff,cbar)
%NOTE: as the figures are generated in the Visible = Off mode, when opening
%them in Matlab, one needs to set(gcf,'Visible','on') in order to see them
%at all!
fh1 = [];
ax1 = [];
hc = [];

write_fig = 0;
%choose font for figure here
fontsize_choice = 16;
%choose to overwrite max_T and min_T for map thresholds 
if cbar.colorbar_override == 1
    overwrite_map_thresholds = 1;
    max_T_choice = cbar.c_max;
    min_T_choice = cbar.c_min;
else
    overwrite_map_thresholds = 0;
end

switch fign
    case 1
        %no threshold
        str_cor = 'all';
        th_z = 0;
        index_over = find(T_map);
    case 2
        %uncorrected
        str_cor = 'unc';
        %threshold
        th_z = spm_invTcdf(1-p_value, erdf);
        index_over = find(T_map > th_z);
    case 3
        %Tube
        str_cor = 'tube';
        %find threshold
        z_value = 1:0.0001:7;
        p_value_tube = ((sum_kappa * gamma(3/2))/(2*(pi^(3/2))))*...
            (1-gammainc((z_value(:).^2)/2, 3/2));
        index_z = [];
        ini_ran = 10^(-10);
        n = 0;
        while isempty(index_z) == 1
            ran = ini_ran * (10^n);
            n = n+1;
            index_z = find(p_value_tube > p_value - ran & p_value_tube < p_value + ran);
        end
        index_z = index_z(end);
        th_z = z_value(index_z);  
        index_over = find(T_map > th_z);
    case 4
        %group 
        str_cor = 'Group';
        th_z = spm_invTcdf(1-p_value, erdf);
        index_over = find(T_map > th_z);
end

if ~isempty(index_over) 
    min_T = min(T_map(index_over));
    max_T = max(T_map(index_over));
    if min_T == max_T %for example if only one data point in index_over
        min_T = min_T - 0.0001;
    end
    if overwrite_map_thresholds
        max_T = max_T_choice;
        min_T = min_T_choice;
    end
    smin_T = max_T - ((max_T - min_T)./63) * 127;
    sbar = linspace(smin_T, max_T, 128);
    T_brain_over = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
    T_brain_over(index_over) = T_map(index_over);
    if overwrite_map_thresholds
        %need to rescale T_brain_over to match the specified scale
        T_brain_over(T_brain_over>max_T) = max_T;
        T_brain_over(1,1) = max_T;
        T_brain_over(end,end) = min_T;
    end
    write_fig = 1;
end
if write_fig
    fh1 = figure('Visible',cbar.visible,'Name',[str_cor '_' contrast_info],'NumberTitle','off');  
    ax1 = axes;
    title([str_cor ' ' contrast_info_for_fig]);
    imagesc(T_brain_over); 
    colormap(split)
    axis(ax1, 'off')
    axis(ax1, 'image')
    hc = colorbar;
    set(hc, 'YLim', [sbar(65) sbar(128)]);    
    y_tick = linspace(sbar(65), sbar(128), 7)';
    set(hc, 'YTick', y_tick);
    set(hc, 'FontSize', fontsize_choice);
    %Customize here number of decimals
    set(hc,'YTickLabel',sprintf('%.1f |',get(hc,'YTick')'));   
    filen1 = fullfile(pathn,[str_cor '_' contrast_info '.fig']);
    if gen_fig
        saveas(fh1,filen1,'fig');
    end
    if gen_tiff
        filen2 = fullfile(pathn,['T_ ' str_cor '_' contrast_info '.tiff']);
        %lower resolution - compressed
        %print(fh1, '-dtiff', filen2);
        %full resolution, 10x larger, no compression:
        print(fh1, '-dtiffn', filen2);
%         try movegui(fh1); end
%         [cdata1,~] = getframe(fh1);    
%         filen2 = fullfile(pathn,['T_ ' str_cor '_' contrast_info '.tiff']);
%         imwrite(cdata1,filen2);
    end
    %try close(fh1); end
end
end