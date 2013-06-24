function fh1 = nirs_make_figure_cine(I,F,W,Z,str_cor,th_z)
try
    Z.GroupColorbars = 0;
    index_over = I.index_over;
    index_over2 = I.index_over2;
    contrast_info = F.contrast_info_both;
    contrast_info_for_fig = F.contrast_info_both_for_fig;
    if isfield(Z,'Idx')
        contrast_info = ['S' gen_num_str(Z.Idx,2) '_' contrast_info];
        contrast_info_for_fig = ['S' gen_num_str(Z.Idx,2) ' ' contrast_info_for_fig];
        if isfield(Z,'subj_id')
            %This is useful if the subject numbers are not preserved due to
            %some subjects being dropped for various reasons
            contrast_info = [Z.subj_id '_' contrast_info];
            contrast_info_for_fig = [Z.subj_id ' ' contrast_info_for_fig];
        end
    end
    %choose font for figure here
    fontsize_choice = 16;
    tick_number = 7;
    th_z_shrink = 0.9; %some value slightly less than 1 required
    min_max_gap = 0.0001; %some small value required -- min resolution
    cmap_res = 3*64;
    cbar = Z.cbar;
    gen_fig = Z.gen_fig;
    gen_tiff = Z.gen_tiff;
    write_fig = 0;
    brain = W.brain;
    pathn = F.pathn;
    split = F.split;
    s_map = F.s_map;
    tstr = F.tstr;
    %fcool = 64; %default, but may get overwritten later
    %choose to overwrite max_T and min_T for map thresholds
    TC = [];
    if cbar.colorbar_override == 1
        overwrite_map_thresholds = 1;
        TC.max_T_choice = cbar.c_max;
        TC.min_T_choice = cbar.c_min;
        TC.max_T_choice2 = cbar.c_max2;
        TC.min_T_choice2 = cbar.c_min2;
    else
        overwrite_map_thresholds = 0;
    end
    fixed_colorbars = 1;
    if ~isempty(index_over) || ~isempty(index_over2)
        %Get T_brain_over
        if fixed_colorbars
            [min_T max_T min_T2 max_T2 T_brain_over fcool fgray fhot sbar] = nirs_get_Tmap_cine_fixed_colorbars(index_over,index_over2,...
                th_z,min_max_gap,s_map,overwrite_map_thresholds,brain,th_z_shrink,TC);
        else
            [min_T max_T min_T2 max_T2 T_brain_over fcool fgray fhot sbar] = nirs_get_Tmap_cine(index_over,index_over2,...
                th_z,min_max_gap,s_map,overwrite_map_thresholds,brain,th_z_shrink,TC);
        end
        write_fig = 1;
    end
    if write_fig
        tmph = figure;
        jet2 = colormap(jet(2*fcool)); %doubling resolution of jet colormap
        close(tmph);
        cbar.visible = 'off'; %c_on = 'on'; %
        fh1 = figure('Visible',cbar.visible,'Name',[str_cor ' ' contrast_info],'NumberTitle','off');
        ax1 = axes;
        title([str_cor ' ' contrast_info_for_fig]);
        if ~isreal(T_brain_over)
            disp('Warning: some values were imaginary. Only real part taken');
        end
        T_brain_over = real(T_brain_over);
        imagesc(T_brain_over);
        if ~fixed_colorbars
            if ~isempty(index_over)
                cbar1 = 1;
            else
                cbar1 = 0;
            end
            if ~isempty(index_over2)
                cbar2 = 1;
            else
                cbar2 = 0;
            end
        else
            %Overwrite
            cbar1 = 1;
            cbar2 = 1;
        end
        %Define split1
        if cbar1 && cbar2
            split1 = [jet2(1:fcool,:);gray(fgray);hot(fhot)];
        else
            if cbar1 && ~cbar2
                split1 = split;
            else
                if ~cbar1 && cbar2
                    split1 = [jet2(1:fcool,:);gray];
                end
            end
        end
        
        colormap(split1);
        axis(ax1, 'off')
        axis(ax1, 'image')
        %Create colorbar(s)
        %if cbar1
            hc1 = colorbar('EastOutside');
        %end
        %if cbar2
            hc2 = colorbar('WestOutside');
        %end
        %Set min and max for colorbar(s)
        if cbar1 && cbar2
            x1 = round((min_T-min_T2)/(max_T-min_T2)*cmap_res);
            if x1 < 1, x1 = 1; end
            hc1_min = sbar(x1);
            hc1_max = sbar(192);
            hc2_min = sbar(1);
            x1 = round((max_T2-min_T2)/(max_T-min_T2)*cmap_res);
            if x1 < 1, x1 = 1; end
            hc2_max = sbar(x1);
        else
            if cbar1 && ~cbar2
                hc1_min = sbar(65);
                hc1_max = sbar(128);
            else
                if ~cbar1 && cbar2
                    hc2_min = sbar(1);
                    hc2_max = sbar(64);
                end
            end
        end
        %Set colorbars limits and fontsize
        if cbar1
            hc1 = nirs_set_colorbar(hc1,hc1_min,hc1_max,tick_number,fontsize_choice);
        end
        if cbar2
            hc2 = nirs_set_colorbar(hc2,hc2_min,hc2_max,tick_number,fontsize_choice);
        end
%         hleg = legend(gca,F.legend); %works only for integer time_resolution
%         set(hleg,'TextColor',[1 1 1])
%         %Write figure
        if gen_fig
            pathfig = fullfile(pathn,'fig');
            if ~exist(pathfig,'dir'),mkdir(pathfig); end
            filen1 = fullfile(pathfig,[tstr '_' str_cor '_' contrast_info '.fig']);
            saveas(fh1,filen1,'fig');
        end
        if gen_tiff
            filen2 = fullfile(pathn,[tstr '_' str_cor '_' contrast_info '.png']);
            print(fh1, '-dpng', filen2,'-r300');
%             filen3 = fullfile(pathn,[tstr '_' str_cor '_' contrast_info '.pdf']);
%             print(fh1, '-dpdf', filen3);
        end
    else 
        fh1 = [];
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in make figure');
end