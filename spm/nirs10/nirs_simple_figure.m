function nirs_simple_figure(F,W,Z)
try
    p_value = Z.p_value;
    erdf = F.erdf;
    eidf = F.eidf;
    tstr = F.tstr;
    T_map = F.T_map;
    switch F.hb
        case 'HbO'
            nchn = length(W.ch_HbO);
        case 'HbT'
            nchn = length(W.ch_HbT);
        case 'HbR'
            nchn = length(W.ch_HbR);
    end
    %threshold
    if tstr == 'T'
        th_z = spm_invTcdf(1-p_value, erdf);
    else
        th_z = spm_invFcdf(1-p_value, eidf,erdf);
    end
    index_over = find(T_map > th_z);
    index_over2 = [find(T_map > th_z); find(-T_map > th_z)];
    I = [];
    I.index_over = index_over;
    I.index_over2 = index_over2;
    Y2 = nirs_make_simple_figure(I,F,W,Z,th_z);
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end

function Y = nirs_make_simple_figure(I,F,W,Z,th_z)
index_over = I.index_over;
index_over2 = I.index_over2;

%choose font for figure here
fontsize_choice = 16;
tick_number = 7;
th_z_shrink = 0.9; %some value slightly less than 1 required
min_max_gap = 0.0001; %some small value required -- min resolution
try
    cbar = Z.cbar;
    write_fig = 0;
    brain = W.brain;
    split = F.split;
    T_map = F.T_map;
    tstr = F.tstr;
    fcool = 64; %default, but may get overwritten later
    
    if ~isempty(index_over) ||  ~isempty(index_over2)
        if ~isempty(index_over)
            min_T = min(T_map(index_over));
            max_T = max(T_map(index_over));
            if min_T == max_T %for example if only one data point in index_over
                min_T = min_T - min_max_gap;
            end
        else
            max_T = -th_z;
            min_T = -th_z-min_max_gap;
        end
        if ~isempty(index_over2)
            min_T2 = min(T_map(index_over2));
            max_T2 = max(T_map(index_over2));
            if min_T2 == max_T2 %for example if only one data point in index_over2
                min_T2 = min_T2 - min_max_gap;
            end
        else
            max_T2 = th_z+min_max_gap;
            min_T2 = th_z;
        end
        
        if ~isempty(index_over) &&  tstr == 'F'
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
        else
            if tstr == 'T'
                %colormap resolution:
                cmap_res = 3*64;
                
                if ~isempty(index_over) && ~isempty(index_over2)
                    sbar = linspace(min_T2, max_T, cmap_res);
                    fcool = round((-th_z-min_T2)/(max_T-min_T2)*cmap_res);
                    fgray = round(2*th_z/(max_T-min_T2)*cmap_res);
                    fhot = round((max_T-th_z)/(max_T-min_T2)*cmap_res);
                    th_z_eff = th_z*th_z_shrink;
                    maxb = max(brain(:)); %minb = min(brain(:)); %expect =0
                    T_brain_over = (2*th_z_eff/maxb).*brain -th_z_eff;
                    %add both negative and positive regions of interest
                    if ~isempty(index_over)
                        T_brain_over(index_over) = T_map(index_over);
                    end
                    if  ~isempty(index_over2)
                        T_brain_over(index_over2) = T_map(index_over2);
                    end
                    if overwrite_map_thresholds
                        T_brain_over(T_brain_over>max_T) = max_T;
                        T_brain_over(T_brain_over<min_T2) = min_T2;
                        T_brain_over(1,1) = max_T;
                        T_brain_over(end,end) = min_T2;
                    end
                else
                    if ~isempty(index_over) && isempty(index_over2)
                        smin_T = max_T - ((max_T - min_T)./63) * 127;
                        sbar = linspace(smin_T, max_T, 128);
                        T_brain_over = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
                        T_brain_over(index_over) = T_map(index_over);
                        if overwrite_map_thresholds
                            T_brain_over(T_brain_over>max_T) = max_T;
                            T_brain_over(1,1) = max_T;
                            T_brain_over(end,end) = min_T;
                        end
                    else
                        if isempty(index_over) && ~isempty(index_over2)
                            fcool = 64;
                            smax_T = min_T2 + ((max_T2 - min_T2)./63) * 127;
                            sbar = linspace(min_T2,smax_T,128);
                            T_brain_over = ((-sbar(65) + sbar(128))/(0.5)).*brain + sbar(65);
                            T_brain_over(index_over2) = T_map(index_over2);
                            if overwrite_map_thresholds
                                T_brain_over(T_brain_over<min_T2) = min_T2;
                                T_brain_over(1,1) = max_T2;
                                T_brain_over(end,end) = min_T2;
                            end
                        end
                    end
                end
                
            end
        end
        write_fig = 1;
    end
    if write_fig
        if tstr == 'T'
            tmph = figure;
            jet2 = colormap(jet(2*fcool)); %doubling resolution of jet colormap
            close(tmph);
        end
        cbar.visible = 0;
        fh1 = figure('Visible',cbar.visible,'Name',F.ftiff,'NumberTitle','off');
        ax1 = axes;
        title(F.ftiff);
        imagesc(T_brain_over);
        if tstr == 'F'
            cbar1 = 1;
            cbar2 = 0;
        else
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
        end
        %Define split1
        if tstr == 'F'
            split1 = split;
        else
            if tstr == 'T'
                
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
            end
        end
        colormap(split1);
        axis(ax1, 'off')
        axis(ax1, 'image')
        %Create colorbar(s)
        if cbar1, hc1 = colorbar('EastOutside'); end
        if cbar2, hc2 = colorbar('WestOutside'); end
        %Set min and max for colorbar(s)
        if tstr == 'F'
            hc1_min = sbar(65);
            hc1_max = sbar(128);
        else
            if tstr == 'T'
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
            end
        end
        %Set colorbars limits and fontsize
        if cbar1
            if hc1_min==hc1_max %quick fix
                hc1_min=hc1_max-0.0001;
            end
            set(hc1, 'YLim', [hc1_min hc1_max]);
            y_tick = linspace(hc1_min, hc1_max, tick_number)';
            set(hc1, 'YTick', y_tick);
            set(hc1, 'FontSize', fontsize_choice);
            set(hc1,'YTickLabel',sprintf('%.1f |',get(hc1,'YTick')'));
        end
        if cbar2
            if hc2_min==hc2_max %quick fix
                hc2_min=hc2_max-0.0001;
            end
            set(hc2, 'YLim', [hc2_min hc2_max]);
            y_tick = linspace(hc2_min, hc2_max, tick_number)';
            set(hc2, 'YTick', y_tick);
            set(hc2, 'FontSize', fontsize_choice);
            set(hc2,'YTickLabel',sprintf('%.1f |',get(hc2,'YTick')'));
        end
        %Write figure
        if Z.gen_fig
            saveas(fh1,F.ffig,'fig');
        end
        if Z.gen_tiff
            print(fh1, '-dtiffn', F.ftiff);
        end
        
        %Fill Y
        Y.split1 = split1;
        Y.fh1 = fh1;
        Y.ax1 = ax1;
        if cbar1
            Y.hc1 = hc1;
            Y.hc1_min = hc1_min;
            Y.hc1_max = hc1_max;
        end
        if cbar2
            Y.hc2 = hc2;
            Y.hc2_min = hc2_min;
            Y.hc2_max = hc2_max;
        end
        Y.tick_number = tick_number;
        Y.fontsize_choice = fontsize_choice;
    else
        Y = [];
    end
    
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end