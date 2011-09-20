function DF = nirs_draw_figure(fign,F,W,Z)
%SORRY!!! THIS IS VERY POORLY CODED UP NOW - too many options and cases
%code is now very obscure and hard to maintain
try
    %NOTE: as the figures are generated in the Visible = Off mode, when opening
    %them in Matlab, one needs to set(gcf,'Visible','on') in order to see them
    %at all!
    % fh1 = [];
    % ax1 = [];
    % hc1 = [];
    GInv = Z.GInv;
    % if GInv
    %     %For figure combining activation and deactivation
    %     fh2 = [];
    %     ax2 = [];
    %     hc2 = [];
    % end
    p_value = Z.p_value;
    erdf = F.erdf;
    eidf = F.eidf;
    tstr = F.tstr;
    T_map = F.T_map;
    try
        sum_kappa = F.sum_kappa;
    end
    if ~(fign == 4)
        switch F.hb
            case 'HbO'
                nchn = length(W.ch_HbO);
            case 'HbT'
                nchn = length(W.ch_HbT);
            case 'HbR'
                nchn = length(W.ch_HbR);
        end
    end
    
    switch fign
        %     case 1
        %         %no threshold
        %         str_cor = 'all';
        %         th_z = 0;
        %         index_over = find(T_map);
        case 2
            %uncorrected
            str_cor = 'unc';
            %threshold
            if tstr == 'T'
                th_z = spm_invTcdf(1-p_value, erdf);
            else
                th_z = spm_invFcdf(1-p_value, eidf,erdf);
            end
            index_over = find(T_map > th_z);
            %if GInv, index_over2 = find(-T_map > th_z); end
            if GInv,
                if Z.GroupColorbars
                    index_over2 = [find(T_map > th_z); find(-T_map > th_z)];
                else
                    index_over2 = find(-T_map > th_z);
                end
            end
        case 3
            if tstr == 'T'
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
                %if GInv, index_over2 = find(-T_map > th_z); end
                if GInv,
                    if Z.GroupColorbars
                        index_over2 = [find(T_map > th_z); find(-T_map > th_z)];
                    else
                        index_over2 = find(-T_map > th_z);
                    end
                end
            else
                %For F-contrast, do a Bonferroni correction for now
                str_cor = 'Bonf';
                p_value = p_value/nchn;
                th_z = spm_invFcdf(1-p_value, eidf,erdf);
                index_over = find(T_map > th_z);
            end
        case 4
            %group
            if strcmp(F.tstr,'F')
                str_cor = 'Anova';
                xth_z = zeros(1,length(erdf));
                index_over = [];
                xerdf = erdf(:);
                udf = unique(xerdf);
                xth = zeros(1,length(udf));
                for i1=2:length(udf) %exclude 0
                    xth(i1) = spm_invFcdf(1-p_value, eidf, udf(i1));
                end
                T_map = T_map(:); %??
                for i1 = 1:length(xerdf)
                    if xerdf(i1) > 0
                        xth_z(i1) = xth(xerdf(i1)==udf);
                        if T_map(i1) > xth_z(i1)
                            index_over = [index_over i1];
                        end
                    end
                end
                th_z = min(xth_z(xth_z>0));
            else
                str_cor = 'Group';
                th_z = spm_invTcdf(1-p_value, erdf);
                index_over = find(T_map > th_z);
                %if GInv, index_over2 = find(-T_map > th_z); end
                if GInv,
                    if Z.GroupColorbars
                        index_over2 = [find(T_map > th_z); find(-T_map > th_z)];
                    else
                        index_over2 = find(-T_map > th_z);
                    end
                end
            end
    end
    I = [];
    I.index_over = index_over;
    Y1 = nirs_make_figure(I,F,W,Z,str_cor,th_z,0);
    if ~isempty(Y1)
        DF.fh1 = Y1.fh1;
        DF.ax1 = Y1.ax1;
        DF.hc1 = Y1.hc1;
        DF.split1 = Y1.split1;
        DF.hc1_min = Y1.hc1_min;
        DF.hc1_max = Y1.hc1_max;
        DF.tick_number = Y1.tick_number;
        DF.fontsize_choice = Y1.fontsize_choice;
        DF.GroupColorbars = Z.GroupColorbars;
    else
        DF = [];
    end
    
    if strfind(F.contrast_info,'Pos')
        run_other_chart = 1;
    else
        run_other_chart = 0;
    end
    if (GInv && Z.GroupColorbars)||(~Z.GroupColorbars && run_other_chart)
        if tstr == 'T'
            I.index_over2 = index_over2;
            Y2 = nirs_make_figure(I,F,W,Z,str_cor,th_z,1);
            if ~isempty(Y2)
                DF.fh2 = Y2.fh1;
                DF.ax2 = Y2.ax1;
                if Z.GroupColorbars
                    DF.hc2 = Y2.hc1;
                    DF.hc2_min = Y2.hc1_min;
                    DF.hc2_max = Y2.hc1_max;
                    DF.split2 = Y2.split1;
                else
                    if isfield(Y2,'hc1')
                        DF.hc1 = Y2.hc1;
                        DF.hc1_min = Y2.hc1_min;
                        DF.hc1_max = Y2.hc1_max;
                    end
                    if isfield(Y2,'hc2')
                        DF.hc2 = Y2.hc2;
                        DF.hc2_min = Y2.hc2_min;
                        DF.hc2_max = Y2.hc2_max;
                    end
                    DF.split1 = Y2.split1;
                end
                DF.tick_number = Y2.tick_number;
                DF.fontsize_choice = Y2.fontsize_choice;
                DF.GroupColorbars = Z.GroupColorbars;
            end
        else %'F' - careful - as this does not create a new figure, need to not close the figure
            %later on
            DF.fh2 = Y1.fh1;
            DF.ax2 = Y1.ax1;
            if isfield(Y1,'hc1')
                DF.hc2 = Y1.hc1;
                DF.hc2_min = Y1.hc1_min;
                DF.hc2_max = Y1.hc1_max;
            end
            DF.split2 = Y1.split1;
            DF.tick_number = Y1.tick_number;
            DF.fontsize_choice = Y1.fontsize_choice;
            DF.GroupColorbars = Z.GroupColorbars;
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end


function Y = nirs_make_figure(I,F,W,Z,str_cor,th_z,combinedfig)
index_over = I.index_over;
if combinedfig
    contrast_info = F.contrast_info_both;
    contrast_info_for_fig = F.contrast_info_both_for_fig;
    index_over2 = I.index_over2;
else
    contrast_info = F.contrast_info;
    contrast_info_for_fig = F.contrast_info_for_fig;
end
%choose font for figure here
fontsize_choice = 16;
tick_number = 7;
th_z_shrink = 0.9; %some value slightly less than 1 required
min_max_gap = 0.0001; %some small value required -- min resolution
try
    cbar = Z.cbar;
    gen_fig = Z.gen_fig;
    gen_tiff = Z.gen_tiff;
    write_fig = 0;
    brain = W.brain;
    pathn = F.pathn;
    split = F.split;
    T_map = F.T_map;
    tstr = F.tstr;
    fcool = 64; %default, but may get overwritten later
    %choose to overwrite max_T and min_T for map thresholds
    if cbar.colorbar_override == 1
        overwrite_map_thresholds = 1;
        max_T_choice = cbar.c_max;
        min_T_choice = cbar.c_min;
        if ~Z.GroupColorbars
            max_T_choice2 = cbar.c_max2;
            min_T_choice2 = cbar.c_min2;
        end
    else
        overwrite_map_thresholds = 0;
    end
    if ~isempty(index_over) || (combinedfig && ~isempty(index_over2))
        if Z.GroupColorbars
            if ~combinedfig
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
            else
                min_T = min(T_map(index_over2));
                max_T = max(T_map(index_over2));
                if min_T == max_T %for example if only one data point in index_over
                    if max_T > 0
                        min_T = -min_T; % - min_max_gap; %careful -- need to include grey
                    else
                        max_T = -min_T;
                    end
                end
            end
        else
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
            if combinedfig && ~isempty(index_over2)
                min_T2 = min(T_map(index_over2));
                max_T2 = max(T_map(index_over2));
                if min_T2 == max_T2 %for example if only one data point in index_over2
                    min_T2 = min_T2 - min_max_gap;
                end
            else
                max_T2 = th_z+min_max_gap;
                min_T2 = th_z;
            end
        end
        if overwrite_map_thresholds
            max_T = max_T_choice;
            min_T = min_T_choice;
            if ~Z.GroupColorbars
                max_T2 = max_T_choice2;
                min_T2 = min_T_choice2;
            end
        end
        if ~isempty(index_over) && (~combinedfig || tstr == 'F')
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
                if Z.GroupColorbars
                    sbar = linspace(min_T, max_T, cmap_res); %for colormap in 1-192 range
                    %fraction of colorbar in each of cool, gray and hot areas:
                    fcool = round((-th_z-min_T)/(max_T-min_T)*cmap_res);
                    fgray = round(2*th_z/(max_T-min_T)*cmap_res);
                    fhot = round((max_T-th_z)/(max_T-min_T)*cmap_res);
                    
                    %Choose location of background brain as a function of th_z:
                    %want it to be in between -th_z and +th_z
                    %To ensure no overlap between background brain and
                    %(de)activations, reduce th_z by a factor
                    th_z_eff = th_z*th_z_shrink;
                    maxb = max(brain(:)); %minb = min(brain(:)); %expect =0
                    T_brain_over = (2*th_z_eff/maxb).*brain -th_z_eff;
                    %add both negative and positive regions of interest
                    if ~isempty(index_over)
                        T_brain_over(index_over) = T_map(index_over);
                    end
                    if combinedfig && ~isempty(index_over2)
                        T_brain_over(index_over2) = T_map(index_over2);
                    end
                    if overwrite_map_thresholds
                        %need to rescale T_brain_over to match the specified scale
                        T_brain_over(T_brain_over>max_T) = max_T;
                        T_brain_over(T_brain_over<min_T2) = min_T2; %PROBLEM:
                        %The lower cut for the positive map, and the upper cut for
                        %the negative map, will not be implemented
                        %required for unclear reason:
                        T_brain_over(1,1) = max_T;
                        T_brain_over(end,end) = min_T2;
                    end
                    
                else
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
                        if combinedfig && ~isempty(index_over2)
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
        end
        write_fig = 1;
    end
    if write_fig
        if combinedfig && tstr == 'T'
            tmph = figure;
            jet2 = colormap(jet(2*fcool)); %doubling resolution of jet colormap
            close(tmph);
        end
        
        fh1 = figure('Visible',cbar.visible,'Name',[str_cor ' ' contrast_info],'NumberTitle','off');
        ax1 = axes;
        title([str_cor ' ' contrast_info_for_fig]);
        imagesc(T_brain_over);
        if Z.GroupColorbars || ~combinedfig || tstr == 'F'
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
        if ~combinedfig  || tstr == 'F'
            split1 = split;
        else
            if tstr == 'T'
                if Z.GroupColorbars
                    %picking lower half of jet colormap;
                    split1 = [jet2(1:fcool,:);gray(fgray);hot(fhot)];
                else
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
        end
        colormap(split1);
        axis(ax1, 'off')
        axis(ax1, 'image')
        %Create colorbar(s)
        if Z.GroupColorbars || ~combinedfig || tstr == 'F'
            hc1 = colorbar;
        else
            if cbar1
                hc1 = colorbar('EastOutside');
            end
            if cbar2
                hc2 = colorbar('WestOutside');
            end
        end
        %Set min and max for colorbar(s)
        if ~combinedfig || tstr == 'F'
            hc1_min = sbar(65);
            hc1_max = sbar(128);
        else
            if tstr == 'T'
                if Z.GroupColorbars
                    hc1_min = sbar(1);
                    hc1_max = sbar(cmap_res);
                else
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
            %Customize here number of decimals
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
            %Customize here number of decimals
            set(hc2,'YTickLabel',sprintf('%.1f |',get(hc2,'YTick')'));
        end
        %cbfreeze(hc1);
        %Write figure
        if Z.write_neg_pos || combinedfig || tstr == 'F'
            if gen_fig
                pathfig = fullfile(pathn,'fig');
                if ~exist(pathfig,'dir'),mkdir(pathfig); end
                filen1 = fullfile(pathfig,[tstr '_' str_cor '_' contrast_info '.fig']);
                saveas(fh1,filen1,'fig');
            end
            if gen_tiff
                filen2 = fullfile(pathn,[tstr '_' str_cor '_' contrast_info '.tiff']);
                if Z. SmallFigures
                    print(fh1, '-dtiff', filen2);
                else
                    print(fh1, '-dtiffn', filen2);
                end
            end
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
    %Save as nifti
    if Z.save_nifti_contrasts
        pathnii = fullfile(pathn,'nii');
        if ~exist(pathnii,'dir'),mkdir(pathnii); end
        filen3 = fullfile(pathnii,[tstr '_' str_cor '_' contrast_info '.nii']);
        %note it is the contrast that should be written, not the T or F-stat maps
        %             switch W.side_hemi
        %                 case 3 %right
        %                     M = [[0 1;-1 0] zeros(2); zeros(2) eye(2)];
        %                 case 4 %left
        M = [[0 1;-1 0] zeros(2); zeros(2) eye(2)];
        %             end
        if strcmp(tstr,'T')            
            V = nirs_create_vol(filen3,...
                [W.s1 W.s2 1], [16,0], [1;0;352],M, F.con);
        else
            V = nirs_create_vol(filen3,...
                [W.s1 W.s2 1], [16,0], [1;0;352],M, F.ess);
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
end