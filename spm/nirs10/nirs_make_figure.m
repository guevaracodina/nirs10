function Y = nirs_make_figure(I,F,W,Z,str_cor,th_z,combinedfig)
index_over = I.index_over;
if combinedfig
    contrast_info = F.contrast_info_both;
    contrast_info_for_fig = F.contrast_info_both_for_fig;
    index_over2 = I.index_over2;
else
    if strcmp(F.tstr,'T')
        contrast_info = F.contrast_info;
        contrast_info_for_fig = F.contrast_info_for_fig;
    else
        contrast_info = F.contrast_info_both;
        contrast_info_for_fig = F.contrast_info_both_for_fig;
    end
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
    s_map = F.s_map;
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
                    min_T = min(s_map(index_over));
                    max_T = max(s_map(index_over));
                    if min_T == max_T %for example if only one data point in index_over
                        min_T = min_T - min_max_gap;
                    end
                else
                    max_T = -th_z;
                    min_T = -th_z-min_max_gap;
                end
            else
                min_T = min(s_map(index_over2));
                max_T = max(s_map(index_over2));
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
                min_T = min(s_map(index_over));
                max_T = max(s_map(index_over));
                if min_T == max_T %for example if only one data point in index_over
                    min_T = min_T - min_max_gap;
                end
            else
                max_T = -th_z;
                min_T = -th_z-min_max_gap;
            end
            if combinedfig && ~isempty(index_over2)
                min_T2 = min(s_map(index_over2));
                max_T2 = max(s_map(index_over2));
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
            T_brain_over(index_over) = s_map(index_over);
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
                        T_brain_over(index_over) = s_map(index_over);
                    end
                    if combinedfig && ~isempty(index_over2)
                        T_brain_over(index_over2) = s_map(index_over2);
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
                            T_brain_over(index_over) = s_map(index_over);
                        end
                        if combinedfig && ~isempty(index_over2)
                            T_brain_over(index_over2) = s_map(index_over2);
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
                            T_brain_over(index_over) = s_map(index_over);
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
                                T_brain_over(index_over2) = s_map(index_over2);
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
        if ~isreal(T_brain_over)
            disp('Warning: some values were imaginary. Only real part taken');
        end
        T_brain_over = real(T_brain_over);
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
                if Z.SmallFigures
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
    if Z.save_nifti_contrasts && ( Z.write_neg_pos || combinedfig || tstr == 'F' )
        pathnii = fullfile(pathn,'nii');
        if ~exist(pathnii,'dir'),mkdir(pathnii); end
        %NP = not permuted
        filen5 = fullfile(pathnii,[tstr '_' str_cor '_' contrast_info 'NP.nii']);
        filen3 = fullfile(pathnii,[tstr '_' str_cor '_' contrast_info '.nii']);
        %note it is the contrast that should be written, not the T or F-stat maps
        M = [[0 1;-1 0] zeros(2); zeros(2) eye(2)];
        if strcmp(tstr,'T')
            V = nirs_create_vol(filen5,...
                [W.s1 W.s2 1], [16,0], [1;0;352],M, F.con);
        else
            V = nirs_create_vol(filen5,...
                [W.s1 W.s2 1], [16,0], [1;0;352],M, F.ess);
        end
        %             end
        % test : on se trompe pas pour le milieu
        % 1:ventral, 2:dorsal, 3:right, 4:left, 5:frontal, 6:occipital
        switch W.side_hemi
            case 2% s2 x vx and s1 z vx
                dim = [W.s2 1 W.s1];
                M = [[0 0 1;1 0 0;0 1 0;0 0 0] zeros(4,1)];
                vecta = -M(1:3,1:3)*[round(dim(1)/2);1;round(dim(3)/2)];
                Fcon = permute(F.con,[2,3,1]);
                Fess = permute(F.ess,[2,3,1]);
                
            case 3% s2 -x vx and s1 y vx
                dim = [W.s2 W.s1 1];
                M = [[0 0 1;1 0 0;0 1 0;0 0 0] zeros(4,1)];
                vecta = -M(1:3,1:3)*[round(dim(1)/2);round(dim(2)/2);1];
                Fcon = permute(F.con,[2,1,3]);
                Fess = permute(F.ess,[2,1,3]);
            case 4% s2 x vx and s1 y vx
                dim = [W.s2 W.s1 1];
                M = [[0 0 1;-1 0 0;0 1 0;0 0 0] zeros(4,1)];%%% matrice de base
                vecta = -M(1:3,1:3)*[round(dim(1)/2);round(dim(2)/2);1];
                Fcon = permute(F.con,[2,1,3]);
                Fess = permute(F.ess,[2,1,3]);
                
            case 5% s2 -z vx and s1 y vx
                dim = [1 W.s1 W.s2];
                M = [[0 0 1;-1 0 0;0 -1 0;0 0 0] zeros(4,1)];
                vecta = -M(1:3,1:3)*[1;round(dim(1)/2);round(dim(2)/2)];
                Fcon = permute(F.con,[3,1,2]);
                Fess = permute(F.ess,[3,1,2]);
                %case 6% s2 z vx and s1 y vx
            otherwise % cas sans interet
        end
        M(:,4) = [vecta;1];
        % test : use the V.mat of the T1
        %         path_T1 = [fileparts(fileparts(pathn)) '\T1'];
        %         V = spm_vol([path_T1 '\T1.nii']);
        %         M(1:3,1:3) = V.mat(1:3,1:3)./norm(V.mat(1:3,1:3));
        if strcmp(tstr,'T')
            V = nirs_create_vol(filen3,...
                dim, [16,0], [1;0;352],M, Fcon);
            %                 [W.s1 W.s2 1], [16,0], [1;0;352],M, F.con);
            
        else
            V = nirs_create_vol(filen3,...
                dim, [16,0], [1;0;352],M, Fess);
            %                 [W.s1 W.s2 1], [16,0], [1;0;352],M, F.ess);
            
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in make figure');
end