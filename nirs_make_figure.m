function Y = nirs_make_figure(I,F,W,Z,str_cor,th_z,combinedfig)
% leak_test1 = 1;
% if leak_test1
%     %return empty -- don't do any figures
%     Y = [];
% else
index_over = I.index_over; 
if isfield(I,'index_over2')
index_over2 = I.index_over2;
else
    index_over2 = [];
end
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
if isfield(Z,'Idx')
    contrast_info = ['S' gen_num_str(Z.Idx,2) '_' contrast_info];
    contrast_info_for_fig = ['S' gen_num_str(Z.Idx,2) ' ' contrast_info_for_fig];
    try 
    if isfield(Z,'subj_id')
        %This is useful if the subject numbers are not preserved due to
        %some subjects being dropped for various reasons
        contrast_info = [Z.subj_id '_' contrast_info]; 
        contrast_info_for_fig = [Z.subj_id ' ' contrast_info_for_fig]; 
    end
    end
end
%choose font for figure here
fontsize_choice = 16;
tick_number = 7;
th_z_shrink = 0.9; %some value slightly less than 1 required
min_max_gap = 0.0001; %some small value required -- min resolution
cmap_res = 3*64;
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
    %fcool = 64; %default, but may get overwritten later
    %choose to overwrite max_T and min_T for map thresholds
    TC = [];
    if cbar.colorbar_override == 1
        overwrite_map_thresholds = 1;
        TC.max_T_choice = cbar.c_max;
        TC.min_T_choice = cbar.c_min;
        if ~Z.GroupColorbars
            TC.max_T_choice2 = cbar.c_max2;
            TC.min_T_choice2 = cbar.c_min2;
        end
    else
        overwrite_map_thresholds = 0;
    end
    if ~isempty(index_over) || (combinedfig && ~isempty(index_over2))
        %Get T_brain_over
        [min_T max_T min_T2 max_T2 T_brain_over fcool fgray fhot sbar] = nirs_get_Tmap(F.tstr,Z.GroupColorbars,combinedfig,index_over,index_over2,...
            th_z,min_max_gap,s_map,overwrite_map_thresholds,brain,th_z_shrink,TC);
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
            hc1 = nirs_set_colorbar(hc1,hc1_min,hc1_max,tick_number,fontsize_choice);
        end
        if cbar2
            hc2 = nirs_set_colorbar(hc2,hc2_min,hc2_max,tick_number,fontsize_choice);
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
    
    nirs_save_fig_as_nifti(Z,W,pathn,combinedfig,tstr,str_cor,contrast_info,F);
    
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in make figure');
end