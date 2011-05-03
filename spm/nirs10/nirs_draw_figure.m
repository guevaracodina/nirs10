function DF = nirs_draw_figure(fign,F,W,Z)
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
    case 1
        %no threshold
        str_cor = 'all';
        th_z = 0;
        index_over = find(T_map);
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
        if GInv, index_over2 = [find(T_map > th_z); find(-T_map > th_z)]; end
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
            if GInv, index_over2 = [find(T_map > th_z); find(-T_map > th_z)]; end
        else
            %For F-contrast, do a Bonferroni correction for now
            str_cor = 'Bonf';
            p_value = p_value/nchn;
            th_z = spm_invFcdf(1-p_value, eidf,erdf);
            index_over = find(T_map > th_z);
        end
    case 4
        %group 
        str_cor = 'Group';
        th_z = spm_invTcdf(1-p_value, erdf);
        index_over = find(T_map > th_z);
        %if GInv, index_over2 = find(-T_map > th_z); end
        if GInv, index_over2 = [find(T_map > th_z); find(-T_map > th_z)]; end
end
[fh1 ax1 hc1 split1 hc1_min hc1_max tick_number fontsize_choice] = nirs_make_figure(index_over,F,W,Z,str_cor,th_z,0);

DF.fh1 = fh1;
DF.ax1 = ax1;
DF.hc1 = hc1;
DF.split1 = split1;
DF.hc1_min = hc1_min;
DF.hc1_max = hc1_max;
DF.tick_number = tick_number;
DF.fontsize_choice = fontsize_choice;
if GInv 
    if tstr == 'T'
        [fh2 ax2 hc2 split2 hc2_min hc2_max tick_number fontsize_choice] = nirs_make_figure(index_over2,F,W,Z,str_cor,th_z,1);
        DF.fh2 = fh2;
        DF.ax2 = ax2;
        DF.hc2 = hc2;
        DF.split2 = split2;
        DF.hc2_min = hc2_min;
        DF.hc2_max = hc2_max;
    else %'F' - careful - as this does not create a new figure, need to not close the figure
        %later on
        DF.fh2 = fh1;
        DF.ax2 = ax1;
        DF.hc2 = hc1;
        DF.split2 = split1;
        DF.hc2_min = hc1_min;
        DF.hc2_max = hc1_max;
    end
end
catch exception
    disp(exception.identifier);
end 
end


function [fh1 ax1 hc1 split1 hc1_min hc1_max tick_number fontsize_choice] = nirs_make_figure(index_over,F,W,Z,str_cor,th_z,combinedfig)
%fh1 = []; ax1 = []; hc1 = []; split1 = []; 
%hc1_min = 0; hc1_max = 1; %Not a good or robust solution to fh1 not being generated when there are no figures...
try 
cbar = Z.cbar;
gen_fig = Z.gen_fig;
gen_tiff = Z.gen_tiff;
write_fig = 0;
brain = W.brain;
if combinedfig
    contrast_info = F.contrast_info_both;
    contrast_info_for_fig = F.contrast_info_both_for_fig;
else
    contrast_info = F.contrast_info;
    contrast_info_for_fig = F.contrast_info_for_fig;
end
pathn = F.pathn;
split = F.split;
T_map = F.T_map;
tstr = F.tstr;
%choose font for figure here
fontsize_choice = 16;
tick_number = 7;
%choose to overwrite max_T and min_T for map thresholds 
if cbar.colorbar_override == 1
    overwrite_map_thresholds = 1;
    max_T_choice = cbar.c_max;
    min_T_choice = cbar.c_min;
else
    overwrite_map_thresholds = 0;
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
    if ~combinedfig || tstr == 'F'
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
            sbar = linspace(min_T, max_T, cmap_res); %for colormap in 1-192 range
            %fraction of colorbar in each of cool, gray and hot areas:
            
            fcool = round((-th_z-min_T)/(max_T-min_T)*cmap_res);
            fgray = round(2*th_z/(max_T-min_T)*cmap_res);
            fhot = round((max_T-th_z)/(max_T-min_T)*cmap_res); 

            %Choose location of background brain as a function of th_z:
            %want it to be in between -th_z and +th_z
            %To ensure no overlap between background brain and
            %(de)activations, reduce th_z by a factor
            th_z_eff = th_z*0.9;
            maxb = max(brain(:)); %minb = min(brain(:)); %expect =0
            T_brain_over = (2*th_z_eff/maxb).*brain -th_z_eff;
            %add both negative and positive regions of interest
            T_brain_over(index_over) = T_map(index_over);
            if overwrite_map_thresholds
                %need to rescale T_brain_over to match the specified scale
                T_brain_over(T_brain_over>max_T) = max_T;
                T_brain_over(T_brain_over<min_T) = min_T;
                %required for unclear reason:
                T_brain_over(1,1) = max_T;
                T_brain_over(end,end) = min_T;
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
    if ~combinedfig  || tstr == 'F'
        split1 = split;
        colormap(split1); 
    else
        if tstr == 'T'
            %picking lower half of jet colormap; 
            split1 = [jet2(1:fcool,:);gray(fgray);hot(fhot)];
            colormap(split1);   
        end
    end
    axis(ax1, 'off')
    axis(ax1, 'image')
    hc1 = colorbar;
    if ~combinedfig || tstr == 'F'
        hc1_min = sbar(65);
        hc1_max = sbar(128);
    else
        if tstr == 'T'
            hc1_min = sbar(1);
            hc1_max = sbar(cmap_res);
        end
    end
    set(hc1, 'YLim', [hc1_min hc1_max]);    
    y_tick = linspace(hc1_min, hc1_max, tick_number)';
    set(hc1, 'YTick', y_tick);
    set(hc1, 'FontSize', fontsize_choice);
    %Customize here number of decimals
    set(hc1,'YTickLabel',sprintf('%.1f |',get(hc1,'YTick')'));  
    %cbfreeze(hc1);
    filen1 = fullfile(pathn,[tstr '_' str_cor '_' contrast_info '.fig']);
    if gen_fig
        saveas(fh1,filen1,'fig');
    end
    if gen_tiff
        filen2 = fullfile(pathn,[tstr '_' str_cor '_' contrast_info '.tiff']);
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
catch  exception
    disp(exception.identifier);
end
end