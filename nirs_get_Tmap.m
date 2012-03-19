function [min_T max_T min_T2 max_T2 T_brain_over fcool fgray fhot sbar] = ...
    nirs_get_Tmap(tstr,GroupColorbars,combinedfig,index_over,index_over2,...
    th_z,min_max_gap,s_map,overwrite_map_thresholds,brain,th_z_shrink,TC)
%Now we always use GroupColorbars = 0 and combinedfig = 1;
fcool = [];
fgray = [];
fhot = [];
min_T2 = 0;
max_T2 = 0;
if GroupColorbars
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
    max_T = TC.max_T_choice;
    min_T = TC.min_T_choice;
    if ~GroupColorbars
        max_T2 = TC.max_T_choice2;
        min_T2 = TC.min_T_choice2;
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
        if GroupColorbars
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