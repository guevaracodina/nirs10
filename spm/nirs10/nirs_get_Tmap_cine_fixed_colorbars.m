function [min_T max_T min_T2 max_T2 T_brain_over fcool fgray fhot sbar] = ...
    nirs_get_Tmap_cine_fixed_colorbars(index_over,index_over2,...
    th_z,min_max_gap,s_map,overwrite_map_thresholds,brain,th_z_shrink,TC)
%Positive responses
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
%Negative responses
if ~isempty(index_over2)
    min_T2 = min(s_map(index_over2));
    max_T2 = max(s_map(index_over2));
    if min_T2 == max_T2 %for example if only one data point in index_over2
        min_T2 = min_T2 - min_max_gap;
    end
else
    max_T2 = th_z+min_max_gap;
    min_T2 = th_z;
end

if overwrite_map_thresholds
    max_T = TC.max_T_choice;
    min_T = TC.min_T_choice;
    max_T2 = TC.max_T_choice2;
    min_T2 = TC.min_T_choice2;
end

%colormap resolution:
cmap_res = 3*64;
sbar = linspace(min_T2, max_T, cmap_res);
fcool = round((-th_z-min_T2)/(max_T-min_T2)*cmap_res);
fgray = round(2*th_z/(max_T-min_T2)*cmap_res);
fhot = round((max_T-th_z)/(max_T-min_T2)*cmap_res);
th_z_eff = th_z*th_z_shrink;
maxb = max(brain(:)); %minb = min(brain(:)); %expect =0
T_brain_over = (2*th_z_eff/maxb).*brain -th_z_eff;
%add both negative and positive regions of interest
%if ~isempty(index_over)
    T_brain_over(index_over) = s_map(index_over);
%end
%if ~isempty(index_over2)
    T_brain_over(index_over2) = s_map(index_over2);
%end
if overwrite_map_thresholds
    T_brain_over(T_brain_over>max_T) = max_T;
    T_brain_over(T_brain_over<min_T2) = min_T2;
    T_brain_over(1,1) = max_T;
    T_brain_over(end,end) = min_T2;
end
