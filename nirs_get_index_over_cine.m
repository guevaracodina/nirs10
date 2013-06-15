function [th_z str_cor index_over index_over2] = nirs_get_index_over_cine(erdf,s_map,p_value)
%uncorrected
str_cor = 'unc';
%threshold
th_z = spm_invTcdf(1-p_value, erdf);
index_over = find(s_map > th_z);
index_over2 = [find(s_map > th_z); find(-s_map > th_z)];