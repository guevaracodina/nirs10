function save_sens(save_path, sens_data)
%**************************************************************************
%Function: this function saves sensitivity matrix in binary flow
%Editor: Ke Peng, Laboratoire d'Imagerie Optique et Moleculaire
%Data: March 05, 2012
%**************************************************************************

try
    fid = fopen(fullfile(save_path,'sens.bin'),'wb');
    fwrite(fid,sens_data,'double');
    fclose(fid);
catch
    disp('Can not save sensitivity matrix in binary flow');
end

