function [sens_load, number] = load_sens(load_path,sens_row_number,sens_column_number)
%**************************************************************************
%funciton: This function is to read the sensitivity matrix in binary mode
%Editor: Ke Peng, Loboratoire d'Imagerie Optique et Moleculaire
%Date: March 5th, 2012
%**************************************************************************

try
    fid = fopen(fullfile(load_path,'sens.bin'),'rb');
    [sens_load, number] = fread(fid,[sens_row_number,sens_column_number],'double');
    fclose(fid);
catch
    disp('Cannot load sensitivity matrix in binary mode');
end




