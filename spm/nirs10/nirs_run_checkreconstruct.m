function out = nirs_run_checkreconstruct(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clement Bonnery

load(job.NIRSmat);
out_Hbs = job.outreconstruct_Hb;

dec_DHbO = zeros(1,size(out_Hbs,1));
dec_DHbR = zeros(1,size(out_Hbs,1));

for iHbs =1:size(out_Hbs,1)
    f = out_Hbs{iHbs,1};
    V = spm_vol(f);
    Y = spm_read_vols(V);
    
    [~,name,~] = fileparts(f);
    sep =strfind(name,'_');
    timee =name(sep(2)+2:sep(3)-1);
    if ~isempty(strfind(name,'HbO'))
        dec_DHbO(1,timee) = Y(5,4,7);
    else
        dec_DHbR(1,timee) = Y(5,4,7);
    end
end
figure;
subplot(2,1,1)
plot(dec_DHbO,'r');
subplot(2,1,2)
plot(dec_DHbR);

out =1;
end