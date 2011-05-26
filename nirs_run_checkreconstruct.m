function out = nirs_run_checkreconstruct(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clement Bonnery
Kas =1;
switch Kas
    case 1
        out_Hbs = job.outreconstruct_Hb;
        
        load(job.NIRSmat{:});
        fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
        
        dec_DHbO = zeros(1,length(fnirs.d));
        dec_DHbR = zeros(1,length(fnirs.d));
        
        for iHbs =1:size(out_Hbs,1)
            f = out_Hbs{iHbs,1};
            V = spm_vol(f);
            Y = spm_read_vols(V);
            
            [~,name,~] = fileparts(f);
            sep =strfind(name,'_');
%             timee =str2num(name(sep(2)+2:sep(3)-1));
            timee =str2num(name(sep(2)+2:end));
            if ~isempty(strfind(name,'HbO'))
                
                dec_DHbO(1,timee:timee+25) = Y(11,8,8);
            else
                dec_DHbR(1,timee:timee+25) = Y(11,8,8);
            end
        end
        figure;
        subplot(2,1,1)
        plot(dec_DHbO,'r');
        subplot(2,1,2)
        plot(dec_DHbR);
        
    case 2
        out_muas = job.outreconstruct_Hb;
        
        load(job.NIRSmat{:});
        fnirs = load(NIRS.Dt.fir.pp.p{:},'-mat');
        
        dec_Dmua = zeros(1,length(fnirs.d));
        
        for imuas =1:size(out_muas,1)
            f = out_muas{imuas,1};
            V = spm_vol(f);
            Y = spm_read_vols(V);
            
            [~,name,~] = fileparts(f);
            sep =strfind(name,'_');
            timee =str2num(name(sep(2)+2:sep(3)-1));
                
            dec_Dmua(1,timee:timee+25) = Y(11,8,7);
        end
        figure;
        plot(dec_Dmua(1,5000:7000));
end
out =1;
end