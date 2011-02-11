function out = nirs_criugm_readtechen(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________


% Channels
f = job.nirs_file;

switch int2str(job.system)
    case '6'
        
        sDtp = job.sDtp;
        load(fullfile(sDtp,'NIRS.mat'));
        
        Cgen = zeros(0);    % [iC iS iD wl gp]
        
        Sn = NIRS.Cf.H.S.n;
        NS = NIRS.Cf.H.S.N;
        
        count = 1;
        for iS =1:NS
            Snum = Sn{1,iS};
            Snum = Snum(2:end);
            C_Snum(1,:) = (f.SD.MeasList(:,1)==str2num(Snum))';
            for i=1:size(C_Snum,2)
                if C_Snum(1,i)==1
                    Cgen(count,1) = i;
                    Cgen(count,2) = iS;
                    Cgen(count,3) = f.SD.MeasList(i,2);
                    Cgen(count,4) = f.SD.MeasList(i,4);
                    % Ne devrait-on pas corriger pour la courbure ?????
                    Cgen(count,5) = sqrt(sum((f.SD.SrcPos(iS,:) - f.SD.DetPos(Cgen(count,3),:)).^2));
                    count = count+1;
                end
            end
        end
        
        % Sort by wavelength
        Cgen_s = sortrows(Cgen,[4 2]);
        
        NIRS.Cf.H.C.N = size(Cgen_s,1); % length(Cid);
        NIRS.Cf.H.C.id = Cgen_s(:,1:3)'; % Cid;
        NIRS.Cf.H.C.wl = Cgen_s(:,4)';   % Cwl;
        NIRS.Cf.H.C.gp = Cgen_s(:,5)';   % Cgp;
        
        % Information about CW system
        NIRS.Cf.dev.n = 'CW6';
        NIRS.Cf.dev.wl = f.SD.Lambda;
        NIRS.Cf.dev.gn = f.systemInfo.gain;
        NIRS.Cf.dev.fs = 1/0.04;
        
        % Information about the Helmet
        NIRS.Cf.H.n = f.systemInfo.SDfilenm;
        NIRS.Cf.H.p = f.systemInfo.SDfilepath;
        
    otherwise
        disp('Read Techen file failed');
end

save(fullfile(NIRS.Dt.s.p,'NIRS.mat'),'NIRS');
out.NIRSmat{1} = fullfile(NIRS.Dt.s.p,'NIRS.mat');
return