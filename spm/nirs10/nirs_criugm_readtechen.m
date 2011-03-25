function out = nirs_criugm_readtechen(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________


% This is the content of the nirs file that will be used for information
% on the experimental setup (measurement channels)
f = job.nirs_file;

switch int2str(job.system) % System used for the acquisition
    case '6'
        
        sDtp = job.sDtp;
        load(fullfile(sDtp,'NIRS.mat'));
        
        % Information about CW system
        NIRS.Cf.dev.n = 'CW6';
        NIRS.Cf.dev.wl = f.SD.Lambda;
        NIRS.Cf.dev.gn = f.systemInfo.gain;
        NIRS.Cf.dev.fs = 1/(f.t(2)-f.t(1));

        % Information about the Helmet
        NIRS.Cf.H.n = f.systemInfo.SDfilenm;
        NIRS.Cf.H.p = f.systemInfo.SDfilepath;
        
        switch job.coregType
            
            case 'Brainsight(c)'
        
                Cgen = zeros(0);    % [iC iS iD wl gp]

                Sn = NIRS.Cf.H.S.n;
                NS = NIRS.Cf.H.S.N;

                count = 1;      
                for iS =1:NS % for each source (in the BS file...)
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

                % Sort by wavelength, then by source number
                Cgen_s = sortrows(Cgen,[4 2]);

                % Number of channels (pairs)
                NIRS.Cf.H.C.N = size(Cgen_s,1); % length(Cid);
                % Measurement list: # of Pair / # of source / # of detector
                NIRS.Cf.H.C.id = Cgen_s(:,1:3)'; % Cid;
                % Wavelength
                NIRS.Cf.H.C.wl = Cgen_s(:,4)';   % Cwl;
                % Source-detector distance (usually 2D)
                NIRS.Cf.H.C.gp = Cgen_s(:,5)';   % Cgp;


                
            %case 'T1_vitamins'
                
            otherwise % In any other case, when no BrainSight file is provided,
                % all the setup info is taken from the nirs file (until
                % enventually some coregistration (e.g., from vitamin
                % markers in the T1) is done).
                
                % Sources
                %NIRS.Cf.H.S.n
                NIRS.Cf.H.S.N = f.SD.nSrcs;
                %NIRS.Cf.H.S.void
                NIRS.Cf.H.S.r.o.mm.p = f.SD.SrcPos*10; % nSrcs x 3, cm->mm
                
                % Detectors
                %NIRS.Cf.H.D.n
                NIRS.Cf.H.D.N = f.SD.nDets;
                %NIRS.Cf.H.D.void
                NIRS.Cf.H.D.r.o.mm.p = f.SD.DetPos*10; % nSrcs x 3, cm->mm
                
                % All points
                NIRS.Cf.H.P.N = f.SD.nSrcs + f.SD.nDets;
                NIRS.Cf.H.P.r.o.mm.p = [f.SD.SrcPos; f.SD.DetPos]*10; % nPts x 3, cm->mm
                
                % Channels (pairs)
                NIRS.Cf.H.C.N = size(f.SD.MeasList,1);
                
                % Measurement list: # of Pair / # of source / # of detector
                ml = f.SD.MeasList;       
                % Add a column for pair ID
                ml = [1:size(ml,1); ml']';
                % Sort by wavelength, then source index
                ml = sortrows(ml,[5 2]);
                NIRS.Cf.H.C.id = ml(:,[1 2 3])'; % 3 x nChannels
                
                % Wavelength
                NIRS.Cf.H.C.wl = ml(:,5)';
                
                % Source-detector distance (usually 2D)
                distSD = zeros(NIRS.Cf.H.P.N,1);
                for iC = 1:NIRS.Cf.H.P.N
                    iS = ml(iC,2);   % # source
                    iD = ml(iC,3);   % # détecteur
                    distSD(iC,1) = sqrt(sum((f.SD.SrcPos(iS,:) - f.SD.DetPos(iD,:)).^2)) ;
                end
                NIRS.Cf.H.C.gp = distSD'; % 1 x nChannels, cm
                
        end
        
    otherwise
        disp('Reading of Techen file (.nirs) failed. For now only CW6 output files are supported.');
        
end

save(fullfile(NIRS.Dt.s.p,'NIRS.mat'),'NIRS');
out.NIRSmat{1} = fullfile(NIRS.Dt.s.p,'NIRS.mat');

return
