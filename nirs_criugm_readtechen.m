function out = nirs_criugm_readtechen(job)
% This is the content of the nirs file that will be used for information
% on the experimental setup (measurement channels)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

f = job.nirs_file;
out = 0;
if ~isfield(f,'systemInfo')
    [dir_spm, dummy1,dummy2] = fileparts(which('spm'));
    load(fullfile(dir_spm, 'toolbox','nirs10','nirs10_templates','systemInfo.mat'));
    f.systemInfo = systemInfo;
    load(fullfile(dir_spm, 'toolbox','nirs10','nirs10_templates','SD.mat'));
    f.SD = SD;
end
switch int2str(job.system) % System used for the acquisition
    case '6'
        
        sDtp = job.sDtp;
        NIRS = job.NIRS;
        
        
        NIRS.Cf.dev.n = 'CW6';
        NIRS.Cf.dev.wl = f.SD.Lambda;
        % Information about auxilaries
        NIRS.Dt.aux = f.aux;
        if isfield(job.biopac, 'choice_biopac')
            biopac.n = job.biopac.choice_biopac;
            max.aux = 4;
            NIRS.Dt.aux.eprime = f.aux(:,1:max.aux-biopac.n);
            NIRS.Dt.aux.biopac = f.aux(:,biopac.n+1:max.aux);
        else
            NIRS.Dt.aux.eprime = f.aux;
        end
        
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
                
                % % Sort by channel number
                % Cgen_s = sortrows(Cgen,[4 2]);
                Cgen_s = sortrows(Cgen,1);
                
                %%%%%%%%%%%%%%%%%%%% ATTENTION, IL FAUT QUE d SOIT AUSSI
                %%%%%%%%%%%%%%%%%%%% CHANGE SINON LES PAIRES NE
                %%%%%%%%%%%%%%%%%%%% CORRESPONDENT PLUS.......
                
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
                
                % % !! DO NOT DO THIS !! KEEP ORDER OF PAIRS !!
                % % Sort by wavelength, then source index
                % ml = sortrows(ml,[5 2]);
                
                NIRS.Cf.H.C.id = ml(:,[1 2 3])'; % 3 x nChannels
                
                % Wavelength
                NIRS.Cf.H.C.wl = ml(:,5)';
                
                % Source-detector distance (usually 2D)
                distSD = zeros(NIRS.Cf.H.C.N,1);
                for iC = 1:NIRS.Cf.H.C.N
                    iS = ml(iC,2);   % # source
                    iD = ml(iC,3);   % # détecteur
                    distSD(iC,1) = sqrt(sum((f.SD.SrcPos(iS,:) - f.SD.DetPos(iD,:)).^2)) ;
                end
                NIRS.Cf.H.C.gp = distSD'; % 1 x nChannels, cm
                
        end
        
    case {'5','1'}  %System 1: Imaginc
        sDtp = job.sDtp;
        NIRS = job.NIRS;
        
        % Information about CW system
        if int2str(job.system) == '1'
            NIRS.Cf.dev.n = 'Imaginc';
            NIRS.Cf.dev.wl = [735 850];
        else
            % Information about CW system
            NIRS.Cf.dev.n = 'CW5';
            NIRS.Cf.dev.wl = f.SD.Lambda;
            %         NIRS.Cf.dev.gn = f.systemInfo.gain;
        end
         NIRS.Cf.dev.fs = 1/(f.t(2)-f.t(1));
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
                
                % % Sort by channel number
                % Cgen_s = sortrows(Cgen,[4 2]);
                Cgen_s = sortrows(Cgen,1);
                
                % Number of channels (pairs)
                NIRS.Cf.H.C.N = size(Cgen_s,1); % length(Cid);
                % Measurement list: # of Pair / # of source / # of detector
                NIRS.Cf.H.C.id = Cgen_s(:,1:3)'; % Cid;
                % Wavelength
                NIRS.Cf.H.C.wl = Cgen_s(:,4)';   % Cwl;
                % Source-detector distance (usually 2D)
                NIRS.Cf.H.C.gp = Cgen_s(:,5)';   % Cgp;
                
            otherwise % In any other case, when no BrainSight file is provided,
                % all the setup info is taken from the nirs file (until
                % enventually some coregistration (e.g., from vitamin
                % markers in the T1) is done).
                
                % Sources
                NIRS.Cf.H.S.N = size(f.SD.SrcPos,1);
                NIRS.Cf.H.S.r.o.mm.p = (f.SD.SrcPos*10)'; % nSrcs x 3, cm->mm
                
                f.SD.DetPos(3,2) =9;
                % Detectors
                NIRS.Cf.H.D.N = size(f.SD.DetPos,1);
                NIRS.Cf.H.D.r.o.mm.p = (f.SD.DetPos*10)'; % nSrcs x 3, cm->mm
                
                % All points
                NIRS.Cf.H.P.N = NIRS.Cf.H.S.N + NIRS.Cf.H.D.N;
                %                 NIRS.Cf.H.P.r.o.mm.p = [NIRS.Cf.H.S.r.o.mm.p NIRS.Cf.H.D.r.o.mm.p]; % nPts x 3, cm->mm
                
                % Channels (pairs)
                Cgen = zeros(0);    % [iC iS iD wl gp]
                NS = NIRS.Cf.H.S.N;
                
                count = 1;
                for iS =1:NS % for each source (in the BS file...)
                    C_Snum(1,:) = (f.ml(:,1)==iS)';
                    for i=1:size(C_Snum,2)
                        if C_Snum(1,i)==1
                            Cgen(count,1) = i;
                            Cgen(count,2) = iS;
                            Cgen(count,3) = f.ml(i,2);
                            Cgen(count,4) = f.ml(i,4);
                            % Ne devrait-on pas corriger pour la courbure ?????
                            indD = mod(i,NIRS.Cf.H.D.N);
                            if indD==0, indD =NIRS.Cf.H.D.N; end
                            Cgen(count,5) = sqrt(sum((f.SD.SrcPos(iS,:) - f.SD.DetPos(indD,:)).^2));
                            count = count+1;
                        end
                    end
                end
                
                % % Sort by channel number
                % Cgen_s = sortrows(Cgen,[4 2]);
                NIRS.Cf.H.C.order = Cgen(:,1)';
                Cgen_s = sortrows(Cgen,1);
                
                % Number of channels (pairs)
                NIRS.Cf.H.C.N = size(Cgen_s,1); % length(Cid);
                % Measurement list: # of Pair / # of source / # of detector
                NIRS.Cf.H.C.id = Cgen_s(:,1:3)'; % Cid;
                % Wavelength
                NIRS.Cf.H.C.wl = Cgen_s(:,4)';   % Cwl;
                % Source-detector distance (usually 2D)
                NIRS.Cf.H.C.gp = Cgen_s(:,5)';   % Cgp;
        end
        
    otherwise
        disp('Reading of Techen file (.nirs) failed. For now only CW6 output files are supported.');
end

% Return updated NIRS matrix (or 0 if fail)
out = NIRS;

return
