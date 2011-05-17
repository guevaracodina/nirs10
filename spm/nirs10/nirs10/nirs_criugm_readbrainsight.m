function out = nirs_criugm_readbrainsight(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

staxp = job.staxp;
sDtp = job.sDtp;

text_brainsight = fopen(staxp,'r'); %open the file for reading
% load(fullfile(sDtp,'NIRS.mat')); % CB: absurde ! la matrice n'a pas
% encore ete creee !!! (si le code brise a cet endroit : en parler avec moi)

On = {};
Sn = {};
Dn = {};
Op_rom = [];
Sp_rom = [];
Dp_rom = [];

line1 = fgetl(text_brainsight);

switch line1
    case '# Coordinate system: Brainsight' % Brainsight file from CW5 room
        %         headerlines = 3; % First lines are dedicated to format and orientation conventions
        %         linesCell = textscan(text_brainsight, '%s', 'delimiter', '\n', 'whitespace', '', 'headerLines', headerlines);
        %         linesCell = linesCell{:};
        %
        %         for i = 1:length(linesCell)
        %             % lecture de la ligne
        %             line = textscan(linesCell{i},'%s %f %f %f');
        %
        %             % Annotation venant de Brainsight de Rogue-Research Inc.
        %             Stringletter = char(line{1});
        %
        %             if Stringletter(1) == 'S'
        %                 % les coordonnees correspondent a une source
        %                 srcNameArray = [srcNameArray ; sprintf('%s', Stringletter(1:end-1))];
        %                 src = [src ; line{2} line{3} line{4}];
        %             elseif Stringletter(1) == 'D'
        %                 % les coordonnees correspondent a un detecteur
        %                 detNameArray = [detNameArray ; sprintf('%s', Stringletter(1:end-1))];
        %                 det = [det ; line{2} line{3} line{4}];
        %             else
        %                 markersNameArray = [markersNameArray ; sprintf('%s', Stringletter(1:end-1))];
        %                 markers = [markers ; line{2} line{3} line{4}];
        %             end
        %         end
        
    case '# Electrodes' % Brainsight file from CW6 room
        headerlines = 3;
        lines_cell = textscan(text_brainsight, '%s', 'delimiter', '\n', 'whitespace', '', 'headerLines', headerlines);
        lines_cell = lines_cell{:};
        
        for i = 1:length(lines_cell)/2
            line = textscan(lines_cell{2*i-1},'%s');
            line = char(line{1});
            
            if line(1) == 'S' % sources
                Sn = [Sn ; sprintf('%s', line)];
                Sp_rom  = [Sp_rom ; textscan(lines_cell{2*i},'%f %f %f')];
            elseif line(1) == 'D' % detectors
                Dn = [Dn ; sprintf('%s', line)];
                Dp_rom  = [Dp_rom ; textscan(lines_cell{2*i},'%f %f %f')];
            else % anatomical markers
                %O
                On = [On ; sprintf('%s', line)];
                Op_rom  = [Op_rom ; textscan(lines_cell{2*i},'%f %f %f')];
            end
        end
        
    otherwise % error('Bad file format: this is not a Brainsight markers file');
        %headerlines = -1;
        msgbox('Format error : the specified file doesn''t fit with markers file from Brainsight.');
        % sortir proprement... voir les erreurs qu'on peut  lancer avec
        % spm%%%%%
end

fclose(text_brainsight);

%fill in NIRS
NIRS=[];
NIRS.Cf.H.O.n = On';
NIRS.Cf.H.S.n = Sn';
NIRS.Cf.H.D.n = Dn';
NIRS.Cf.H.O.r.o.mm.p = Op_rom';
NIRS.Cf.H.S.r.o.mm.p = Sp_rom';
NIRS.Cf.H.D.r.o.mm.p = Dp_rom';
NIRS.Cf.H.O.N = size(NIRS.Cf.H.O.n,2);
NIRS.Cf.H.S.N = size(NIRS.Cf.H.S.n,2);
NIRS.Cf.H.D.N = size(NIRS.Cf.H.D.n,2);
NIRS.Cf.H.P.N = NIRS.Cf.H.S.N + NIRS.Cf.H.D.N + NIRS.Cf.H.O.N;

% save([sDtp 'NIRS'],'NIRS');
out = NIRS;
% out.NIRSpath{1} = fullfile(sDtp,'NIRS.mat'); %why not fullfile(NIRS.subj_path,'NIRS.mat');?
return