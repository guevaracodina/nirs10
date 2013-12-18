%% Clustering/linkage script with NIRS data
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
%% Read a random subject
clear; clc
% Le cas 104MAL est extrêmement inhabituel. Il ne pourra pas être utilisé
% pour les analyses.
% epi127SD Good optodes coverage
NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi127SD\dataSPMa\coreg\NIRS.mat';
% epi143MDM  Good optodes coverage
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS\epi143MDM\dataSPMa\coreg\NIRS.mat';
% Control subject 1
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet001\dataSPMa\coreg\NIRS.mat';
% Control subject 2
% NIRSmat = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS\Sujet002\dataSPMa\coreg\NIRS.mat';
load(NIRSmat)

%% Read all files from the 4th processing level ODtoHbOHbR
% 1) readBOXY
% 2) remove_chn_stdev
% 3) normalize_baseline
% 4) ODtoHbOHbR
preProcStep = 4;
% switch chromophore
%     case 1
%         hb = 'HbO';
%     case 2
%         hb = 'HbR';
%     case 3
%         hb = 'HbT';
% end
hb = 2;
% nSessions = numel(NIRS.Dt.fir.pp(preProcStep).p);
nSessions = numel(NIRS.Dt.fir.Sess);
dataNIRS = [];
for iFiles = 1:nSessions,
    % currentData [nChannels x nTimePoints]
    currentData = fopen_NIR(NIRS.Dt.fir.pp(preProcStep).p{iFiles}, NIRS.Cf.H.C.N);
    % dataNIRS [nChannels x (nTimePoints x nSessions)]
    dataNIRS = [dataNIRS, currentData];
end
dataNIRS = dataNIRS(NIRS.Cf.H.C.wl == hb,:);
% Time vector
t = linspace(0, (size(dataNIRS,2)-1)/NIRS.Cf.dev.fs , size(dataNIRS,2));

%% Display data
h = figure; set(h,'color','w')
set(h,'Name',sprintf('Pre-processing step: %s', NIRS.Dt.fir.pp(preProcStep).pre))
plot(t, dataNIRS');
ylabel('[\muM]')
xlabel('t [s]')
switch hb
    case 1
        titleString = 'HbO';
    case 2
        titleString = 'HbR';
    case 3
        titleString = 'HbT';
end
title(titleString)

%% Using clusterdata (hierarchical clustering)
[PERM1 T2 corrMat] = nirs_clustergram(dataNIRS);

%% Print clustergram
% set(gcf,'color',[1 1 1],'paperpositionmode','auto');
% [pathName,~,~] = fileparts(NIRSmat);
% [~, subjectName, ~] = fileparts(NIRS.Dt.s.p);
% print(gcf,'-dpng',fullfile(pathName,[subjectName '_corrMat_RAW']));
% EOF
