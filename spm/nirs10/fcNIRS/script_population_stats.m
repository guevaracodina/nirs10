%% script_population_stats
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
OP.path0 = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS';% data path
OP.path1 = 'F:\Edgar\Data\NIRS\epiNIRS_data\epiNIRS_Processed\epiNIRS';% analysis path
OP.pathbatch = 'F:\Edgar\Data\NIRS\epiNIRS_data\Batch_July2013'; %where to save the batches

%1) List of all subjects
subj{1} = 'epi101LH';
subj{2} = 'epi102FA';
subj{3} = 'epi103GJ';
% subj{4} = 'epi104MAL';
subj{5} = 'epi105MER';
subj{6} = 'epi106VLL';
subj{7} = 'epi107GC';
subj{8} = 'epi108DD';
subj{9} = 'epi109OC';
subj{10} = 'epiLFCDB'; %PP put this subject here for convenience
subj{11} = 'epi111ML';
subj{12} = 'epi112VV';
subj{13} = 'epi113AG';
subj{14} = 'epi114MCG';
subj{15} = 'epi115MD';
subj{16} = 'epi116GA';
subj{17} = 'epi117FDH';
subj{18} = 'epi118MC';
subj{19} = 'epi119MG';
subj{21} = 'epi121FB';
subj{22} = 'epi122PEV';
subj{23} = 'epi123JR';
subj{24} = 'epi124YP';
subj{25} = 'epi125YD';
subj{26} = 'epi126AB';
subj{27} = 'epi127SD';
subj{28} = 'epi128RW';
subj{29} = 'epi129DP';
subj{30} = 'epi130CO';
subj{31} = 'epi131AMA';
subj{32} = 'epi132FG';
subj{33} = 'epi133EA';
subj{34} = 'epi134ASC';
subj{35} = 'epi135MHT';
subj{36} = 'epi136JSL';
subj{37} = 'epi137SR';
subj{38} = 'epi138BG';
subj{39} = 'epi139JB';
subj{40} = 'epi140GG';
subj{41} = 'epi141JP';
subj{42} = 'epi142EG';
subj{43} = 'epi143MDM';
subj{44} = 'epi144TL';
subj{45} = 'epi145FS';
subj{46} = 'epi146MM'; %No T1 available, used epi127SD to run script, but only project on template
subj{47} = 'epi147LL';
subj{48} = 'epi148MHC';
subj{49} = 'epi149AHL';
n = length(subj);
OP.path_coreg = 'dataSPMa\coreg\';

%%
clc
epiAge      = nan([n 1]);
nSessions   = nan([n 1]);
nSources    = nan([n 1]);
nDetectors  = nan([n 1]);
for iSubjects = 1:n
    NIRSmat = fullfile(OP.path1, subj{iSubjects});
    NIRSmat = fullfile(NIRSmat, OP.path_coreg);
    NIRSmat = fullfile(NIRSmat,'NIRS.mat');
    % fprintf('%s\n',NIRSmat)
    try
        load(NIRSmat)
        epiAge(iSubjects)       = NIRS.Dt.s.age;
        nSessions(iSubjects)    = numel(NIRS.Dt.fir.pp(preProcStep).p);
        nSources(iSubjects)     =  NIRS.Cf.H.S.N;
        nDetectors(iSubjects)   =  NIRS.Cf.H.D.N;
    catch
        fprintf('%s not available\n', NIRSmat)
    end
end
fprintf('Mean age = %0.2fy, std. dev. = %0.2fy\n', nanmean(epiAge), nanstd(epiAge));
fprintf('Mean sessions = %0.2f, std. dev. = %0.2f\n', nanmean(nSessions), nanstd(nSessions));
fprintf('Mean Sources = %0.2f, std. dev. = %0.2f\n', nanmean(nSources), nanstd(nSources));
fprintf('Mean Detectors = %0.2f, std. dev. = %0.2f\n', nanmean(nDetectors), nanstd(nDetectors));
% EOF
