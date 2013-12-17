%% Main script for controls in epilepsy fc study
clc
OP.path0 = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS';% data path
OP.path1 = 'F:\Edgar\Data\NIRS\epiNIRS_data\ControlNIRS_Processed\ControlNIRS';% analysis path
OP.pathbatch = 'F:\Edgar\Data\NIRS\epiNIRS_data\Batch_December2013'; %where to save the batches
%1) List of all subjects

subj{1} = 'Sujet001';
subj{2} = 'Sujet002';
% subj{3} = 'epi103GJ';
% subj{4} = 'epi104MAL';
% subj{5} = 'epi105MER';
% subj{6} = 'epi106VLL';
% subj{7} = 'epi107GC';
% subj{8} = 'epi108DD';
% subj{9} = 'epi109OC';
% subj{10} = 'epiLFCDB'; %PP put this subject here for convenience
% subj{11} = 'epi111ML';
% subj{12} = 'epi112VV';
% subj{13} = 'epi113AG';
% subj{14} = 'epi114MCG';
% subj{15} = 'epi115MD';
% subj{16} = 'epi116GA';
% subj{17} = 'epi117FDH';
% subj{18} = 'epi118MC';
% subj{19} = 'epi119MG';
% subj{21} = 'epi121FB';
% subj{22} = 'epi122PEV';
% subj{23} = 'epi123JR';
% subj{24} = 'epi124YP';
% subj{25} = 'epi125YD';
% subj{26} = 'epi126AB';
% subj{27} = 'epi127SD';
% subj{28} = 'epi128RW';
% subj{29} = 'epi129DP';
% subj{30} = 'epi130CO';
% subj{31} = 'epi131AMA';
% subj{32} = 'epi132FG';
% subj{33} = 'epi133EA';
% subj{34} = 'epi134ASC';
% subj{35} = 'epi135MHT';
% subj{36} = 'epi136JSL';
% subj{37} = 'epi137SR';
% subj{38} = 'epi138BG';
% subj{39} = 'epi139JB';
% subj{40} = 'epi140GG';
% subj{41} = 'epi141JP';
% subj{42} = 'epi142EG';
% subj{43} = 'epi143MDM';
% subj{44} = 'epi144TL';
% subj{45} = 'epi145FS';
% subj{46} = 'epi146MM'; %No T1 available, used epi127SD to run script, but only project on template
% subj{47} = 'epi147LL';
% subj{48} = 'epi148MHC';
% subj{49} = 'epi149AHL';
n = length(subj);
OP.run_label = 'a';

%35 subjects with identified focus
% focus{1} = [20.5000000000000,55.8000000000000,26.5000000000000]; %PP: changed 1st coordinate to Right side 
% focus{3} = [33.5000000000000,-10.5000000000000,72.5000000000000]; %PP: changed 1st coordinate to Right side
% %focus{4} = occipital
% focus{5} = [-8,7.50000000000000,114]; %PP: changed 1st coordinate to Left side 
% focus{6} = [-40.9000000000000,41.5000000000000,77.1000000000000]; %OK
% focus{7} = [-51.5000000000000,-42.8000000000000,28.5000000000000]; %OK
% focus{8} = [16,6.50000000000000,52]; %multifoci
% focus{9} = [45.5000000000000,17.9000000000000,35.5000000000000]; %PP: changed 1st coordinate to Right side 
% focus{10} = [-35.9 -0.7 25.6]; %PP -- could check with Dang
% focus{11} = [45.5000000000000,30.9000000000000,56.5000000000000]; %OK
% focus{12} = [34.5000000000000,10,4.50000000000000]; %PP: changed 1st coordinate to Right side 
% focus{13} = [38.8000000000000,47.1000000000000,7]; %PP: changed 1st coordinate to Right side 
% focus{14} = [-47.5000000000000,-18.9000000000000,25.5000000000000]; %OK
% focus{15} = [-36.8 -10.0 7.5]; %PP -- could check with Dang
% focus{17} = [-16.9 5.0 23.6] ; % [83,6.50000000000000,14];  %this looks wrong -- needs to be redone
% focus{18} = [47.8000000000000,2.80000000000000,16.4000000000000]; %PP: changed 1st coordinate to Right side 
% focus{19} = [29.5000000000000,-14.6000000000000,54.4000000000000]; %PP: changed 1st coordinate to Right side 
% focus{21} = [31,-66.2000000000000,84];
% focus{22} = [-48.5000000000000,-24.3000000000000,-5.20000000000000]; %PP: changed 1st coordinate to Left side 
% focus{24} = [-42.7000000000000,-11.2000000000000,-6.10000000000000]; %PP: changed 1st coordinate to Left side 
% focus{25} = [41.3 14.2 50.2]; %PP -- could check with Dang
% focus{27} = [-33.7 -0.1 54.2]; %PP modified, was: [-31.4000000000000,-21.7000000000000,23.9000000000000];
% focus{28} = [54.6000000000000,-12.1000000000000,18.9000000000000]; %PP: changed 1st coordinate to Right side 
% focus{29} = [-35.1000000000000,15.7000000000000,75.4000000000000];
% focus{30} = [-51.5000000000000,-28,20.2000000000000]; %PP: changed 1st coordinate to Left side 
% focus{31} = [47.8000000000000,-27.4000000000000,-9.30000000000000]; %PP: changed 1st coordinate to Right side 
% focus{32} = [27.9000000000000,34,-4.40000000000000]; %PP: changed 1st coordinate to Right side 
% focus{33} = [38.8000000000000,48.6000000000000,23.5000000000000];
% focus{34} = [-35.5000000000000,-52.9000000000000,26.2000000000000];
% focus{35} = [39.5000000000000,25.9000000000000,2.50000000000000]; %PP: changed 1st coordinate to Right side 
% focus{36} = [34.6000000000000,5.20000000000000,34.3000000000000]; %PP: changed 1st coordinate to Right side 
% focus{37} = [40.3000000000000,34.2000000000000,9.80000000000000]; %PP: changed 1st coordinate to Right side 
% focus{39} = [51.5000000000000,-9.50000000000000,33.5000000000000]; %PP: changed 1st coordinate to Right side 
% focus{40} = [-31.8000000000000,26.4000000000000,10.3000000000000];
% focus{42} = [38.5000000000000,-38.8000000000000,16.5000000000000]; %PP: changed 1st coordinate to Right side 
% focus{43} = [34.7 12.8 43.4]; %PP -- could check with Dang
%focus{44} = occipital
%focus{45} = occipital
%focus{46} = occipital
%focus{47} = occipital
%focus{48} = occipital
%focus{49} = occipital
% OP.focus = focus;
OP.coreg_size = 20; % Change to 30 or something bigger to have a better coregistration
OP.test_hrf = 0;
% Loop over subjects
for i=1:n, 
    mux{i} = 32; 
    views{i} = [4 3 2 5];
    %Boolean for specifying fiducials of subject (unnormalized) -- otherwise use fiducials on normalized template 
    subj_fiducials{i} = 1; 
%     focus_size{i} = 18;
    channel_pca{i} = 1;
    nComponents{i} = 5;
    volt{i} = 1; % 2;
    project_subject{i} = 1;
end
%for occipitals
% occipital_only = 1;
% if occipital_only
%     for i = [34 42 44:49]
%         views{i} = 6;
%     end
%     views{4} = 2;
% else
%     for i = [42 44:49]
%         views{i} = [6 4 3 2 5];
%     end
%     views{4} = [2 4 3 5];
% end

% mux{1} = 16;
% focus_size{27} = 36; %big focus
% 
% subj_fiducials{5} = 0; project_subject{5} = 0;
% subj_fiducials{8} = 0; project_subject{8} = 0;
% subj_fiducials{17} = 0; project_subject{17} = 0;
% subj_fiducials{21} = 0; project_subject{21} = 0; %For those subjects whose fiducial coordinates are not normal, we project them onto the template
% 
% subj_fiducials{46} = 0; project_subject{46} = 0;
% subj_fiducials{48} = 0;

% Do not use custom fiducials if 0
subj_fiducials{1} = 0;
subj_fiducials{2} = 0;

% Fcoord{1} = [-2.7 82.1 -27.7; -67.9 -13.4 -28.5; 64.4 -9.8 -28.5];
% Fcoord{2} = [-4.9 102.7 -22.3; -84.9 -6.7 -22.3; 77.8 -16.4 -36.6];
% Fcoord{1} = [-2.7 82.1 -27.7; -67.9 -13.4 -28.5; 64.4 -9.8 -28.5];
% Fcoord{2} = [-4.9 102.7 -22.3; -84.9 -6.7 -22.3; 77.8 -16.4 -36.6];
% Fcoord{3} = [0.9 95.4 -41.0; -61.7 -5.4 -27.7; 58.1 -8.0 -25.0];
% Fcoord{4} = [-0.7 137.9 -10.2; -79.1 30.7 -18.2; 77.7 24.1 -11.1];
% Fcoord{5} = [46.3 83.6 80.6; -26.6 -3.1 33.4; 109.4 -12.4 32.5]; %changed slightly
% Fcoord{6} = [-2.1 114.4 -9.4; -69.8 2.3 -16.7; 72.9 11.4 -15.8]; %[-2.6 119.6 -4.2; -76.2 3.4 9.4; 81.2 4.7 8.2];
% Fcoord{7} = [-1.3 87.4 -24.1; -66.4 -16.9 -24.1; 60.2 -16.9 -30.3]; %[-2.4 87.4 -17.0; -69.2 -19.4 -9.7; 66.8 -26.7 -4.9];
% Fcoord{8} = [44.8 87.5 35.0; -13.8 -10.5 32.6; 108.3 -8.1 24.0]; %[44.3 89.1 41.2; -17.6 -16.5 32.7; 118.2 -16.5 31.5]; %OK, but could still move down slightly
% Fcoord{9} = [0.0 94.6 -1.8; -69.8 -0.0 -66.9; 68.0 2.7 -66.9]; %Large change %[-0.6 93.4 -1.2; -72.6 -17.0 -21.8; 75.1 -31.5 -21.8];
% Fcoord{11} = [2.2 104.4 -0.9; -76.9 -3.6 -13.4; 79.6 -1.8 -18.7]; %small change %[2.4 105.6 -6.1; -79.5 -16.6 -12.1; 81.9 -12.8 -12.1];
% Fcoord{12} = [0.4 92.8 -32.1; -73.3 -18.2 -52.6; 76.0 -15.6 -41.0]; %large change %[0.0 93.4 -32.8; -76.6 -24.3 -26.7; 76.5 -21.8 -18.2];
% Fcoord{13} = [0.0 95.3 -25.9; -68.7 -15.1 -51.0; 69.6 -8.4 -53.5]; %large change %[2.3 96.7 -26.2; -73.0 -15.9 -26.2; 77.6 -14.8 -27.3];
% Fcoord{14} = [-1.8 105.3 -25.9; -65.3 -5.4 -43.7; 68.8 -5.4 -43.7]; %large change %[0.0 108.0 -21.8; -68.6 -9.7 -24.3; 67.3 -12.1 -21.8];
% Fcoord{15} = [0.6 87.4 11.1; -71.3 -18.8 -43.1; 69.1 -12.2 -40.9];
% Fcoord{16} = [-5.6 93.0 -11.3; -63.0 -6.1 -46.1; 53.4 -2.6 -53.9]; %large change %[-3.5 91.0 -9.5; -67.3 -14.2 -22.5; 64.9 -8.3 -27.2];
% Fcoord{17} = [35.5 112.1 16.4; -41.8 -1.3 -38.0; 119.9 8.1 -24.6]; %very large change -- might be better to display on template since segmentation is very bad %[35.9 114.5 19.9; -45.5 2.9 5.3; 123.3 2.9 1.7];
% Fcoord{18} = [0.0 89.5 -50.2; -68.1 -6.7 -47.7; 66.4 -5.0 -45.2]; %large change %[-1.7 89.9 -43.2; -71.3 -2.9 -25.0; 70.2 -5.1 -21.6];
% Fcoord{19} = [4.4 85.7 7.0; -69.2 -11.4 -26.2; 70.1 -11.4 -32.4]; %small change %[5.4 89.0 8.3; -70.9 -9.0 -23.8; 76.8 -16.1 -36.9];
% Fcoord{21} = [42.6 104.9 20.9; -34.7 -4.4 -17.5; 118.2 -2.6 -8.6]; %large change %[39.5 101.2 25.9; -37.0 -9.2 15.0; 120.8 -6.8 13.8];
% Fcoord{22} = [0.0 96.2 -5.9; -74.7 -19.2 -60.2; 73.1 -15.1 -66.1]; %large change %[0.0 95.5 -2.3; -76.5 -20.0 -40.9; 75.4 -22.3 -33.0];
% Fcoord{23} = [2.1 92.8 -13.4; -75.4 -15.1 -52.7; 77.9 -11.7 -52.7]; %large change %[0.0 94.4 -9.1; -77.6 -11.4 -38.7; 79.9 -15.9 -35.3];
% Fcoord{24} = [0.4 94.5 -31.8; -69.7 -9.2 -46.0; 69.7 -10.0 -46.0]; %small change %[1.1 97.8 -31.8; -72.7 -15.9 -28.4; 71.6 -13.6 -35.3];
% Fcoord{25} = [1.3 96.2 -5.9; -67.1 -8.4 -46.8; 72.1 -10.9 -46.8]; %small change %[2.3 95.5 -1.1; -73.0 -12.5 -29.6; 75.3 -12.5 -25.0];
% Fcoord{27} = [-2.1 94.5 -19.2; -71.4 -15.9 -57.7; 73.9 -10.0 -59.4]; %[0.9 103.6 -3.1; -77.2 -11.2 -25.1; 77.2 -5.8 -25.1];
% Fcoord{28} = [0.0 89.1 -7.2; -67.8 -2.4 -67.2; 67.0 -2.0 -67.2]; %very large change %[-1.6 93.3 -9.9; -71.7 -20.2 -49.3; 72.8 -15.8 -47.1];
% Fcoord{29} = [-1.2 88.8 -8.1; -70.5 -16.9 -25.4; 69.0 -13.9 -25.4]; %[-0.6 86.7 -7.5; -74.3 -13.6 -15.8; 74.5 -13.6 -20.0];
% Fcoord{30} = [-4.6 92.0 -6.7; -73.1 -15.9 -27.6; 68.9 -6.7 -25.1]; %[-4.0 92.1 -8.0; -71.3 -1.7 -15.9; 71.3 -2.9 -14.8];
% Fcoord{31} = [0.0 89.5 -25.1; -64.7 -10.0 -60.2; 67.2 -5.9 -60.2]; %[0.6 92.1 -26.2; -72.4 -4.0 -23.9; 71.3 -8.6 -25.0];
% Fcoord{32} = [0.0 94.9 -6.1; -79.1 -25.3 -36.6; 77.3 -34.8 -37.5]; %[0.0 98.2 -8.3; -78.7 -22.5 -27.3; 78.7 -29.2 -21.3];
% Fcoord{33} = [1.8 96.3 -20.5; -62.6 -5.4 -46.4; 67.0 -7.1 -46.4]; %[-5 102 -10; -63 -14 -27; 72 -5 -17]; 
% Fcoord{34} = [0.0 86.1 -39.3; -70.4 -22.6 -49.3; 62.9 -14.2 -55.2]; %[0.0 88.7 -34.1; -73.0 -17.1 -29.6; 73.0 -18.2 -30.7];
% Fcoord{35} = [0.0 84.7 -10.7; -67.9 -17.8 -41.9; 63.5 -23.2 -40.1]; %[-0.6 89.8 -14.6; -70.9 -23.1 -29.1; 69.7 -26.7 -29.1];
% Fcoord{36} = [0.0 91.1 -21.7; -68.1 -12.5 -36.8; 63.0 -9.2 -41.8]; %[-0.6 93.3 -19.3; -70.2 -18.9 -25.0; 69.0 -12.0 -27.3];
% Fcoord{37} = [0.0 99.5 -11.7; -73.1 -10.0 -48.5; 75.6 -0.8 -48.5]; %[0.6 102.4 -10.2; -74.7 -7.4 -23.9; 78.1 -9.7 -22.7];
% Fcoord{39} = [0.0 114.2 4.5; -83.2 4.3 -2.4-25; 79.5 -13.3 2.4-25]; %[1.2 116.5 -6.1; -83.2 4.3 -2.4; 79.5 -13.3 2.4];
% Fcoord{40} = [-2.0 92.1 -36.2; -73.0 -28.0 -56.7; 75.8 -24.6 -52.6]; %[-1.4 98.8 -34.9; -75.6 -26.4 -36.8; 77.4 -26.4 -34.0];
% Fcoord{41} = [0.0 85.6 -25.9; -62.0 -8.0 -33.9; 59.3 -7.1 -33.9]; %[0.0 88.6 -23.1; -63.3 -10.9 -10.9; 65.7 -17.0 -19.4];
% Fcoord{42} = [4.5 96.3 8.0; -69.7 -13.4 -33.9; 66.1 -16.1 -35.7]; %[1.8 105.6 2.4; -74.5 -20.6 -14.6; 74.5 -21.8 -18.2];
% Fcoord{43} = [11.7 96.1 -39.2; -68.5 -10.2 -59.7; 71.4 -11.0 -61.4]; %[9.7 94.4 -32.7; -68.3 -9.4 -40.8; 75.0 -15.7 -43.1];
% Fcoord{44} = [0.1 74.0 -17.9; -66.7 -25.3 -29.0; 71.3 -19.8 -32.5];
% Fcoord{45} = [-1.6 117.3 -18.3; -70.3 0.9 -31.8; 70.1 9.2 -45.3];
% Fcoord{47} = [-4.7 98.3 -37.0; -67.6 -6.0 -30.2; 66.6 -4.8 -30.2]; 
% Fcoord{49} = [-0.2 99.4 -5.8; -62.3 7.8 -44.8; 71.2 7.8 -37.1];
% Fcoord{10} = [3.8 101.2 -38.4; -75.5 -6.1 -80.8; 69.5 -5.2 -86.2]; %[1.7 106.4 -36.3; -75.6 -9.4 -53.5; 74.1 -11.2 -53.5]; %LFCDB

%2) Subjects to run
OP.subjects = [2]; %49; %[49 47 46 42 4]; %[ 48];
% volt{49} = 2; %run a: volt = 2; run b: volt = 1
% volt{42} = 2; %run a: volt = 2; run b: volt = 1
OP.includeHR = 1; %include heart rate as a regressor
OP.LPFhrf = 1;
%Weak high pass filter
OP.weakHPF = 0;
if OP.weakHPF
    OP.HPForder = 5;
    OP.HPFfreq = 0.004;
else
    OP.HPForder = 4;
    OP.HPFfreq = 0.01;
end
OP.do_coreg = 1;
OP.force_coreg = 0;
OP.onlyCoregandPreprocess = 1;
OP.intrasubject_average = 0;
OP.stat_dir = 'Stat';
OP.main_dir = ['dataSPM' OP.run_label];
OP.coreg_dir = 'coreg';
%3) call main_module
OP.subj = subj;
OP.mux = mux;
OP.views = views;
OP.subj_fiducials = subj_fiducials;
% OP.Fcoord = Fcoord;
% OP.focus_size = focus_size;
OP.project_subject = project_subject;
OP.first_session = 0; %only run first session, and thus no group stats
OP.channel_pca = channel_pca;
OP.nComponents = nComponents;
OP.volt = volt;
OP.downsize = 1;
%all the above is the base case
OP.run = 0;
switch OP.run 
    case 1
        %description of the difference with the base case or the previous run
end
main_moduleEG(OP);
% OP.copy_figures = 0;
% %OP.path_coreg = 'coreg2\extra_coreg';
% OP.path_coreg = 'extra_coreg';
% copy_coreg(OP); 
%main_Projection(OP);
