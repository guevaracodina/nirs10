function epifMRI
%Batch for fMRI studies - for epilepsy
%INSTALLATION
%SPM8 needs to be present on the local machine
%It is available at http://www.fil.ion.ucl.ac.uk/spm/ and more precisely at
%http://www.fil.ion.ucl.ac.uk/spm/software/spm8/
%Take the LATEST version (currently July 2010)
%A path to spm8 need to be added in Matlab, for example by running a line
%such as: addpath('your_spm8_path'); where your_spm8_path is the location
%on your disk of your version of spm8. Careful: unzipping might give an
%spm8 folder within an spm8 folder; the outermost spm8 folder should be
%removed.

%USAGE:
%For each patient to be analyzed,
% 0- If you have MINC or DICOM files, they first need to be converted to
%    Nifti format. This can be done with SPM for MINC and MRIconvert for
%    DICOM.
% 1- create a folder with a name that begins with Raw_, such as Raw_epi127SD
%    The files in that folder will not be modified by the code; a new folder
%    with name a_, such as a_epi127SD, will be generated to put the analysis
% 2- In this Raw_ folder, put the Nifti .nii images (functional and anatomical)
%    and the files of markers that were generated with the "Export Markers"
%    function in Brain Vision Analyzer 2. This file of markers must contain
%    the epileptic spike markers and Scan Start markers. It is OK to have
%    both Scan Start markers and enumerated Scan Start markers.
%    Note that functionality of the script is limited when there is more
%    than one type of spikes. Optionally, cardiac peak markers can be
%    added (but does not work well if there are gaps in the ECG peaks)
%    Important: do not put additional files or folders in this folger, as
%    this will confuse the program which is trying to interpret all the
%    content of the folder
% 3- Make sure that the right functional file is associated with the right
%    markers file (also called "onsets" later). After launching the script
%    a message will be written in the Matlab command window and in the log
%    report detailing the association that the code made. This can easily
%    go wrong if the numbering of the files is inconsistent. For example, note
%    that files 10 and 11 will be listed before files 2, 3,... resulting
%    in a wrong order. For safety, relabel functional and onset files as
%    something like Session1, Session2, etc.
% 4- Set up parameters below - for recent acquisitions (with TR=3seconds)
%    no adjustment should be necessary.
% 5- Launch the script by pressing F5 in this window or the Green button
% 6- Select the Raw_ folders to be analyzed
% 7- Select folder where to put the analysis - the code should now run,
%    reporting on what it did and if errors were found
% 8- Look at the results:
%    a) Check the log file that the code ran correctly
%    b) Check the Ghostscript results reports, in a folder called
%    "ResultsAll"
%    c) Confirm results and produce diagnostic on each patient by taking
%    snapshots in SPM. This begins by clicking "Display Results" in SPM,
%    selecting SPM.mat in the desired folder of statistics for that patient
%    and creating a Powerpoint document with snapshots of regions of
%    (de)activation.

%Description of what the script does
% A- Check on the data and processing of onset files
% B- Preprocessing of the functional and anatomical images
%    1- Realign and unwarp (spatial)
%    2- Slice align (temporal) - perhaps best not to do it
%    3- Coregistration between functional and anatomical images
%    4- Normalize to Atlas - not required for individual epilepsy patients
%    5- Segmentation - not required
%    6- Smooth - sufficient smoothing is required for correct statistics
% C- Statistical analysis by General Linear Model (GLM)
%    1- Set up the GLM
%    2- Estimate the GLM
%    3- Generate contrasts
%    4- Produce result reports

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for user to set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Labels for statistic folders:
% a: slice timing correction done
% C: cardiac regressor included
% D: time and dispersion derivatives included in the GLM
% G: onsets were grouped and called "AllSpikes"
% u: unwrap was done (more precise realignment of functional images)
% V: Volterra (nonlinear analysis)
% w: normalization to Talairach Tournoux atlas was done
%
%To Do
% - MINC files
% - DICOM files
% - fix results reports
% - lags of movement parameters/economical movement parameters
% - convolve with cardiac HRF
% - Tyvaert analysis - need to generate F-stat triplets
% - check that normalize doesn't need to be called twice
% - save batch for later editing
% - analyse some of the sessions only
% - less copying of files? deleting files?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%When there are several types of stimuli, this option groups them into the same type
%The code will enforce grouping onsets if onset names are inconsistent
%between sessions. If onsets have been grouped, the stats folder(s) will have
%a "G" (for group) to help indicate that fact
Rest_no_onsets = 1; %Boolean to allow preprocessing only, no stats, for rest fMRI data
force_group_onsets = 1; %Boolean: 0: off; 1: on
allow_unequal_onset_numbers = 0; %Boolean: 0: off; 1: on
%onsets to remove
%Code allows up to 3 onset types to be removed - this is useful if onset
%files contains onsets that we want to exclude.
rem_onsets{1} = 'spk biF'; %'spkLFT';
rem_onsets{2} = 'eES';
rem_onsets{3} = 'sES';
%wait time in between each result report, adjustable so user can slow down
%the display if looking at the results during generation - set to zero for
%no delay
%waittime = 0.5; %in seconds
TR0 = 3.0; %13; %default TR value, in case code is unable to calculate it from
nslices = 47;
%the onsets files
%Run GLM without time and dispersion derivatives
noDerivs = 2; %Integer: 0: off (Derivs included); 1: on (noDerivs); 2: both

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advanced options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%skip slice timing correction - SPM recommends not doing slice timing
%correction - if slice timing correction is done, an "A" label will be
%added to the stats folder(s)
skip_slice_timing = 0; %Boolean: 0 (do not skip): off; 1: on (skip)

%Prefix put before name of patient folder, replacing 'Raw'
AnalysisPrefix = 'a';
%Create folder where functional files are put
EPIlabel = 'EPI';
%Created folder where EEG marker files are added
EEGlabel = 'EEG';
%Created folder where anatomical files are put
T1label = 'T1';
%Created folder where stat results are put
Statslabel = 'Stat_';
%results directory
DirResults = 'ResultsAll';
%epilepsy root name
RootName = 'epi';
DirLog = 'log';
%Run without graphical output - Unfortunately not implemented yet in SPM8
%Run_nogui = 1; %Boolean: 0: off; 1: on
%Analyze specific sessions - specify an array of sessions to analyze
%- if empty, code will analyze all sessions
Analyze_sessions = [];

SaveStatsBatch = 1; %Boolean: 0: off; 1: on
SaveReportBatch = 1; %Boolean: 0: off; 1: on
SkipUncorrected = 0; %Integer: 0: FWE & none: 1: FWE; 2: none
SkipGroupOnsets = 0; %Boolean: 0: (do not skip): off; 1: on (skip)
SkipNegativeBOLD = 0; %Boolean: 0: (do not skip): off; 1: on (skip)

%Prefixes for the various output files
pfx_temporal = 'a';
pfx_normalise = 'w';
pfx_realign = 'r';
pfx_unwrap = 'u';
pfx_smooth = 's';
Dlabel = {'','D'}; %No derivatives; Derivatives
Vlabel = {'','V'}; %No Volterra; Volterra
Wlabel = {'','W'}; %Not normalized; Normalized
Glabel = {'','G'}; %Canonical HRF (+derivs optional) or Gamma function(s)
%Size of smoothing kernel, in mm
skernel = 8;
skernelz = 8;
%Code will not run if there are no spikes in one of the sessions -
%By specifying a minimum number of spikes, sessions will be excluded
%and the code will run for the remaining sessions.
minimum_number_spikes = 0;
NewSegmentOn = 0; %Boolean: 0: off; 1: on
%normalize to TT atlas using Dartel - very long to run and not currently used
DartelOn = 0; %Boolean: 0: off; 1: on
%Normalize to Talairach-Tournoux atlas using old SPM normalise procedure
%Do not normalize to an atlas for a typical fMRI study in epilepsy
normalizeOn = 0; %Integer: 0: off; 1: on; 2: both
%During output of results to a Ghostscript file, results can be displayed
%by session
reportContrastsBySession = 0; %Boolean: 0: off; 1: on
%Generate stats with HRF peaking at different times with respect to onsets
McGilldelaysOn = 0; %Boolean: 0: off; 1: on
%Delays in seconds - to use with
delay = [-4 -2]; % -2 -1 1 3 6]; % [-9 -6 -3 -2 -1 1 2 3 6 9]; %[-7 -8 -5 -4]; %[-9 -6 -3 3 6 9];
%remove negative onsets entirely - negative onsets may arise when using
%a negative delay; it is best to remove such onsets
McGill_remove_negative = 1; %Boolean: 0: off; 1: on
%Generate stats with square or first derivative of movement parameters
MovementOn = 1; %1; %Boolean: 0: off; 1: on
%Adjust the onset to account that the middle slice starts TR/2 later
STC_middle_slice_adjust_onsets = 1; %Boolean: 0: off; 1: on
%Do slice timing correction centered on the middle slice
STC_middle_slice = 1; %Boolean: 0: off; 1: on

%Regenerate Result Reports even if stats have already been calculated
regenerate_reports = 1; %Boolean: 0: off; 1: on
%force running unwarp anyway even if movement is small
force_unwarp = 0; %Boolean: 0: off; 1: on
%threshold to unwarp in mm
unwarp_threshold = 0.3; %Boolean: 0: off; 1: on
%whether to include time derivative and dispersion of canonical HRF
%inc_derivs = 0;
add_pulse_regressor = 0; %Boolean: 0: off; 1: on
%Volterra - nonlinearities
VolterraOn = 0; %Integer: 0: off; 1: on; 2: both
%Window size for Gamma function HRF
GammaOn = 0; %Integer: 0: off; 1: on; 2: both
gamma_window = 20; %in seconds
gamma_order = 1; %number of gamma function bases
default_analysis_dir = 1; %Boolean: 0: let user choose analysis dir; 1: automatic default
AnovaOn = 0; %Boolean
threshold = [0.05 0.001]; %FWE, uncorrected thresholds
mask_threshold = 0.05;
extent = 1;
%Tyvaert Analysis - using groups of 2 second windows
TyvOn = 0; %Boolean: 0: off; 1: on
TyvWindow = 4; %in seconds %Careful, this should not be less than TR, to
%ensure that generated regressors are independent
%Time to begin constructing GLM bases before TyvStartSz
TyvBefore = 9; %in seconds
%Time to stop constructing GLM base after TyvEndSz
TyvAfter = 10; %in seconds
%Labels in onset files to look for
TyvStartSz = 'sES';
TyvEndSz = 'eES';
%number of regressors per group of F-stats
TyvPerGroup = 3;
TyvWriteFname = 1; %Boolean: 0: off; 1: on
TyvAdjustMultipleFtests = 1; %1;  %Boolean: 0: off; 1: on
TyvStatThreshold = 0.05;
TyvConjunctionMask = 0;  %Boolean: 0: off; 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    % Initialise SPM
    %--------------------------------------------------------------------------
    %spm('Defaults','fMRI');
    spm fmri;
    %spm_jobman('initcfg');
    %spm_check_installation;
    [t,sts] = spm_select([1,inf],'dir','Select folder(s) of data files to analyze (anatomical, functional and onsets)');
    if ~sts, return; end %fatal
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    disp(strvcat('SPM not found!',...
        'Solution: 1- Install SPM8 with latest update and',...
        '          2- Add spm8 to Matlab path'));
    return; %fatal
end

temp = which('spm');
sep = filesep;
idx = findstr(sep,temp);
DirSPM = temp(1:idx(end));

%number of subjects -- nominal (might be less if no good data found)
nomNsubj = size(t,1);
%actual number of subjects; initialize to zero
aNsubj = 0;

[dir2,dummy] = fileparts(t(1,:));
idx = findstr(sep,dir2);
temp_dir = dir2(1:idx(end));
if default_analysis_dir
    tA = temp_dir;
else
    try
        [tA,sts] = spm_select(1,'dir','Select folder to put all the analysis',{temp_dir});
        if ~sts
            try
                tA = spm_input('Enter analysis directory',1,'s',temp_dir);
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                disp('Problem with directory selection -- Aborting');
                return %fatal
            end
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        disp('Problem with directory selection -- Aborting');
        return %fatal
    end
end
if ~exist(tA,'dir'), mkdir(tA); end
DirAnalysis = tA;

DirStart = pwd; %keep track of initial directory
cd(DirAnalysis);

%Start time counter
tStart = tic;

%Check if log file exists
if ~exist(DirLog,'dir'), mkdir(DirLog); end
temp_log = fullfile(DirAnalysis,DirLog,['epifMRI_log' date]);
k = 0;
if exist(temp_log,'file')
    k = 2;
    %test = 1;
    while 1
        temp_log2 = [temp_log '_' int2str(k)];
        if exist(temp_log2,'file')
            k = k+1;
        else
            break;
        end
    end
end

if k
    log_file = temp_log2;
else
    log_file = temp_log;
end

%Open log file
try
    flog = fopen(log_file,'wt');
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    disp('Could not open log file in current directory -- exiting');
    return; %fatal
end

%Tyvaert
if TyvOn
    %Note FWHM of standard SPM8 gamma function response is 4.12 seconds and
    %it peaks at 3.01 seconds after the onset.
    tyv.TyvWindow = TyvWindow;
    tyv.TyvBefore = TyvBefore;
    tyv.TyvAfter = TyvAfter;
    tyv.TyvStartSz = TyvStartSz;
    tyv.TyvEndSz = TyvEndSz;
    tyv.TyvPerGroup = TyvPerGroup;
    tyv.TyvWriteFname = TyvWriteFname;
    tyv.TyvAdjustMultipleFtests = TyvAdjustMultipleFtests;
    tyv.TyvStatThreshold = TyvStatThreshold;
    tyv.TyvConjunctionMask = TyvConjunctionMask;
else
    tyv = [];
end
% if Run_nogui %Not implemented yet in SPM8
%     RunNoGui = '_nogui';
% else
%     RunNoGui = '';
% end
generate_patient_level = 1; %flag used to check if onsets allow generating patient level results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check user chosen data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Result: create three structures containing info per subject and per session:
%f{i}.fT1
%f{i}.fEPI{j}
%f{i}.fOnset{j}
%dStats{i} Boolean if GLM stats should be calculated
f{nomNsubj}.good_session = [];
for i=1:nomNsubj
    try %try to find data for this subject
        %Select all files recursively
        temp_dir = deblank(t(i,:));
        [filesRec,DirsRec{i}] = spm_select('FPListRec',temp_dir(1:end-1),'.*'); %need to remove a '\'
        
        %Initialize counters :
        % 1- nT1 of anatomical images (should be 1)
        % 2- nSess of functional images (number of sessions)
        % 3- nOnset of onset files (should match number of sessions)
        %There ought also to be a check that the onset files are in the correct
        %one-to-one correspondance with the functional images
        nT1 = 0;
        nSess = 0; %counter for EPI sessions
        nSess_onset = 0; %counter for onsets: final values should match with nSess
        %tentatively create a subject -- will remove it if unsuccessful
        aNsubj = aNsubj + 1;
        data_load_error = 0;
        %Look for files and assign them to specific directories
        
        disp(['Checking data... Subject ' int2str(i)]);
        f{i}.good_session = [];
        for j=1:size(filesRec,1)
            try
                %check if there are nifti images
                vols = spm_vol(filesRec{j,:}); %vols = vols{1,1};
                nImg = size(vols,1);
                if nImg == 1
                    %Anatomical image
                    nT1 = nT1 + 1;
                    if nT1 == 1
                        %add location of image for this subject
                        f{aNsubj}.fT1 = vols;
                    else
                        write_log(flog,'Too many images looking like anatomical images found.');
                        data_load_error = 1;
                    end
                else
                    if nImg > 1
                        %functional images found -- add a session
                        nSess = nSess + 1;
                        f{aNsubj}.fEPI{nSess} = vols;
                        if Rest_no_onsets
                            f{i}.good_session = [f{i}.good_session nSess];
                        end
                    else
                        %no good images
                        write_log(flog,'No good images. Code should not be coming here!');
                    end
                end %end if nImg == 1
                
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                %not images - then check they are onsets
                %two cases: raw files from Analyzer2 or .mat structure of
                %names, onsets, durations
                try
                    %check if this is a .mat file
                    load(filesRec{j,:},'names','onsets','durations');
                    %found a good onset file
                    nSess_onset = nSess_onset + 1;
                    too_few_spikes = 0;
                    if isempty(tyv)
                        for i1 =1:length(onsets)
                            if minimum_number_spikes > length(onsets{i1})
                                too_few_spikes = 1;
                            end
                        end
                    end
                    if ~too_few_spikes  || ~isempty(tyv) 
                        f{i}.good_session = [f{i}.good_session nSess_onset];
                        f{aNsubj}.fOnset{nSess_onset} = filesRec{j,:};
                    end
                    
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    try
                        %check there is not a .mat file already with the
                        %same name
                        onset_mat_file = 0; %assume onset .mat file does not already exist
                        for j2=1:size(filesRec,1)
                            if j2 ~= j
                                if strcmp(filesRec{j2,:},[filesRec{j,:} '.mat'])
                                    onset_mat_file = 1; %onset .mat file found
                                    write_log(flog,['2 onset files found; only '...
                                        [filesRec{j,:} '.mat'] ' will be used']);
                                    break;
                                end
                            end
                        end
                        %try to read the onsets as a textfile, extract also
                        %the repetition time of EPI images from Scan Start
                        %markers
                        if ~onset_mat_file
                            [TR,temp_file,too_few_spikes] = readOnsets(...
                                filesRec{j,:},rem_onsets,add_pulse_regressor,minimum_number_spikes,tyv,flog);
                            nSess_onset = nSess_onset + 1;
                            if ~too_few_spikes || ~isempty(tyv)
                                f{i}.good_session = [f{i}.good_session nSess_onset];
                                f{aNsubj}.fOnset{nSess_onset} = temp_file;
                            end
                            if TR < 0.5 || TR > 7
                                write_log(flog,['Warning: unusual TR value: ' num2str(TR)...
                                    ' for subject ' int2str(i)]);
                            end
                        end
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))%#ok<*CTCH>
                        write_log(flog,['Could not read onset file ' filesRec{j,:}]);
                    end
                end
            end
            
        end %end for j=1:size(filesRec,1)
        
        %check for fatal errors for this subject
        if data_load_error
            %remove this subject
            aNsubj = aNsubj-1;
            if aNsubj > 0
                f=f{1:aNsubj};
            else
                clear f; %if no subject left at this stage
            end
        else
            %check that there is at least a match between the number of EPI
            %sessions and onset files
            dStats{aNsubj} = 1; %data for GLM statistics
            if nSess_onset ~= nSess %Still do preprocessing, but don't do stats
                dStats{aNsubj} = 0;
                write_log(flog,strvcat(['For subject ' int2str(i)...
                    ': Mismatch between the number of EPI sessions (' int2str(nSess) ')'],...
                    ['and the number of onset files (' int2str(nSess_onset)...
                    ') -- only preprocessing will be done for this subject (no GLM statistics)']));
            end %end if nSess_onset ~= nSess
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))%try to find data for this subject
        write_log(flog,strvcat([t(i,:) ' is not a directory'],...
            ['This subject (' int2str(i) ') will be skipped']));
    end
end
%clean up
clear tA t vols filesRec data_load_error nT1 nSess nSess_onset nImg
clear names onsets durations onset_mat_file j2 nomNsubj sts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List data ready to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~aNsubj
    write_log(flog,'No subject found -- aborting');
    fclose(flog);
    return
end

write_log(flog,strvcat(['Found ' int2str(aNsubj) ' subject(s):']));
for i=1:aNsubj
    write_log(flog,strvcat(['Subject ' int2str(i) ': found ' int2str(size(f{i}.fEPI,2)) ' session(s),' ...
        ' with anatomical image: ' f{i}.fT1.fname]));
    if length(f{i}.fEPI) == length(f{i}.good_session)
        f{i}.RenameFiles = 0;
        write_log(flog,strvcat('of which all had enough onsets'));
    else
        f{i}.RenameFiles = 1;
        write_log(flog,strvcat(['of which only ' int2str(length(f{i}.good_session)) ...
            ' had enough onsets - sessions will be relabeled and renumbered']));
    end
    if dStats{i}
        for j=1:size(f{i}.fEPI,2)
            if any(j==f{i}.good_session)
                write_log(flog,strvcat(['Session ' int2str(j) ': ' ...
                    int2str(size(f{i}.fEPI{j},1)) ' volumes in ' ...
                    f{i}.fEPI{j}(1).fname ' matched with ' f{i}.fOnset{j}]));
            else
                write_log(flog,strvcat(['Session ' int2str(j) ': ' ...
                    'too few onsets']));
            end
        end
    else
        %only do preprocessing
        write_log(flog,'Problem with onsets or no onsets found -- only preprocessing will be done');
        for j=1:size(f{i}.fEPI,2)
            try
                f{i}.fOnset{j};
                write_log(flog,strvcat(['Session ' int2str(j) ': ' ...
                    int2str(size(f{i}.fEPI{j},1)) ' volumes in ' ...
                    f{i}.fEPI{j}(1).fname ' matched with ' f{i}.fOnset{j}]));
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,strvcat(['Session ' int2str(j) ': ' ...
                    int2str(size(f{i}.fEPI{j},1)) ' volumes in ' ...
                    f{i}.fEPI{j}(1).fname ' ; no onsets file']));
            end
        end
    end
end

write_log(flog,'Beginning data processing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create directory structure and copy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_log(flog,'Creating file structure and copying data files');
%working from analysis directory
cd(DirAnalysis);
DirResultsAll = [DirAnalysis sep DirResults sep];
if ~exist(DirResultsAll,'dir'), mkdir(DirResultsAll); end
try
    for i=1:aNsubj
        write_log(flog,['Begin processing Subject ' int2str(i)]);
        [dir2 filT1 ext1] = fileparts(f{i}.fT1.fname);
        idx = findstr(sep,dir2);
        
        %temp_dir = deblank(t(i,:));
        %[~,DirsRec] = spm_select('FPListRec',temp_dir(1:end-1),'.*');
        %try to find a sensible name for the subject root directory
        if isempty(DirsRec{i})
            tempDirSubj = dir2(idx(end)+1:end); % 'a_' for analysis -- could be changed
        else
            %Assume directory depth for this subject is one
            tempDirSubj = dir2(idx(end-1)+1:idx(end)-1);
        end
        idx = findstr(RootName,tempDirSubj);
        if isempty(idx)
            DirSubj = [AnalysisPrefix '_' tempDirSubj];
        else
            DirSubj = [AnalysisPrefix '_' tempDirSubj(idx:end)];
        end
        f{i}.DirSubj = DirSubj;
        if ~exist(DirSubj,'dir'), mkdir(DirSubj); end
        DirT1 = [DirSubj sep T1label sep];
        DirEPI = [DirSubj sep EPIlabel sep];
        DirEEG = [DirSubj sep EEGlabel sep];
        if ~exist(DirT1,'dir'), mkdir(DirT1); end
        if ~exist(DirEPI,'dir'), mkdir(DirEPI); end
        if ~exist(DirEEG,'dir'),
            mkdir(DirEEG);
        else
            %Clean up! Otherwise, might compound delays for onsets
            k = 1;
            while 1
                temp_dirEEG = [DirSubj sep EEGlabel '_' int2str(k)];
                if exist(temp_dirEEG,'file')
                    k = k+1;
                else
                    break;
                end
            end
            movefile(DirEEG,temp_dirEEG);
            %rmdir(DirEEG,'s');
            mkdir(DirEEG);
        end
        
        %Copy T1 file(s) if required
        temp_file = fullfile(DirT1,[ filT1 ext1]);
        if ~exist(temp_file,'file'),
            copyfile(f{i}.fT1.fname,temp_file);
            if strcmp(ext1,'.img')
                %copy the header file too
                copyfile([f{i}.fT1.fname(1:end-4) '.hdr'],fullfile(DirT1,[ filT1 '.hdr']));
            end
        end
        %Update the location of the file
        f{i}.fT1.fname = fullfile(DirAnalysis, temp_file);
        if ~f{i}.RenameFiles && isempty(Analyze_sessions)
            %Copy EPI and onset files if required -- looping over sessions
            for j=1:size(f{i}.fEPI,2)
                %check file name of first EPI file for this session
                [dummy, filEPI ext2] = fileparts(f{i}.fEPI{j}(1).fname);
                temp_file = fullfile(DirEPI,[filEPI ext2]);
                if ~exist(temp_file,'file'), %use spm_existfile instead?
                    if strcmp(ext2,'.img')
                        %loop over time
                        for k=1:size(f{i}.fEPI{j},1)
                            [dummy, filEPI ext2] = fileparts(f{i}.fEPI{j}(k).fname);
                            temp_file = fullfile(DirEPI,[filEPI ext2]);
                            copyfile(f{i}.fEPI{j}(k).fname,temp_file);
                            copyfile([f{i}.fEPI{j}(k).fname(1:end-4),fullfile(DirEPI,[filEPI '.hdr'])]);
                        end
                    else %.nii
                        copyfile(f{i}.fEPI{j}(1).fname,temp_file);
                    end
                end
                %Update the location of the file
                for k=1:size(f{i}.fEPI{j},1)
                    f{i}.fEPI{j}(k).fname = [fullfile(DirAnalysis,DirEPI,[filEPI ext2]) ',' int2str(k)];
                end
                
                %for onsets
                if dStats{i}
                    [dummy, filOnset ext2] = fileparts(f{i}.fOnset{j});
                    temp_file = fullfile(DirEEG,[filOnset ext2]);
                    if ~exist(temp_file,'file'),
                        copyfile(f{i}.fOnset{j},temp_file);
                    end
                    %Update the location of the file
                    f{i}.fOnset{j} = fullfile(DirAnalysis, temp_file);
                end
            end
        else
            %Create new directory for EPI files
            %To Do...
            
            DirEPI = [DirSubj sep EPIlabel sep];
            if ~exist(DirEPI,'dir'), mkdir(DirEPI); end
            nf = [];
            %T1
            nf.fT1.fname = fullfile(DirAnalysis, temp_file);
            %try
                nf.fEPI{length(f{i}.good_session)} = [];
                nf.fOnset{length(f{i}.good_session)} = [];
%             catch exception
%                 disp(exception.identifier)
%                 disp(exception.stack(1))
%             end
            %Copy and rename EPI and onset files if required -- looping over sessions
            n_sess = 0;
            for j=1:size(f{i}.fEPI,2)
                %only if this is a good session
                if (any(j==f{i}.good_session) && isempty(Analyze_sessions)) || ...
                        ( ~isempty(Analyze_sessions) && any(j==Analyze_sessions))
                    n_sess = n_sess+1;
                    
                    %check file name of first EPI file for this session
                    [dummy, dummy, ext2] = fileparts(f{i}.fEPI{j}(1).fname);
                    nfilEPI = ['Session' int2str(n_sess)]; %Rename
                    temp_file = fullfile(DirEPI,[nfilEPI ext2]);
                    if ~exist(temp_file,'file'), %use spm_existfile instead?
                        if strcmp(ext2,'.img')
                            %loop over time
                            for k=1:size(f{i}.fEPI{j},1)
                                [dummy, dummy, ext2] = fileparts(f{i}.fEPI{j}(k).fname);
                                temp_file = fullfile(DirEPI,[nfilEPI ext2]);
                                copyfile(f{i}.fEPI{j}(k).fname,temp_file);
                                copyfile([f{i}.fEPI{j}(k).fname(1:end-4),fullfile(DirEPI,[nfilEPI '.hdr'])]);
                            end
                        else %.nii
                            copyfile(f{i}.fEPI{j}(1).fname,temp_file);
                        end
                    end
                    %Update the location of the file
                    for k=1:size(f{i}.fEPI{j},1)
                        nf.fEPI{n_sess}(k).fname = [fullfile(DirAnalysis,DirEPI,[nfilEPI ext2]) ',' int2str(k)];
                    end
                    
                    %for onsets
                    if dStats{i} %must be true at this stage
                        [dummy, dummy, ext2] = fileparts(f{i}.fOnset{j});
                        nfilOnset = ['Session' int2str(n_sess)];
                        temp_file = fullfile(DirEEG,[nfilOnset ext2]);
                        if ~exist(temp_file,'file'),
                            copyfile(f{i}.fOnset{j},temp_file);
                        end
                        %Update the location of the file
                        nf.fOnset{n_sess} = fullfile(DirAnalysis, temp_file);
                        
                    end
                end
            end
            %replace f{i}
            f{i} = nf;
            f{i}.DirSubj = DirSubj;
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    write_log(flog,strvcat('Corrupted file structure in analysis directory',...
        ' -- Suggestion: remove content of analysis directory and start again'));
    return
end
clear DirEEG DirT1 DirEPI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data processing: big loop over subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:aNsubj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPM preprocessing: realign, unwarp, coreg, normalize, smooth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check if module has already been run
    [dummy, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
    cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
    
    if spm_existfile(['r' fil ext(1:end-2)])
        %Only check for first session
        write_log(flog,'Realign module already run -- skipping');
    else
        try
            %Generate file names required
            data = {};
            for j=1:size(f{i}.fEPI,2)
                %[dummy, fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                %fn = [DirAnalysis f{i}.DirSubj sep 'EPI' sep fil ext];
                data_sess = {};
                for k=1:size(f{i}.fEPI{j},1)
                    data_sess= [data_sess; f{i}.fEPI{j}(k).fname];
                end
                data= [data; {data_sess}];
            end
            clear matlabbatch
            matlabbatch{1}.spm.spatial.realign.estwrite.data = data';
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; %0: register to first image; 1: register to mean
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = pfx_realign;
            spm_jobman('run',matlabbatch);
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'Realign module failed');
        end
    end %end if spm_existfile(f{i}.fEPI{1}(1).fname)
    
    %Check movement parameters
    fMVT1 = spm_select('List',pwd,'^rp_');
    if size(fMVT1,1) == 1
        f{i}.fMVT1{1} = fMVT1;
    else
        for j=1:size(fMVT1,1)
            f{i}.fMVT1{j} = fMVT1(j,:);
        end
    end
    
    for j=1:size(f{i}.fMVT1,2)
        p1 = load(f{i}.fMVT1{j});
        %convert radians to degrees
        p1temp = [p1(:,1:3) p1(:,4:6)*180/pi];
        f{i}.unwarp = 0; %Boolean whether to use unwarp -- triggers to use it for all sessions
        if max(abs(p1temp(:))) > unwarp_threshold %set threshold to run unwarp and to include derivative of movement parameters
            f{i}.unwarp = 1;
        end
        dp1a = [zeros(1,6); diff(p1,1,1)];
        %Combine derivatives of motion parameters with motion parameters themselves
        dp1 = [p1 dp1a];
        f{i}.fMVT2{j} = ['rpd' f{i}.fMVT1{j}(3:end)];
        save(f{i}.fMVT2{j},'dp1','-ascii');
        %squares of parameters
        f{i}.fMVT22{j} = ['rp2' f{i}.fMVT1{j}(3:end)];
        p2a = p1 .* p1;
        p2 = [p1 p2a];
        save(f{i}.fMVT22{j},'p2','-ascii');
        %Combine derivatives of motion parameters with motion parameters
        %themselves and squares of parameters
        f{i}.fMVT222{j} = ['rpd2' f{i}.fMVT1{j}(3:end)];
        prd = [p1 dp1a p2a];
        save(f{i}.fMVT222{j},'prd','-ascii');
    end
    
    %Run unwarp if required
    if spm_existfile(['u' fil ext(1:end-2)])
        %Only check for first session
        write_log(flog,'Unwarp module already run -- skipping');
    else
        try
            if f{i}.unwarp || force_unwarp
                clear matlabbatch
                for j=1:size(f{i}.fEPI,2)
                    data_sess = {};
                    for k=1:size(f{i}.fEPI{j},1)
                        data_sess= [data_sess; f{i}.fEPI{j}(k).fname];
                    end
                    matlabbatch{1}.spm.spatial.realignunwarp.data(j).scans = data_sess;
                    matlabbatch{1}.spm.spatial.realignunwarp.data(j).pmscan = '';
                end
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0; %0: register to first image; 1: register to mean
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = {''};
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = pfx_unwrap;
                spm_jobman('run',matlabbatch);
            else
                write_log(flog,'Little movement on all sessions; no need to run unwrap');
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'Unwarp module failed');
        end
    end
    
    if f{i}.unwarp || force_unwarp
        Ulabel = pfx_unwrap;
        typ = pfx_unwrap;
    else
        Ulabel = '';
        typ = pfx_realign;
    end
    
    %New Segment
    if NewSegmentOn
        cd([DirAnalysis f{i}.DirSubj sep T1label]);
        [dummy, fil ext] = fileparts(f{i}.fT1.fname);
        if spm_existfile(['c1' fil ext(1:end)])
            write_log(flog,'New Segment already run -- skipping');
        else
            try
                clear matlabbatch
                matlabbatch{1}.spm.tools.preproc8.channel.vols = {f{i}.fT1.fname}; %{[f{i}.fT1.fname ',1']};
                matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
                matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
                matlabbatch{1}.spm.tools.preproc8.channel.write = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,1']};
                matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
                matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,2']};
                matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
                matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,3']};
                matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
                matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,4']};
                matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
                matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,5']};
                matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
                matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[DirSPM 'toolbox\Seg\TPM.nii,6']};
                matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
                matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [1 1];
                matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [1 1];
                matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
                matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
                matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
                matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
                spm_jobman('run',matlabbatch);
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'New Segment failed to run');
            end
        end
    end
    %Return to EPI directory
    cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
    
    try
        load('coreg_OK');
        %Coregistration already done
        write_log(flog,'Coregistration already done -- skipping');
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        %Coregistration
        %Select meanu or meanr image
        temp_fm = spm_select('List',pwd,'^meanu');
        if isempty(temp_fm)
            temp_fm = spm_select('List',pwd,'^mean');
            if isempty(temp_fm)
                write_log(flog,'Cannot find mean functional image - coregistration will fail');
            else
                %if there are several mean images, take the first
                if size(temp_fm,1) > 1
                    fm = {temp_fm(1,:)};
                else
                    fm = {temp_fm}; %{[temp_fm ',1']};
                end
            end
        else
            %if there are several mean images, take the first
            if size(temp_fm,1) > 1
                fm = {temp_fm(1,:)};
            else
                fm = {temp_fm};
            end
        end
        
        try
            %Generate file names required
            data = {};
            for j=1:size(f{i}.fEPI,2)
                [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                for k=1:size(f{i}.fEPI{j},1)
                    [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                    data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];
                end
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'File specification error at Coregistration step');
        end
        
        try
            clear matlabbatch
            matlabbatch{1}.spm.spatial.coreg.estimate.ref = {f{i}.fT1.fname};
            matlabbatch{1}.spm.spatial.coreg.estimate.source = fm;
            matlabbatch{1}.spm.spatial.coreg.estimate.other = data;
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            spm_jobman('run',matlabbatch);
            write_log(flog,'Coregistration module successful - writing coreg_OK file to avoid rerunning coreg');
            coreg_OK = 1;
            cd([DirAnalysis f{i}.DirSubj sep EPIlabel]);
            save('coreg_OK','coreg_OK');
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'Coregistration module failed to run');
        end
    end
    
    %Slice timing correction
    [dummy, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
    if spm_existfile([pfx_temporal typ fil ext(1:end-2)])
        write_log(flog,'Slice Timing Correction already run -- skipping');
    else
        if ~skip_slice_timing
            try
                %Generate file names required
                data = {};
                for j=1:size(f{i}.fEPI,2)
                    [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                    data_sess = {};
                    for k=1:size(f{i}.fEPI{j},1)
                        [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                        data_sess= [data_sess; [DirAnalysis f{i}.DirSubj sep EPIlabel sep typ fil ext]];
                    end
                    data= [data; {data_sess}];
                end
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'File specification error at Slice Timing Correction step');
            end
            try
                %load calculated TR
                TR = load(f{i}.fOnset{1},'TR');
                TR = TR.TR;
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,strvcat('No TR value found in first onsets file, or no onsets file found',...
                    'Using TR value of 3.013 for slice timing correction'));
                TR = TR0;
            end
            try
                clear matlabbatch
                matlabbatch{1}.spm.temporal.st.scans = data';
                matlabbatch{1}.spm.temporal.st.nslices = nslices;
                matlabbatch{1}.spm.temporal.st.tr = TR;
                matlabbatch{1}.spm.temporal.st.ta = TR*(1-1/nslices);
                matlabbatch{1}.spm.temporal.st.so = 1:nslices;
                if STC_middle_slice
                    matlabbatch{1}.spm.temporal.st.refslice = round(nslices/2);
                else
                    matlabbatch{1}.spm.temporal.st.refslice = 1;
                end
                matlabbatch{1}.spm.temporal.st.prefix = pfx_temporal;
                spm_jobman('run',matlabbatch);
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'Slice Timing Correction module failed to run');
            end
        end
    end
    if skip_slice_timing
        Alabel = '';
        Atyp = typ;
    else
        Alabel = pfx_temporal;
        Atyp = [pfx_temporal typ];
    end
    
    if normalizeOn
        if DartelOn
            %Return to T1 directory
            cd([DirAnalysis f{i}.DirSubj sep T1label]);
            [dummy, fil ext] = fileparts(f{i}.fT1.fname);
            %DARTEL
            if spm_existfile('Template_6.nii')
                write_log(flog,'DARTEL already run -- skipping');
            else
                try
                    %Generate file names required
                    data = {};
                    for j=1:size(f{i}.fEPI,2)
                        [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                        for k=1:size(f{i}.fEPI{j},1)
                            [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                            data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep Atyp fil ext]];
                        end
                    end
                    %Normalize some of the anatomical images
                    data = {[data; f{i}.fT1.fname ',1']};
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    write_log(flog,'File specification error for Dartel');
                end
                
                try
                    [dummy, fil ext] = fileparts(f{i}.fT1.fname);
                    data2 = {{[DirAnalysis f{i}.DirSubj sep T1label sep 'rc1' fil ext]}
                        {[DirAnalysis f{i}.DirSubj sep T1label sep 'rc2' fil ext]}};
                    clear matlabbatch
                    matlabbatch{1}.spm.tools.dartel.warp.images = data2';
                    matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
                    matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-006];
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
                    matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;
                    %Need to break into several modules and check if template
                    %was generated, as this takes a long time to run (almost an
                    %hour)
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1) = cfg_dep;
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).tname = 'DARTEL Template';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).tgt_spec{1}.name = 'filter';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).tgt_spec{1}.value = 'nifti';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).sname = 'Run DARTEL (create Templates): Template (Iteration 6)';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                    matlabbatch{2}.spm.tools.dartel.mni_norm.template(1).src_output = substruct('.','template', '()',{7});
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1) = cfg_dep;
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).tname = 'Flow fields';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).tgt_spec{1}.name = 'filter';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).tgt_spec{1}.value = 'nifti';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).sname = 'Run DARTEL (create Templates): Flow Fields';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1).src_output = substruct('.','files', '()',{':'});
                    matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.images = data';
                    matlabbatch{2}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
                    matlabbatch{2}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN
                        NaN NaN NaN];
                    matlabbatch{2}.spm.tools.dartel.mni_norm.preserve = 0;
                    matlabbatch{2}.spm.tools.dartel.mni_norm.fwhm = [skernel skernel skernelz];
                    spm_jobman('run',matlabbatch);
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    %Run Normalise instead
                    try
                        %Generate file names required
                        data = {};
                        for j=1:size(f{i}.fEPI,2)
                            [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                            for k=1:size(f{i}.fEPI{j},1)
                                [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                                data = [data; [DirAnalysis f{i}.DirSubj sep EPIlabel sep Atyp fil ext]];
                            end
                        end
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))
                        write_log(flog,'File specification error for Normalise');
                    end
                    %Normalize some of the anatomical images
                    data = {[data; f{i}.fT1.fname]};
                    call_normalize(data,f{i}.fT1.fname,dirSPM,flog,pfx_normalise);
                end
            end
        else %if DartelOn
            %Run Normalise instead
            [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
            tlabel = [pfx_normalise Atyp];
            if spm_existfile([tlabel fil ext(1:end-2)])
                write_log(flog,'Normalise already run -- skipping');
            else
                try
                    %Generate file names required
                    data = {};
                    for j=1:size(f{i}.fEPI,2)
                        %[dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                        for k=1:size(f{i}.fEPI{j},1)
                            [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                            data = [data; [DirAnalysis f{i}.DirSubj sep EPIlabel sep Atyp fil ext]];
                        end
                    end
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    write_log(flog,'File specification error for Normalise');
                end
                data = {[data; f{i}.fT1.fname]};
                call_normalize(data,f{i}.fT1.fname,dirSPM,flog,pfx_normalise);
                %Normalize anatomical T1 image
                %data = {f{i}.fT1.fname};
                %call_normalize(data,f{i}.fT1.fname,dirSPM,flog);
            end
        end %end if DartelOn
    end %end if normalizeOn
    
    
    %Smooth: several runs
    %1 - Smooth ar or au files (or directly r or u files)
    %2 - Smooth war or wau files
    %Return to EPI directory
    cd([DirAnalysis f{i}.DirSubj sep EPIlabel]);
    [dummy, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
    if spm_existfile([pfx_smooth Atyp fil ext(1:end-2)])
        write_log(flog,'Smooth already run -- skipping');
    else
        try
            %Generate file names required
            data = {};
            for j=1:size(f{i}.fEPI,2)
                %[dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                for k=1:size(f{i}.fEPI{j},1)
                    [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                    data = [data; [DirAnalysis f{i}.DirSubj sep EPIlabel sep Atyp fil ext]];
                end
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'File specification error at Smooth step');
        end
        
        try
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = data;
            matlabbatch{1}.spm.spatial.smooth.fwhm = [skernel skernel skernelz];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = pfx_smooth;
            spm_jobman('run',matlabbatch);
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            write_log(flog,'Smooth module on ar/au files failed to run');
        end
    end
    
    %Smooth #2 (need to call a function rather than repeat code here
    if normalizeOn
        Wtyp = [pfx_normalise Atyp];
        cd([DirAnalysis f{i}.DirSubj sep EPIlabel]);
        [dummy, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
        if spm_existfile([pfx_smooth Wtyp fil ext(1:end-2)])
            write_log(flog,'Smooth already run -- skipping');
        else
            try
                %Generate file names required
                data = {};
                for j=1:size(f{i}.fEPI,2)
                    [dummy,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                    for k=1:size(f{i}.fEPI{j},1)
                        [dummy,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                        data = [data; [DirAnalysis f{i}.DirSubj sep EPIlabel sep Wtyp fil ext]];
                    end
                end
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'File specification error at Smooth step');
            end
            
            try
                clear matlabbatch
                matlabbatch{1}.spm.spatial.smooth.data = data;
                matlabbatch{1}.spm.spatial.smooth.fwhm = [skernel skernel skernelz];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = pfx_smooth;
                spm_jobman('run',matlabbatch);
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'Smooth failed to run');
            end
        end
    end
    
    %end of preprocessing: copy ghostview files
    try
        cd([DirAnalysis f{i}.DirSubj sep EPIlabel]);
        temp_file = spm_select('List',pwd,'^spm_');
        for nf=1:size(temp_file,1)
            [dummy, fil ext] = fileparts(temp_file(nf,:));
            if exist(temp_file(nf,:),'file')
                copyfile(temp_file(nf,:),fullfile(DirResultsAll,[f{i}.DirSubj '_PreProcess_' fil ext]));
            end
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        write_log(flog,'Directory error: could not copy preprocessing .gs file');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPM statistics: model specification, model estimation, contrats, results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Adjust onsets
    try
        if STC_middle_slice_adjust_onsets && STC_middle_slice
            for j=1:size(f{i}.fEPI,2)
                load(f{i}.fOnset{j});
                M = size(names,2);
                try
                    TR;
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    TR = TR0;
                    write_log(flog,'TR was not found in onset file at onset adjust step');
                end
                for m=1:M
                    onsets{m} = onsets{m} - TR/2;
                    if McGill_remove_negative
                        p=0;
                        clear temp_onsets temp_durations
                        for k=1:size(onsets{m},2)
                            if onsets{m}(k) >= 0 % keep only positive onsets
                                p = p + 1;
                                temp_onsets{m}(p) = onsets{m}(k);
                                temp_durations{m}(p) = durations{m}(k);
                            end
                        end
                        onsets{m} = temp_onsets{m};
                        durations{m} = temp_durations{m};
                    else
                        for k=1:size(onsets{m},2)
                            %catch potentially negative onsets
                            if onsets{m}(k) < 0, onsets{m}(k) = 0; end
                        end
                    end
                end
                save(f{i}.fOnset{j},'names','onsets','durations','TR');
            end
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        write_log(flog,'Could not adjust onsets for STC_middle_slide');
    end
    
    %Check if there are the same number and types of onsets in each onset file
    try
        permute_onsets = 0;
        %cd([DirAnalysis DirSubj]);
        DirOtherOnsets = fullfile(DirAnalysis, f{i}.DirSubj, EEGlabel,'OtherOnsets');
        if ~exist(DirOtherOnsets,'dir'), mkdir(DirOtherOnsets); end
        %Onset file, Movement
        load(f{i}.fOnset{1});
        names2 = names;
        f{i}.GroupOnsets = 0;
        for j=2:size(f{i}.fEPI,2)
            load(f{i}.fOnset{j});
            M = size(names,2);
            if size(names2,2) == M;
                for m=1:M
                    if ~strcmpi(names2{m},names{m})
                        %try
                        %Try to permute the onsets
                        permute_onsets = 1;
                        %                     catch
                        %                         write_log(flog,strvcat(['Subject ' int2str(i) ': Same number but different types of onsets in session ' int2str(j)],...
                        %                             'Onsets will be grouped into one of the same type'));
                        %                         f{i}.GroupOnsets = 1;
                        %                     end
                    end
                end
            else
                write_log(flog,strvcat(['Subject ' int2str(i) ': Different number of onsets in session ' int2str(j)]));
                if ~allow_unequal_onset_numbers || force_group_onsets
                    write_log(flog,'Onsets will be grouped into the same type');
                    f{i}.GroupOnsets = 1;
                else
                    write_log(flog,'Only session-level results will be generated');
                    generate_patient_level = 0;
                end
                
                %Don't bother trying to permute onsets if number of onsets
                %doesn't match
                permute_onsets = 0;
            end
        end
        
        if permute_onsets && ~f{i}.GroupOnsets
            %Onset file, Movement
            load(f{i}.fOnset{1});
            names2 = names;
            clear onsets durations
            %generate all permutations
            perm1 = perms(1:size(names2,2));
            %loop over sessions
            for j=2:size(f{i}.fEPI,2)
                load(f{i}.fOnset{j});
                %loop over permutations
                for k=1:size(perm1,1)
                    %start assuming this is the good permutation
                    good_perm = 1;
                    %loop over onset types
                    for m=1:size(names2,2)
                        %find the permutation for which all the onset names
                        %match the ones of the first session
                        if ~strcmpi(names{perm1(k,m)},names2{m})
                            good_perm = 0; %not a good permutation
                        end
                    end
                    if good_perm == 1, break; end
                end
                %If we get here, that's because for this value of k,
                %there was no mismatch of names, so this is the correct
                %permutation
                names = names2;
                for m=1:size(names2,2)
                    temp_onsets{m} = onsets{perm1(k,m)};
                    temp_durations{m} = durations{perm1(k,m)};
                end
                onsets = temp_onsets;
                durations = temp_durations;
                %save over the old onsets, don't bother to move them
                save(f{i}.fOnset{j},'names','onsets','durations','TR');
            end
            
        end
        %force group onsets
        if force_group_onsets, f{i}.GroupOnsets = 1; end
        
        %Create new onset files
        if f{i}.GroupOnsets
            for j=1:size(f{i}.fEPI,2)
                load(f{i}.fOnset{j});
                names = {'AllSpikes'};
                if size(onsets,2) > 1
                    c = onsets{1};
                    d = durations{1};
                    for m=2:size(onsets,2)
                        [c ix] = sort([c onsets{m}]);
                        d = [d durations{m}];
                        d = d(ix); %permute the elements
                    end
                else
                    c = onsets{1};
                    d = durations{1};
                end
                clear onsets durations
                onsets{1} = c;
                durations{1} = d;
                [dummy,fil ext] = fileparts(f{i}.fOnset{j});
                save(fullfile(DirOtherOnsets,[fil '_grouped' ext]),'names','onsets','durations','TR');
            end
        end
        
        %Move onset files that will not be used
        if f{i}.GroupOnsets
            DirNotUsed = fullfile(DirAnalysis,f{i}.DirSubj,EEGlabel,'NotUsed');
            if ~exist(DirNotUsed,'dir'), mkdir(DirNotUsed); end
            for j=1:size(f{i}.fEPI,2)
                [dummy,fil ext] = fileparts(f{i}.fOnset{j});
                movefile(f{i}.fOnset{j},[DirNotUsed '\' fil ext]);
                f{i}.fOnset{j} = fullfile(DirOtherOnsets,[fil '_grouped' ext]);
            end
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        write_log(flog,'Could not group onset files');
    end
    
    if f{i}.GroupOnsets
        glabel = 'g';
    else
        glabel = '';
    end
    
    %Add pulse regressor to rp_files
    if add_pulse_regressor
        Clabel = 'C';
        try
            for j=1:size(f{i}.fMVT1,2)
                %add to movement parameters
                p1 = load(f{i}.fMVT1{j});
                load(f{i}.fOnset{j});
                p1 = [p1 vRi'];
                [dir1 fil1 ext1] = fileparts(f{i}.fMVT1{j});
                new_rp = fullfile(dir1,['c' fil1 ext1]);
                save(new_rp,'p1','-ascii');
                f{i}.fMVT1{j} = new_rp;
                %add to square of movement parameters
                p1 = load(f{i}.fMVT22{j});
                p1 = [p1 vRi'];
                [dir1 fil1 ext1] = fileparts(f{i}.fMVT22{j});
                new_rp = fullfile(dir1,['c' fil1 ext1]);
                save(new_rp,'p1','-ascii');
                f{i}.fMVT22{j} = new_rp;
                %add to derivative of movement parameters
                p1 = load(f{i}.fMVT2{j});
                p1 = [p1 vRi'];
                [dir1 fil1 ext1] = fileparts(f{i}.fMVT2{j});
                new_rp = fullfile(dir1,['c' fil1 ext1]);
                save(new_rp,'p1','-ascii');
                f{i}.fMVT2{j} = new_rp;
                %add both derivative and square of movement parameters
                p1 = load(f{i}.fMVT222{j});
                p1 = [p1 vRi'];
                [dir1 fil1 ext1] = fileparts(f{i}.fMVT222{j});
                new_rp = fullfile(dir1,['c' fil1 ext1]);
                save(new_rp,'p1','-ascii');
                f{i}.fMVT222{j} = new_rp;
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
        end
    else
        Clabel = '';
    end
    
    GLM = [];
    GLM.reportContrastsBySession = reportContrastsBySession;
    GLM.flog = flog;
    GLM.DirAnalysis = DirAnalysis;
    GLM.DirResultsAll = DirResultsAll;
    %GLM.waittime = waittime;
    GLM.f = f{i};
    GLM.dStats = dStats{i};
    GLM.regenerate_reports = regenerate_reports;
    GLM.TR0 = TR0;
    GLM.Statslabel = Statslabel;
    GLM.Atyp = Atyp;
    GLM.Alabel = Alabel; %careful, not Atyp
    GLM.Clabel = Clabel;
    GLM.Dlabel = Dlabel{1};
    GLM.Wlabel = Wlabel{1};
    GLM.Ulabel = Ulabel;
    GLM.Vlabel = Vlabel{1};
    GLM.Mlabel = '';
    GLM.Plabel = '';
    GLM.glabel = glabel;
    GLM.Glabel = Glabel{1};
    GLM.Slabel = pfx_smooth;
    GLM.AnovaOn = AnovaOn;
    GLM.extent = extent;
    GLM.threshold = threshold;
    GLM.mask_threshold = mask_threshold;
    GLM.gamma_window = gamma_window;
    GLM.gamma_order = gamma_order;
    %GLM.TyvOn = TyvOn;
    GLM.tyv = tyv;
    GLM.SaveStatsBatch = SaveStatsBatch;
    GLM.SaveReportBatch = SaveReportBatch;
    GLM.SkipUncorrected = SkipUncorrected;
    GLM.SkipGroupOnsets = SkipGroupOnsets;
    GLM.SkipNegativeBOLD = SkipNegativeBOLD;
    GLM.generate_patient_level = generate_patient_level;
    if ~generate_patient_level
        GLM.reportContrastsBySession = 1; %otherwise there is nothing left to output!
    end
    switch noDerivs
        case 0
            DerivsRun = {2};
        case 1
            DerivsRun = {1};
        case 2
            DerivsRun = {1,2};
    end
    switch VolterraOn
        case 0
            VolterraRun = {1};
        case 1
            VolterraRun = {2};
        case 2
            VolterraRun = {1,2};
    end
    switch normalizeOn
        case 0
            normalizeRun = {1};
        case 1
            normalizeRun = {2};
        case 2
            normalizeRun = {1,2};
    end
    switch GammaOn
        case 0
            GammaRun = {1};
        case 1
            GammaRun = {2};
        case 2
            GammaRun = {1,2};
    end
    
    if ~TyvOn
        for i1=1:length(DerivsRun)
            GLM.Dlabel = Dlabel{DerivsRun{i1}};
            for i2=1:length(VolterraRun)
                GLM.Vlabel = Vlabel{VolterraRun{i2}};
                for i3=1:length(normalizeRun)
                    GLM.Wlabel = Wlabel{normalizeRun{i3}};
                    for i4=1:length(GammaRun)
                        GLM.f = f{i};
                        GLM.Glabel = Glabel{GammaRun{i4}};
                        run_stats_and_results(GLM);
                        
                        try
                            if McGilldelaysOn && i4 == 1
                                %GLM.Glabel = Glabel{2};
                                %delays
                                for d=1:size(delay,2)
                                    g = f{i};
                                    for j=1:size(f{i}.fEPI,2)
                                        load(f{i}.fOnset{j});
                                        M = size(names,2);
                                        for m=1:M
                                            onsets{m} = onsets{m} + delay(d); %subtract delay
                                            if McGill_remove_negative
                                                p=0;
                                                clear temp_onsets temp_durations
                                                for k=1:size(onsets{m},2)
                                                    if onsets{m}(k) >= 0 % keep only positive onsets
                                                        p = p + 1;
                                                        temp_onsets{m}(p) = onsets{m}(k);
                                                        temp_durations{m}(p) = durations{m}(k);
                                                    end
                                                end
                                                onsets{m} = temp_onsets{m};
                                                durations{m} = temp_durations{m};
                                            else
                                                for k=1:size(onsets{m},2)
                                                    %catch potentially negative onsets
                                                    if onsets{m}(k) < 0, onsets{m}(k) = 0; end
                                                end
                                            end
                                        end
                                        [dummy,fil ext] = fileparts(f{i}.fOnset{j});
                                        g.fOnset{j} = fullfile(DirOtherOnsets,[ fil '_' int2str(delay(d)) 's' ext]);
                                        save(g.fOnset{j},'names','onsets','durations','TR');
                                    end
                                    GLM.f = g;
                                    GLM.Mlabel = [int2str(delay(d)) 's delay'];
                                    %Statistics with onset delays, to move peak of HRF
                                    run_stats_and_results(GLM);
                                    
                                end %end for d=1:size(delays,2)
                            end
                        catch exception
                            disp(exception.identifier);
                            disp(exception.stack(1));
                            write_log(flog,'McGill delays stats failed to run');
                        end
                        GLM.Mlabel = '';
                        
                        %Movement
                        try
                            if MovementOn
                                %1st Derivative of movement parameters
                                g = f{i};
                                for j=1:size(f{i}.fMVT1,2)
                                    g.fMVT1{j} = f{i}.fMVT2{j};
                                end
                                GLM.Plabel = 'DerMVT';
                                GLM.f = g;
                                run_stats_and_results(GLM);
                                
                                %square of movement parameters
                                g = f{i};
                                for j=1:size(f{i}.fMVT1,2)
                                    g.fMVT1{j} = f{i}.fMVT22{j};
                                end
                                GLM.f = g;
                                GLM.Plabel = 'sqrMVT';
                                run_stats_and_results(GLM);
                                
                                %square of movement parameters and
                                %1st Derivative of movement parameters
                                g = f{i};
                                for j=1:size(f{i}.fMVT1,2)
                                    g.fMVT1{j} = f{i}.fMVT222{j};
                                end
                                GLM.f = g;
                                GLM.Plabel = 'DersqrMVT';
                                run_stats_and_results(GLM);
                            end
                        catch exception
                            disp(exception.identifier)
                            disp(exception.stack(1))
                            write_log(flog,'Movement stats failed to run');
                        end
                        GLM.f = f{i};
                        GLM.Plabel = '';
                    end
                end
            end
        end
        
        
        
    else
        %Tyvaert analysis
        GLM.Glabel = Glabel{2}; %Gamma function
        GLM.Dlabel = ''; %No derivatives
        GLM.Vlabel = ''; %No Volterra
        GLM.SkipGroupOnsets  = 1; %Do not group onsets
        run_stats_and_results(GLM);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data processing: end big loop over subjects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end for i=1:aNsubj

%return to initial directory
cd(DirStart);

%end time counter
tElapsed = toc(tStart);

%Run was successful
write_log(flog,strvcat('Run completed',['Time taken: ' num2str(tElapsed/3600) ' hours'])); %#ok<*VCAT>
%Close log file
fclose(flog);
end

function write_log(flog,mlog)
%write to screen and to log file
disp(mlog);
for i=1:size(mlog,1);
    fprintf(flog,'%s\n',mlog(i,:));
end
end

function [TR,fOnset,too_few_spikes] = readOnsets(...
    file,rem_onsets,add_pulse_regressor,minimum_number_spikes, tyv,flog)
%Reads Markers file exported by Analyzer, to extract 'multiple conditions'
%data required by SPM: names, onsets, durations
%Returns TR, the EPI repetition time, and the name of the new onset file

fp = fopen(file,'r');

%Extract frequency information from first line
temp = fgetl(fp);
i1 = findstr(':',temp); i2 = findstr('ms',temp);
%Delta t between points in seconds
dt = str2double(temp(i1(2)+1:i2-1))/1000;

%find first Scan Start -- ignore onsets that would arrive before the
%first Scan Start
while ~feof(fp)
    temp = fgetl(fp);
    indx = findstr('Scan Start',temp);
    if ~isempty(indx)
        indx = findstr(',',temp);
        start_scan = dt*str2double(temp(indx(2):indx(3)));
        break
    end
end

%Boolean to keep track if TR calculation has been done
TRdone = 0;
%Counters for pulse and movement markers
iR = 0;
iMVT = 0;
%vectors for pulse and movement regressors
vR = [];
vMVT = [];

%Number of types of onsets
NTonset = 0;
names = {};
onsets = {};
durations = {};
%counter of onsets
k = {};

while ~feof(fp)
    temp = fgetl(fp);
    if ~TRdone
        %find next Scan Start
        indx = findstr('Scan Start',temp);
        if ~isempty(indx)
            indx = findstr(',',temp);
            next_start_scan = dt*str2double(temp(indx(2):indx(3)));
            TR_temp = next_start_scan - start_scan;
            if TR_temp > 0
                TR = TR_temp;
                TRdone = 1;
            end
        end
    end
    indx = findstr(',',temp);
    tag = strtrim(temp(indx(1)+1:indx(2)-1));
    if ((~strcmpi(tag,'Scan Start') && ~strcmpi(tag,'mvt') && ~strncmp(tag,'R',1) ...
            && ~strcmpi(tag,'Bad Interval') && ~strncmp(tag,'B',1) ...
            && (~strncmp(tag,'S',1) || strncmpi(tag,'Spk',3) || strncmpi(tag,'SWD',3)) && ~strncmp(tag,'T',1) ...
            && (~strcmp(tag,rem_onsets{1})) ...
            && (~strcmp(tag,rem_onsets{2})) ...
            && (~strcmp(tag,rem_onsets{3}))) && isempty(tyv)) ...
            || (~isempty(tyv) && (strcmpi(tag,tyv.TyvStartSz) || strcmpi(tag,tyv.TyvEndSz)))
        %Possible new onset
        not_found = 1; %true
        for i=1:NTonset
            %check whether onset is one of those already found
            if strcmpi(names{i},tag) %case insensitive
                if isempty(tyv)
                    not_found = 0; %found
                    k{i} = k{i}+1; %#ok<*AGROW>
                    onsets{i}(k{i}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                    durations{i}(k{i}) = dt*(str2double(temp(indx(3):indx(4)))-1); %subtract 1 so point onsets have zero duration
                else
                    if strcmpi(tag,tyv.TyvStartSz)
                        not_found = 0; %found
                        k{i} = k{i}+1; %#ok<*AGROW>
                        onsets{i}(k{i}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                        durations{i}(k{i}) = 0;
                    end
                end
                break
            end
        end
        if not_found && (isempty(tyv) || (~isempty(tyv) && (strcmpi(tag,tyv.TyvStartSz) || strcmpi(tag,tyv.TyvEndSz))))
            if isempty(tyv)
                %new type of onset
                NTonset = NTonset + 1;
                k{NTonset} = 1;
                names{NTonset} = tag;
                onsets{NTonset}(k{NTonset}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                durations{NTonset}(k{NTonset}) = dt*(str2double(temp(indx(3):indx(4)))-1); %subtract 1 so point onsets have zero duration
            else
                if strcmpi(tag,tyv.TyvStartSz)
                    %new type of onset
                    NTonset = NTonset + 1;
                    k{NTonset} = 1;
                    names{NTonset} = tag;
                    onsets{NTonset}(k{NTonset}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                    durations{NTonset}(k{NTonset}) = 0;
                else
                    if strcmpi(tag,tyv.TyvEndSz)
                        %found an end seizure marker, so the duration is
                        %its position less the previous start seizure
                        %marker position
                        durations{NTonset}(k{NTonset}) = str2double(temp(indx(2):indx(3)))*dt-start_scan-onsets{NTonset}(k{NTonset});
                    end
                end
            end
        end
    else
        %ignore pulses, etc., before first Scan Start
        %Look for movement or pulse markers
        if strncmp(tag,'R',1) && length(tag) == 1 %pulse
            iR = iR + 1;
            %Time of pulse
            vR(iR) =  str2double(temp(indx(2):indx(3)))*dt-start_scan;
            
        end
        if strncmp(tag,'mvt',1) %movement
            iMVT = iMVT + 1;
            %Start time of movement
            vMVT(iMVT) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
            %Durations
            dMVT(iMVT) = dt*(str2double(temp(indx(3):indx(4)))-1);
        end
        if strncmp(tag,'Scan Start',9)
            %find last volume
            idx = findstr(tag,'Scan Start');
            tmp_last = tag(idx+10:end);
            if ~isempty(tmp_last)
                last_volume = tmp_last;
            end
        end
    end
end

fclose(fp);

try
    %Calculate pulse regressor
    
    %Get number of volumes
    last_volume = str2double(last_volume);
    lpi = linspace(0,last_volume*TR, last_volume*TR/dt); %*res_factor);
    
    %Complete missing pulses
    dvR = diff(vR);
    %median pulse rate
    medPR = median(dvR);
    var_allowed = 0.5;
    %List of intervals with missing pulses
    [dummy,idx] = find(dvR > medPR*(1+var_allowed)); %??
    %dvRbad = dvR(idx);
    for i1=1:length(idx)
        if idx(i1) > 1 && idx(i1) < length(dvR)
            %new pulses
            tmp_PR = (dvR(idx(i1)-1)+dvR(idx(i1)+1))/2;
        else
            if idx(i1) == 1
                tmp_PR = dvR(idx(i1)+1);
            else
                if idx(i1) == length(dvR)
                    tmp_PR = dvR(idx(i1)-1);
                end
            end
        end
        %Round to nearest integer to get good approximation to number
        %of missing pulses
        n_P = round(dvR(idx(i1))/tmp_PR);
        %start and end values of vR
        int1 = dvR(idx(i1))/n_P;
        s1 = vR(idx(i1));
        %e1 = vR(idx(i1)+1);
        add_P{i1} = linspace(s1+int1,s1+(n_P-1)*int1,n_P-1); %??
        %tmp_vR = [vR(1:1)]; %??
    end
    
    %new pulse rate series
    if ~isempty(idx)
        %Distinguish if first gap occurs at start of dvR or in generic
        %place
        if idx(1) == 1
            tmp_vR = add_P{1};
        else
            tmp_vR = [vR(1:idx(1)) add_P{1}];
        end
        %generic gaps
        for i1=2:length(idx)
            tmp_vR = [tmp_vR vR(idx(i1-1)+1:idx(i1)) add_P{i1}]; %??
        end
        %Distinguish if last gap occurs at end of dvR or in generic place
        if idx(end) == length(dvR)
            %nothing to do
        else
            tmp_vR = [tmp_vR vR(idx(end)+1:end)];
        end
    end
    vR = tmp_vR;
    dvR = diff(vR);
    vRi = interp1(vR(1:end),[dvR dvR(end)],(1:last_volume)*TR,'nearest');
    %derivative of pulse rate - too noisy
    %ddvR = diff(dvR);
    %vRi2 = interp1(vR(2:end-1),ddvR,lpi,'spline');
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end

%Tyvaert
if ~isempty(tyv)
    try
        OldOnsets = onsets{1};
        OldDurations = durations{1};
        clear names onsets durations
        %number of seizures
        nSz = length(OldOnsets);
        if nSz > 1
            write_log(flog,'Tyvaert analysis: make sure there is only one seizure per session');
        end
        %nb(nSz) = 0;
        %for i1=1:nSz
        %number of basis functions for seizure 1
        nb = round((OldDurations+ min(OldOnsets,tyv.TyvBefore)+tyv.TyvAfter)/tyv.TyvWindow);
        for i2=1:nb
            if i2<10
                str = '0';
            else
                str = '';
            end
            if tyv.TyvWriteFname %Note that this gives incorrect name for T contrasts, but
                %correct name for F contrast
                ts = round((i2-1)*tyv.TyvWindow-tyv.TyvBefore);
                te = round((i2+tyv.TyvPerGroup-1)*tyv.TyvWindow-tyv.TyvBefore);
                Tte = round(i2*tyv.TyvWindow-tyv.TyvBefore);
                
                %                 if i2 < nb-round(tyv.TyvAfter/tyv.TyvWindow)
                %                     tstr = '';
                %                 else
                %                     tstr = '';
                %                 end
                names{i2} = ['Sz' str int2str(i2) ' ' int2str(ts) ' to ' int2str(te) ' s'];
                Tnames{i2} = ['Sz' str int2str(i2) ' ' int2str(ts) ' to ' int2str(Tte) ' s'];
            else
                names{i2} = ['Sz' str int2str(i2)];
                Tnames{i2} = ['Sz' str int2str(i2)];
            end
            
            onsets{i2} = max(OldOnsets-tyv.TyvBefore,0)+(i2-1)*tyv.TyvWindow;
            durations{i2} = tyv.TyvWindow;
        end
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        
    end
end

%Save as .mat
fOnset = [file '.mat'];
too_few_spikes = 0;
if isempty(tyv)
    for i1 =1:length(onsets)
        if minimum_number_spikes >length(onsets{i1})
            too_few_spikes = 1;
        end
    end
end
if ~too_few_spikes
    if add_pulse_regressor
        try
            save(fOnset,'names','onsets','durations','TR','vRi')
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            save(fOnset,'names','onsets','durations','TR')
        end
    else
        save(fOnset,'names','onsets','durations','TR')
    end
end
if ~isempty(tyv)
    try
        save(fOnset,'names','Tnames','onsets','durations','TR')
    end
end
end

function run_stats_and_results(GLM)
tyv = GLM.tyv; %Tyvaert analysis
DirAnalysis = GLM.DirAnalysis;
DirOut = [GLM.Statslabel GLM.Vlabel GLM.Wlabel GLM.Clabel GLM.Glabel ...
    GLM.glabel GLM.Alabel GLM.Ulabel GLM.Dlabel GLM.Mlabel GLM.Plabel];
if ~isempty(tyv)
    DirOut = [DirOut '_' int2str(tyv.TyvPerGroup) 'perGroup'  int2str(tyv.TyvWindow) 'sPerWindow'];
end
typI = [GLM.Slabel GLM.Wlabel GLM.Atyp];
reportContrastsBySession = GLM.reportContrastsBySession;
flog = GLM.flog;
switch GLM.Dlabel
    case 'D'
        switch GLM.Glabel
            case 'G'
                nD = GLM.gamma_order;
            otherwise
                nD = 3;
        end
    otherwise
        nD = 1;
end
switch GLM.Vlabel %Volterra
    case 'V'
        nV = 2;
    otherwise
        nV = 1;
end
AnovaOn = GLM.AnovaOn;
extent = GLM.extent;
threshold = GLM.threshold;
mask_threshold = GLM.mask_threshold;
TR0 = GLM.TR0;
regenerate_reports = GLM.regenerate_reports;
DirResultsAll = GLM.DirResultsAll;
dStats = GLM.dStats;
%waittime = GLM.waittime;
f = GLM.f;
nE = length(f.fEPI); %number of sessions

str_noavg = 'No Avg'; %no average over sessions
str_replna = 'replna';
%str_noOnsAvg = 'No Ons Avg';
name_str = 'All Onsets';
str_pos = 'Positive';
str_neg = 'Negative';

SaveStatsBatch = GLM.SaveStatsBatch;
SaveReportBatch = GLM.SaveReportBatch;
SkipUncorrected = GLM.SkipUncorrected;
SkipGroupOnsets = GLM.SkipGroupOnsets;
SkipNegativeBOLD = GLM.SkipNegativeBOLD;
generate_patient_level = GLM.generate_patient_level;
%Notation introduced for results reports on F contrasts:
%a contrast is labeled in order by
%1- onset type = m = {A,N,m,0} (All (No average), iterate over several, doesn't apply because
%only one onset)
%2- Volterra = v = {A,v,0} (same: All (No average), iterate, not applicable, i.e. no Volterra)
%3- Derivatives = d = {A,1,0} (All (No average),first regressor, NA)
%4- response = r = {p,n,0} {positive, negative,0:both for Fcontrast)
%5- Session = s = {A (No average), G (All, average),s,0}

%results directory
run_report = 1;
if dStats
    sep = filesep;
    DirStats = [DirAnalysis f.DirSubj sep DirOut sep];
    if ~exist(DirStats,'dir'), mkdir(DirStats); end
    cd(DirStats);
    if exist('SPM.mat','file') && exist('beta_0001.img','file')
        write_log(flog,'Standard stats already run -- skipping (Delete SPM.mat in order to force rerun)');
        if ~regenerate_reports, run_report = 0; end
    else
        try
            try
                %load calculated TR
                TR = load(f.fOnset{1},'TR');
                TR = TR.TR;
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,strvcat('No TR value found in first onsets file, or no onsets file found',...
                    'Using TR value of 3.013 for Stats'));
                TR = TR0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %GLM Specification
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {DirStats};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
            for j=1:nE
                try
                    %Generate file names required
                    data = {};
                    for k=1:size(f.fEPI{j},1)
                        [dummy,fil ext] = fileparts(f.fEPI{j}(k).fname);
                        data = [data; [DirAnalysis f.DirSubj sep 'EPI' sep typI fil ext]];
                    end
                catch exception
                    disp(exception.identifier)
                    disp(exception.stack(1))
                    write_log(flog,'File specification error at Stats step');
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).scans = data;
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = {f.fOnset{j}};
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg = {[DirAnalysis f.DirSubj sep 'EPI' sep f.fMVT1{j}]};
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).hpf = 128;
            end
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            
            switch GLM.Dlabel
                case 'D'
                    switch GLM.Glabel
                        case 'G'
                            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = GLM.gamma_window;
                            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order = GLM.gamma_order;
                        otherwise
                            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
                    end
                otherwise
                    switch GLM.Glabel
                        case 'G'
                            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = GLM.gamma_window;
                            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order = 1;
                        otherwise
                            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    end
            end
            matlabbatch{1}.spm.stats.fmri_spec.volt = nV;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %GLM Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %GLM Contrast definition
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q = 1; %index to count contrast types
            if isempty(GLM.tyv)
                if generate_patient_level
                    SessRep = 'both'; %'sess';
                else
                    SessRep = 'sess';
                end
                
                matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                %Don't delete contrasts!
                matlabbatch{3}.spm.stats.con.delete = 0;
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
                matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                %Need to generalize this to a session-dependent number of
                %onset types -- not trivial!
                load(f.fOnset{1}); %just to have list of onset names
                M = size(names,2); %number of onset types
                if M>1 && ~SkipGroupOnsets
                    Nreg = nD*M+(nV-1)*nD^2*M*(M+1)/2;
                    cvtmpP = zeros(nV,Nreg); cvtmpN = cvtmpP; %for positive and negative T contrasts for All Onsets, for each Volterra
                    if nD > 1
                        cvtmpF = zeros(nD,Nreg); %F contrast averaged over onsets
                        if nV == 2
                            cvtmpF2 = zeros(nD*(nD+1)/2,Nreg);
                            cvtmpFV = zeros(nD+nD*(nD+1)/2,Nreg);
                        end
                    end
                    if nV == 2
                        cvtmpV = zeros(2,Nreg);
                    end
                    %Not averaged over onsets -
                    %c_mvdp
                    %c_Av0p = zeros(M,Nreg);
                    %c_Av0n = zeros(M,Nreg);
                    %c_A
                end
                %Loop over contrast types
                for m=1:M
                    %Loop over Volterra kernels
                    for v=1:nV
                        if nV == 2
                            Vstr = ['V' int2str(v)];
                        else
                            Vstr = '';
                        end
                        
                        %T contrast for activations: number = nE * M * Volt
                        reg = nD*(m-1)*(2-v)+(v-1)*(nD*M)+...
                            (v-1)*nD^2*( (m-1)*M - (m-1)*(m-2)/2); %requires some mental gymnastics!
                        Frow = nD*(2-v)+(v-1)*nD*(nD+1)/2;
                        if v == 1
                            Vmat = eye(nD);
                            reg1 = reg; %store position for use later when v == 2
                            Vmat1 = Vmat;
                        else %v == 2; for nD == 3
                            switch nD
                                case 1
                                    Vmat = 1;
                                case 3
                                    temp_mat = [zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
                                    Vmat = [eye(3) zeros(3,6); zeros(3,3) temp_mat];
                            end
                        end
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [names{m} Vstr str_pos];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,reg) 1];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                        %T contrast for deactivations: number = nE * M * Volt
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [names{m} Vstr str_neg];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,reg) -1];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                        if nD > 1
                            %F contrast for canonical + derivatives: one contrast for each Volterra
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [names{m} Vstr];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {[zeros(Frow,reg) Vmat]}';
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1; %average over sessions
                            %replicate
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [names{m} Vstr str_noavg];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1; %no average over sessions
                        end
                        if v == 2
                            Vstr2 = 'VA';
                            %combined Volterra  - position:
                            tmp_mat = [Vmat1; zeros(Frow,nD)];
                            tmp_mat2 = [zeros(nD,size(Vmat,2)); Vmat];
                            fVmat = [zeros(Frow+nD,reg1) tmp_mat zeros(Frow+nD,reg-reg1-nD) tmp_mat2]; %added nD as quick fix, needs to be thought through
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [names{m} Vstr2];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {fVmat}';
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1; %no average over V and average over sessions
                            %replicate
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [names{m} Vstr2 str_noavg];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1; %no average over V and no average over sessions
                        end
                        %Combine all spikes types into one
                        cvtmpP(v,reg+1) = 1; cvtmpN(v,reg+1) = -1; %Positive and negative T contrasts
                        if nD > 1
                            if v == 1
                                for j=1:nD, cvtmpF(j,reg+j) = 1; end %F contrast over derivatives
                                for j=1:nD, cvtmpFV(j,reg+j) = 1; end
                            else
                                for j=1:nD, cvtmpF2(j,reg+j) = 1; cvtmpFV(j+3,reg+j) = 1; end %F contrast over Volterra
                                cvtmpF2(4,reg+5) = 1; cvtmpF2(5,reg+6) = 1; cvtmpF2(6,reg+9) = 1;
                                cvtmpFV(7,reg+5) = 1; cvtmpFV(8,reg+6) = 1; cvtmpFV(9,reg+9) = 1;
                            end
                        end
                        if nV == 2
                            cvtmpV(v,reg+1) = 1; %F contrast ex derivatives for each Volterra
                        end
                    end
                end
                %Contrasts for combining all spike types into one
                if M > 1 && ~SkipGroupOnsets
                    for v=1:nV
                        if nV == 2
                            Vstr = ['V' int2str(v)];
                        else
                            Vstr = '';
                        end
                        %T contrast for activations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [name_str Vstr str_pos];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = cvtmpP(v,:); %Average over onsets
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                        %T contrast for deactivations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [name_str Vstr str_neg];
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = cvtmpN(v,:);
                        matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                        
                        if nD > 1
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {cvtmpF}'; %Avg onsets and sessions
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1;
                            %replicate
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr str_noavg];
                            matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1;
                            %str_noOnsAvg
                            
                            if v == 2
                                Vstr2 = 'VA';
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr2];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {cvtmpF2}';
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1;
                                %replicate
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr2 str_noavg];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1;
                                
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str ' Full F'];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {cvtmpFV}';
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1;
                                %replicate
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr str_noavg ' Full F'];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1;
                                
                            end
                        else
                            if v == 2
                                Vstr2 = 'VA';
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr2];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {cvtmpV}';
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = SessRep; q=q+1;
                                %replicate
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon = matlabbatch{3}.spm.stats.con.consess{q-1}.fcon;
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [name_str Vstr2 str_noavg];
                                matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = str_replna; q=q+1;
                            end
                        end
                    end
                end
                if AnovaOn
                    for m1=1:M-1
                        name1 = names{m1};
                        reg1 = nD*(m1-1);
                        for m2=m1+1:M
                            name2 = names{m2};
                            reg2 = nD*(m2-1);
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [name2 ' less ' name1 str_pos];
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,reg1) -1 zeros(1,reg2-reg1-1) 1];
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [name2 ' less ' name1 str_neg];
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,reg1) -1 zeros(1,reg2-reg1-1) 1];
                            matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = SessRep; q=q+1;
                        end
                    end
                end
            else
                %Tyvaert Analysis
                %NEED TO MAKE SURE THERE IS ONLY ONE SEIZURE PER SESSION AND ONLY ONE SESSION
                %If there are several sessions, they need to be treated as separate patients
                %hence there is only one onset type, for the one seizure
                matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                %Don't delete contrasts!
                matlabbatch{3}.spm.stats.con.delete = 0;
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
                matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                
                load(f.fOnset{1}); %for list of onset names
                M = size(names,2);
                %Generate T-stats for (de)activation
                for m=1:M
                    %T contrast for activations: (N+1) contrasts
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [names{m} 'Positive'];
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,m-1) 1];
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = 'none'; q=q+1;
                    %T contrast for deactivations: (N+1) contrasts
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.name = [names{m} 'Negative'];
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.convec = [zeros(1,m-1) -1];
                    matlabbatch{3}.spm.stats.con.consess{q}.tcon.sessrep = 'none'; q=q+1;
                end
                %Generate F-stats masked by (de)activation
                for m=1:(M-tyv.TyvPerGroup+1)
                    matlabbatch{3}.spm.stats.con.consess{q}.fcon.name = [names{m} 'F'];
                    matlabbatch{3}.spm.stats.con.consess{q}.fcon.convec = {[zeros(tyv.TyvPerGroup,m-1) eye(tyv.TyvPerGroup)]}';
                    matlabbatch{3}.spm.stats.con.consess{q}.fcon.sessrep = 'none'; q=q+1;
                end
            end
            if SaveStatsBatch
                save('StatsBatch','matlabbatch');
                x1 = matlabbatch{3}.spm.stats.con.consess;
                write_log(flog,'Specified Contrasts');
                for y1=1:length(x1)
                    try
                        write_log(flog,['T-Contrast ' int2str(y1) ': ' x1{y1}.tcon.name]);
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))
                        write_log(flog,['F-Contrast ' int2str(y1) ': ' x1{y1}.fcon.name]);
                    end
                end
            end
            spm_jobman('run',matlabbatch);
            
            %Redisplay for convenience
            if SaveStatsBatch
                x1 = matlabbatch{3}.spm.stats.con.consess;
                write_log(flog,'Specified Contrasts');
                for y1=1:length(x1)
                    try
                        write_log(flog,['T-Contrast ' int2str(y1) ': ' x1{y1}.tcon.name]);
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))
                        write_log(flog,['F-Contrast ' int2str(y1) ': ' x1{y1}.fcon.name]);
                    end
                end
                
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            if isempty(GLM.tyv)
                write_log(flog,'Standard Stats failed to run');
            else
                write_log(flog,'Tyvaert Stats failed to run');
            end
        end
    end %end if exist('SPM.mat','file')
end %end if dStats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLM Results Reports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dStats
    cd(DirStats);
    if run_report
        if exist('SPM.mat','file')
            spmmat = spm_select('List',pwd,'^SPM.mat');
            try
                clear matlabbatch
                %Results report
                load(spmmat); %load SPM
                x1 = SPM.xCon;
                write_log(flog,'Generated Contrasts');
                for y1=1:length(x1)
                    write_log(flog,[x1(y1).STAT '-Contrast ' int2str(y1) ': ' x1(y1).name]);
                end
                matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
                load(f.fOnset{1});
                M = size(names,2);
                %Results for group of sessions
                mm = 0;
                RO = []; %Report Options
                
                switch SkipUncorrected
                    case 1
                        RO.threshdesc = {'FWE'};
                        threshold = threshold(1);
                    case 2
                        RO.threshdesc = {'none'};
                        threshold = threshold(2);
                    case 0
                        RO.threshdesc = {'FWE' 'none'};
                end
                %Loop over onset types
                if M > 1 && ~SkipGroupOnsets
                    loop_start = 0;
                else
                    loop_start = 1;
                end
                if isempty(GLM.tyv)
                    if SkipNegativeBOLD
                        RO.mid_str = {'Positive'};
                        step_outreport = 1;
                    else
                        RO.mid_str = {'Positive' 'Negative'};
                        step_outreport = 2;
                    end
                    for m=loop_start:M
                        if m == 0
                            name = 'All Onsets';
                        else
                            name = names{m};
                        end
                        %Loop over Volterra kernels
                        for v=1:nV
                            if nV == 2
                                Vstr = ['V' int2str(v)];
                            else
                                Vstr = '';
                            end
                            RO.start_str = [name Vstr];
                            RO.sess_str = ' - All Sessions';
                            
                            RO.thresh = threshold;
                            RO.extent = extent;
                            RO.mask = [];
                            RO.F_str = [];
                            if generate_patient_level
                                if nE > 1
                                    %T-stats
                                    %loop over contrasts to find the desired one
                                    [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                    %Masked F-contrasts
                                    if nD > 1
                                        RO.F_str = str_noavg;
                                        RO.mask.thresh = mask_threshold;
                                        RO.mask.mtype = 0; %inclusive
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        RO.F_str = ''; %F of averaging type
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        if v==2
                                            Vstr2 = 'VA';
                                            Vstr = 'V1';
                                            RO.start_str = [name Vstr2];
                                            RO.mask.pretitlestr = [name Vstr];
                                            RO.F_str = str_noavg;
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                            RO.F_str = ''; %F of averaging type
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        end
                                    else
                                        if v==2
                                            Vstr2 = 'VA';
                                            Vstr = 'V1';
                                            RO.start_str = [name Vstr2];
                                            RO.mask.pretitlestr = [name Vstr];
                                            RO.mask.thresh = mask_threshold;
                                            RO.mask.mtype = 0;
                                            RO.F_str = str_noavg;
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                            RO.F_str = ''; %F of averaging type
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        end
                                    end
                                end
                            end
                            %by session
                            if (reportContrastsBySession && nE > 1) || nE == 1
                                for s1=1:nE
                                    RO.mask = [];
                                    RO.start_str = [name Vstr];
                                    RO.sess_str = [' - Session ' int2str(s1)]; %careful about spaces!
                                    [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                    %Masked F-contrasts
                                    if nD > 1
                                        RO.F_str = str_noavg;
                                        RO.mask.thresh = mask_threshold;
                                        RO.mask.mtype = 0; %inclusive
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        RO.F_str = ''; %F of averaging type
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        if v==2
                                            Vstr2 = 'VA';
                                            Vstr = 'V1';
                                            RO.start_str = [name Vstr2];
                                            RO.mask.pretitlestr = [name Vstr];
                                            RO.F_str = str_noavg;
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                            RO.F_str = ''; %F of averaging type
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        end
                                    else
                                        if v==2
                                            Vstr2 = 'VA';
                                            Vstr = 'V1';
                                            RO.mask.thresh = mask_threshold;
                                            RO.mask.mtype = 0;
                                            RO.start_str = [name Vstr2];
                                            RO.mask.pretitlestr = [name Vstr];
                                            RO.F_str = str_noavg;
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                            RO.F_str = ''; %F of averaging type
                                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if AnovaOn
                        RO.thresh = threshold(1);
                        RO.extent = extent;
                        RO.mask = [];
                        RO.F_str = [];
                        for m1=1:M-1
                            name1 = names{m1};
                            RO.sess_str = ' - All Sessions';
                            for m2=m1+1:M
                                name2 = names{m2};
                                RO.start_str = [name2 ' less ' name1];
                                if generate_patient_level
                                    if nE > 1
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                    end
                                end
                                if (reportContrastsBySession && nE > 1) || nE == 1
                                    for s1=1:nE
                                        RO.sess_str = [' - Session ' int2str(s1)]; %careful about spaces!
                                        [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                                    end
                                end
                            end
                        end
                    end
                else
                    %Tyvaert analysis
                    mid_str2 = {'Positive' 'Negative'};
                    RO.threshdesc = {'FWE'};% 'none'};
                    threshold = threshold(1);
                    for r1=1:length(mid_str2)
                        RO.mid_str = {mid_str2{r1}};
                        for m=1:(M-tyv.TyvPerGroup+1)
                            %F contrasts, masked
                            RO.start_str = [names{m} 'F'];
                            RO.sess_str = '';%' - Session 1';
                            RO.F_str = '';
                            if tyv.TyvAdjustMultipleFtests
                                RO.thresh = tyv.TyvStatThreshold/(M-tyv.TyvPerGroup+1); %adjustment for multiple F tests
                            else
                                RO.thresh = tyv.TyvStatThreshold; %no adjustment
                            end
                            %
                            RO.extent = extent;
                            RO.mask = [];
                            RO.F_str = [];
                            RO.mask.pretitlestr = names{m};
                            RO.mask.thresh = tyv.TyvStatThreshold; %no adjustment
                            RO.mask.mtype = 0; %inclusive
                            [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                            if tyv.TyvConjunctionMask
                                cm0 = matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts;
                                %Add other T contrasts to the mask as
                                %conjunction (intersection) of the contrasts
                                for cm1=1:(tyv.TyvPerGroup-1)
                                    %Careful: factor of 2 required when
                                    %generating both positive
                                    %and negative contrasts
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = ...
                                        [matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts cm0+cm1*2];
                                end
                            end
                        end
                        %                     for m=loop_start:M
                        %                         RO.mask = [];
                        %                         %try
                        %                         %    RO.start_str = Tnames{m};
                        %                         %catch
                        %                             RO.start_str = names{m};
                        %                         %end
                        %                         %T-contrasts
                        %                         [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO);
                        %                     end
                    end
                end
                matlabbatch{1}.spm.stats.results.units = 1;
                matlabbatch{1}.spm.stats.results.print = true;
                %             if SaveStatsBatch
                %                 x1 = SPM.xCon;
                %                 write_log(flog,'Generated Contrasts');
                %                 for y1=1:length(x1)
                %                     write_log(flog,[x1(y1).STAT '-Contrast ' int2str(y1) ': ' x1(y1).name]);
                %                 end
                %             end
                if SaveReportBatch
                    save('ReportBatch','matlabbatch');
                    try
                        x1 = matlabbatch{1}.spm.stats.results.conspec;
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))
                        disp('Found no contrast for reports -- there is a bug');
                    end
                    write_log(flog,'Generated Reports');
                    for y1=1:length(x1)
                        write_log(flog,['Report ' int2str(y1) ': ' x1(y1).titlestr]);
                    end
                end
                %end specification of Results Report
                spm_jobman('run',matlabbatch);
                %pause(waittime);
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,'Results Report failed to run');
            end
            
            %end of stats - copy Results Report
            try
                temp_file = spm_select('List',pwd,'^spm_');
                for nf=1:size(temp_file,1)
                    [dummy, fil ext] = fileparts(temp_file(nf,:));
                    if exist(temp_file(nf,:),'file')
                        copyfile(temp_file(nf,:),fullfile(DirResultsAll,[f.DirSubj '_' DirOut '_' fil ext]));
                    end
                end
            catch exception
                disp(exception.identifier)
                disp(exception.stack(1))
                write_log(flog,['Could not copy ' DirOut ' Results Report']);
            end
        end %end if exist('SPM.mat','file')
    end
end
end

%Normalize to Talairach-Tournoux atlas
function call_normalize(data,T1,DirSPM,flog,pfx_normalise)
try
    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {T1};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = data;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {[DirSPM 'templates\T1.nii,1']};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [1.5 1.5 1.5];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = pfx_normalise;
    spm_jobman('run',matlabbatch);
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    write_log(flog,'Normalize module failed to run');
end
end

function [matlabbatch mm] = search_con(SPM,matlabbatch,mm,RO)
%Loop over mid_str
for p1=1:length(RO.mid_str)
    if isempty(RO.mask)
        titlestr = [RO.start_str RO.mid_str{p1} RO.sess_str];
    else
        if isfield(RO.mask,'pretitlestr')
            titlestr = [RO.start_str RO.F_str RO.sess_str];
            mask.titlestr = [RO.mask.pretitlestr RO.mid_str{p1} RO.sess_str];
        else
            titlestr = [RO.start_str RO.F_str RO.sess_str];
            mask.titlestr = [RO.start_str RO.mid_str{p1} RO.sess_str];
        end
    end
    
    fc = 0; %found contrast
    fmc = 0; %found mask contrast
    %Loop over contrasts
    for c1=1:length(SPM.xCon)
        if strcmp(SPM.xCon(c1).name,titlestr)
            fc = c1;
            break
        end
    end
    if ~isempty(RO.mask)
        for c1=1:length(SPM.xCon)
            if strcmp(SPM.xCon(c1).name,mask.titlestr)
                fmc = c1;
                break
            end
        end
        if fmc==0
            disp(['Could not find mask contrast ' mask.titlestr])
            disp(['Therefore cannnot generate masked report batch for ' titlestr]);
            fc = 0;
        end
    end
    
    if fc
        try
            for d1=1:length(RO.threshdesc)
                mm = mm+1;
                if ~isempty(RO.mask)
                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [titlestr ' masked by ' RO.mid_str{p1}];
                else
                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = titlestr;
                end
                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = fc;
                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = RO.threshdesc{d1};
                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = RO.thresh(d1);
                matlabbatch{1}.spm.stats.results.conspec(mm).extent = RO.extent;
                if ~isempty(RO.mask)
                    matlabbatch{1}.spm.stats.results.conspec(mm).mask = RO.mask;
                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = fmc;
                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = RO.mask.thresh;
                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = RO.mask.mtype;
                else
                    matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                end
            end
        catch exception
            disp(exception.identifier)
            disp(exception.stack(1))
            mm = mm-1;
            disp(['Could not generate report batch for ' titlestr]);
        end
    end
end
end