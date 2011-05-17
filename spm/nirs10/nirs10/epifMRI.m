function epifMRI
%Batch for IRMf studies
%onsets to remove
%Code allows up to 3 onset types to be removed. If want to keep all onsets,
%make sure all rem_onsets are set to ''
rem_onsets{1} = 'spk biF'; %'spkLFT';
rem_onsets{2} = '';
rem_onsets{3} = '';

%normalize to TT atlas using Dartel
DartelOn = 0;
%Normalize to Talairach-Tournoux atlas using old SPM normalise procedure
normalizeOn = 1;
reportContrastsBySession = 1;
%Generate stats with HRF peaking at different times wrt onsets
McGilldelaysOn = 0;
%Delays in seconds
delay = [-2 -1 1]; % -2 -1 1 3 6]; % [-9 -6 -3 -2 -1 1 2 3 6 9]; %[-7 -8 -5 -4]; %[-9 -6 -3 3 6 9];
%remove negative onsets entirely
McGill_remove_negative = 1;
%Generate stats with square or first derivative of movement parameters
MovementOn = 0; %1; %1; %1;
%Adjust the onset to account that the middle slice starts TR/2 later
STC_middle_slice_adjust_onsets = 1;
%Do slice timing correction centered on the middle slice
STC_middle_slice = 1;
%When there are several types of stimuli, group them into the same type
force_group_onsets = 0;
%Regenerate Result Reports even if stats have already been calculated
regenerate_reports = 0;
%force running unwarp anyway even if movement is small
force_unwarp = 1;
%threshold to unwarp in mm
unwarp_threshold = 0.3;
%whether to include time derivative and dispersion of canonical HRF
%inc_derivs = 0;
add_pulse_regressor = 0;



try
    % Initialise SPM
    %--------------------------------------------------------------------------
    %spm('Defaults','fMRI');
    spm fmri;
    %spm_jobman('initcfg');
    %spm_check_installation;
    [t,sts] = spm_select([1,inf],'dir','Select folder(s) of data files to analyze (anatomical, functional and onsets)');
    if ~sts, return; end %fatal
catch
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

[dir2,~] = fileparts(t(1,:));
idx = findstr(sep,dir2);
temp_dir = dir2(1:idx(end));
try 
    [tA,sts] = spm_select(1,'dir','Select folder to put all the analysis',{temp_dir}); 
    if ~sts
        try
            tA = spm_input('Enter analysis directory',1,'s',temp_dir);
        catch
            disp('Problem with directory selection -- Aborting');
            return %fatal
        end 
    end
catch
    disp('Problem with directory selection -- Aborting');
    return %fatal
end
if ~exist(tA,'dir'), mkdir(tA); end
DirAnalysis = tA;

DirStart = pwd; %keep track of initial directory
cd(DirAnalysis);

%Start time counter
tStart = tic;

%Check if log file exists
temp_log = ['epifMRI_log' date];
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
catch
    disp('Could not open log file in current directory -- exiting');
    return; %fatal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check user chosen data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Result: create three structures containing info per subject and per session:
%f{i}.fT1
%f{i}.fEPI{j}
%f{i}.fOnset{j}
%dStats{i} Boolean if GLM stats should be calculated
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
                    else
                        %no good images
                        write_log(flog,'No good images. Code should not be coming here!');
                    end
                end %end if nImg == 1

            catch
                %not images - then check they are onsets
                %two cases: raw files from Analyzer2 or .mat structure of
                %names, onsets, durations
                try
                    %check if this is a .mat file  
                    load(filesRec{j,:},'names','onsets','durations');
                    %found a good onset file
                    nSess_onset = nSess_onset + 1;
                    f{aNsubj}.fOnset{nSess_onset} = filesRec{j,:};
                catch
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
                            [TR,temp_file] = readOnsets(filesRec{j,:},rem_onsets,add_pulse_regressor);
                            nSess_onset = nSess_onset + 1;
                            f{aNsubj}.fOnset{nSess_onset} = temp_file;
                            if TR < 0.5 || TR > 7
                                write_log(flog,['Warning: unusual TR value: ' num2str(TR)...
                                    ' for subject ' int2str(i)]);
                            end
                        end
                    catch %#ok<*CTCH>
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
    catch %try to find data for this subject
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
    if dStats{i}
        for j=1:size(f{i}.fEPI,2)
            write_log(flog,strvcat(['Session ' int2str(j) ': ' ...
                int2str(size(f{i}.fEPI{j},1)) ' volumes in ' ...
                f{i}.fEPI{j}(1).fname ' matched with ' f{i}.fOnset{j}]));
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
            catch
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
DirResultsAll = [DirAnalysis sep 'ResultsAll' sep];
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
    idx = findstr('epi',tempDirSubj);
    if isempty(idx)
        DirSubj = ['a_' tempDirSubj];
    else
        DirSubj = ['a_' tempDirSubj(idx:end)];
    end
    f{i}.DirSubj = DirSubj;
    if ~exist(DirSubj,'dir'), mkdir(DirSubj); end
    DirT1 = [DirSubj sep 'T1' sep]; 
    DirEPI = [DirSubj sep 'EPI' sep];
    DirEEG = [DirSubj sep 'EEG' sep];
    if ~exist(DirT1,'dir'), mkdir(DirT1); end
    if ~exist(DirEPI,'dir'), mkdir(DirEPI); end
    if ~exist(DirEEG,'dir'), 
        mkdir(DirEEG);
    else
        %Clean up! Otherwise, might compound delays for onsets      
        k = 1;
        while 1
            temp_dirEEG = [DirSubj sep 'EEG' '_' int2str(k)];
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
    %DirResults = [DirSubj sep 'Results' sep];
    %DirScripts = [DirSubj sep 'Scripts' sep];
    %DirDCM = [DirSubj sep 'DCM' sep];
    %create directories for each subject
    %if ~exist(DirResults,'dir'), mkdir(DirResults); end
    %if ~exist(DirScripts,'dir'), mkdir(DirScripts); end
    %if ~exist(DirDCM,'dir'), mkdir(DirDCM); end
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
    
    %Copy EPI and onset files if required -- looping over sessions
    for j=1:size(f{i}.fEPI,2) 
        %check file name of first EPI file for this session
        [~, filEPI ext2] = fileparts(f{i}.fEPI{j}(1).fname);        
        temp_file = fullfile(DirEPI,[filEPI ext2]);        
        if ~exist(temp_file,'file'), %use spm_existfile instead?
            if strcmp(ext2,'.img')
                 %loop over time
                for k=1:size(f{i}.fEPI{j},1)
                    [~, filEPI ext2] = fileparts(f{i}.fEPI{j}(k).fname);  
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
            [~, filOnset ext2] = fileparts(f{i}.fOnset{j}); 
            temp_file = fullfile(DirEEG,[filOnset ext2]);
            if ~exist(temp_file,'file'), 
                copyfile(f{i}.fOnset{j},temp_file); 
            end  
            %Update the location of the file
            f{i}.fOnset{j} = fullfile(DirAnalysis, temp_file);
        end                
    end       
end
catch
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
[~, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
cd([DirAnalysis f{i}.DirSubj sep 'EPI']);

if spm_existfile(['r' fil ext(1:end-2)])
    %Only check for first session
    write_log(flog,'Realign module already run -- skipping');
else
    try
        %Generate file names required 
        data = {};        
        for j=1:size(f{i}.fEPI,2)
            %[~, fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
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
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm_jobman('run',matlabbatch);
    catch
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
    dp1 = [zeros(1,6); diff(p1,1,1)];
    %Combine derivatives of motion parameters with motion parameters themselves 
    dp1 = [p1 dp1];

    f{i}.fMVT2{j} = ['rpd' f{i}.fMVT1{j}(3:end)];
    save(f{i}.fMVT2{j},'dp1','-ascii');
    %squares of parameters
    f{i}.fMVT22{j} = ['rp2' f{i}.fMVT1{j}(3:end)];
    p2 = p1 .* p1;
    p2 = [p1 p2];
    save(f{i}.fMVT22{j},'p2','-ascii');          
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
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
            spm_jobman('run',matlabbatch);            
        else
            write_log(flog,'Little movement on all sessions; no need to run unwarp');
        end
    catch
        write_log(flog,'Unwarp module failed');
    end
end

%New Segment
cd([DirAnalysis f{i}.DirSubj sep 'T1']);
[~, fil ext] = fileparts(f{i}.fT1.fname);
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
    catch
        write_log(flog,'New Segment failed to run');
    end
end

%Return to EPI directory
cd([DirAnalysis f{i}.DirSubj sep 'EPI']);

try 
    load('coreg_OK');
    %Coregistration already done
    write_log(flog,'Coregistration already done -- skipping');
catch
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
            [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
            if spm_existfile(['u' fil ext(1:end-2)])
                %run coregistration on 'u' files 
                typ = 'u';
            else
                %run coregistration on 'r' files
                typ = 'r';
            end           
            for k=1:size(f{i}.fEPI{j},1)   
                [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
            end
        end
    catch
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
        cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
        save('coreg_OK','coreg_OK');
    catch
        write_log(flog,'Coregistration module failed to run');
    end
end

%Slice timing correction
[~, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
if spm_existfile(['ar' fil ext(1:end-2)]) || spm_existfile(['au' fil ext(1:end-2)])
    write_log(flog,'Slice Timing Correction already run -- skipping');
else
    try
        %Generate file names required 
        data = {};        
        for j=1:size(f{i}.fEPI,2)
            [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
            if spm_existfile(['u' fil ext(1:end-2)])
                %run slice timing correction on 'u' files 
                typ = 'u';
            else
                %run slice timing correction on 'r' files
                typ = 'r';
            end
            data_sess = {};
            for k=1:size(f{i}.fEPI{j},1)   
                [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                data_sess= [data_sess; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
            end
            data= [data; {data_sess}];
        end
    catch
        write_log(flog,'File specification error at Slice Timing Correction step');
    end
    try
        %load calculated TR
        TR = load(f{i}.fOnset{1},'TR');
        TR = TR.TR;
    catch
        write_log(flog,strvcat('No TR value found in first onsets file, or no onsets file found',...
            'Using TR value of 3.013 for slice timing correction'));
        TR = 3.013;
    end
    nslices = 47;
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
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        spm_jobman('run',matlabbatch);
    catch
        write_log(flog,'Slice Timing Correction module failed to run');
    end
end


if normalizeOn
    if DartelOn
        %Return to T1 directory
        cd([DirAnalysis f{i}.DirSubj sep 'T1']);
        [~, fil ext] = fileparts(f{i}.fT1.fname);
        %DARTEL
        if spm_existfile('Template_6.nii')
            write_log(flog,'DARTEL already run -- skipping');
        else        
            try
                %Generate file names required 
                data = {};        
                for j=1:size(f{i}.fEPI,2)
                    [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                    if spm_existfile(['au' fil ext(1:end-2)])
                        %run smooth on 'au' files 
                        typ = 'au';
                    else
                        %run smooth on 'ar' files
                        typ = 'ar';
                    end           
                    for k=1:size(f{i}.fEPI{j},1)   
                        [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                        data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
                    end
                end
                %Normalize some of the anatomical images
                data = {[data; f{i}.fT1.fname ',1']};
            catch
                write_log(flog,'File specification error for Dartel');
            end

            try
                [~, fil ext] = fileparts(f{i}.fT1.fname);
                data2 = {{[DirAnalysis f{i}.DirSubj sep 'T1' sep 'rc1' fil ext]}
                        {[DirAnalysis f{i}.DirSubj sep 'T1' sep 'rc2' fil ext]}};
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
                matlabbatch{2}.spm.tools.dartel.mni_norm.fwhm = [8 8 8];  
                spm_jobman('run',matlabbatch);
            catch
                %Run Normalise instead

                try
                    %Generate file names required 
                    data = {};        
                    for j=1:size(f{i}.fEPI,2)
                        [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                        if spm_existfile(['au' fil ext(1:end-2)])
                            %run smooth on 'au' files 
                            typ = 'au';
                        else
                            %run smooth on 'ar' files
                            typ = 'ar';
                        end           
                        for k=1:size(f{i}.fEPI{j},1)   
                            [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                            data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
                        end
                    end
                    %Normalize some of the anatomical images
                    data = {[data; f{i}.fT1.fname]};

                catch
                    write_log(flog,'File specification error for Normalise');
                end

                try
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {f{i}.fT1.fname};
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
                    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
                    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
                    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
                    spm_jobman('run',matlabbatch);
                catch
                    write_log(flog,'Normalize module failed to run on EPI images after DARTEL failed');
                end
                
                %Normalize anatomical T1 image
                data = {f{i}.fT1.fname};
                try
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {f{i}.fT1.fname};
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
                    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
                    spm_jobman('run',matlabbatch);
                catch
                    write_log(flog,'Normalize module failed to run for T1 file after DARTEL failed');
                end
        
            end
        end
    else %if DartelOn
        
        %Run Normalise instead

        
        [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
        if spm_existfile(['wau' fil ext(1:end-2)]) || spm_existfile(['war' fil ext(1:end-2)])
            write_log(flog,'Normalized on war or wau files already run -- skipping');
        else
            try
                %Generate file names required 
                data = {};        
                for j=1:size(f{i}.fEPI,2)
                    [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                    if spm_existfile(['au' fil ext(1:end-2)])
                        %run smooth on 'au' files 
                        typ = 'au';
                    else
                        %run smooth on 'ar' files
                        typ = 'ar';
                    end           
                    for k=1:size(f{i}.fEPI{j},1)   
                        [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                        data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
                    end
                end
            

            catch
                write_log(flog,'File specification error for Normalise');
            end

            try
                clear matlabbatch
                matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {f{i}.fT1.fname};
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
                matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [3 3 3];
                matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
                matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,'Normalize module failed to run for EPI files');
            end
        
            %Normalize anatomical T1 image
            data = {f{i}.fT1.fname};
            try
                clear matlabbatch
                matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {f{i}.fT1.fname};
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
                matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,'Normalize module failed to run for T1 file');
            end
        end
    end %end if DartelOn

end %end if normalizeOn


%Smooth: several runs
%1 - Smooth ar or au files (or directly r or u files)
%2 - Smooth war or wau files
%Return to EPI directory
cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
[~, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
if spm_existfile(['sar' fil ext(1:end-2)]) || spm_existfile(['sau' fil ext(1:end-2)])
    write_log(flog,'Smooth on sar files already run -- skipping');
else
    try
        %Generate file names required 
        data = {};        
        for j=1:size(f{i}.fEPI,2)
            [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
            if spm_existfile(['au' fil ext(1:end-2)])
                %run smooth on 'au' files 
                typ = 'au';
            else
                %run smooth on 'ar' files
                typ = 'ar';
            end           
            for k=1:size(f{i}.fEPI{j},1)   
                [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
            end
        end
    catch
        write_log(flog,'File specification error at Smooth step');
    end
    
    try
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = data;
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);
    catch
        write_log(flog,'Smooth module on ar/au files failed to run');
    end
end

%Smooth #2 (need to call a function rather than repeat code here
if normalizeOn
    cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
    [~, fil ext] = fileparts(f{i}.fEPI{1}(1).fname);
    if spm_existfile(['swar' fil ext(1:end-2)]) || spm_existfile(['swau' fil ext(1:end-2)])
        write_log(flog,'Smooth on swar files already run -- skipping');
    else
        try
            %Generate file names required 
            data = {};        
            for j=1:size(f{i}.fEPI,2)
                [~,fil ext] = fileparts(f{i}.fEPI{j}(1).fname);
                if spm_existfile(['wau' fil ext(1:end-2)])
                    %run smooth on 'au' files 
                    typ = 'wau';
                else
                    %run smooth on 'ar' files
                    typ = 'war';
                end           
                for k=1:size(f{i}.fEPI{j},1)   
                    [~,fil ext] = fileparts(f{i}.fEPI{j}(k).fname);
                    data = [data; [DirAnalysis f{i}.DirSubj sep 'EPI' sep typ fil ext]];               
                end
            end
        catch
            write_log(flog,'File specification error at Smooth step');
        end

        try
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = data;
            matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            spm_jobman('run',matlabbatch);
        catch
            write_log(flog,'Smooth module on war/wau files failed to run');
        end
    end
end

%end of preprocessing: copy ghostview files
try 
    cd([DirAnalysis f{i}.DirSubj sep 'EPI']);
    temp_file = spm_select('List',pwd,'^spm_');   
    for nf=1:size(temp_file,1)
        [~, fil ext] = fileparts(temp_file(nf,:));
        if exist(temp_file(nf,:),'file')
            copyfile(temp_file(nf,:),fullfile(DirResultsAll,[f{i}.DirSubj '_PreProcess_' fil ext]));
        end
    end
catch
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
            try 
                TR;
            catch
                TR = 3;
                write_log(flog,'TR was not found in onset file at onset adjust step');
            end
            for m=1:size(names,2)
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
catch
    write_log(flog,'Could not adjust onsets for STC_middle_slide');
end

%Check if there are the same number and types of onsets in each onset file
try
    permute_onsets = 0;
    %cd([DirAnalysis DirSubj]);
    DirOtherOnsets = fullfile(DirAnalysis, f{i}.DirSubj, 'EEG','OtherOnsets');
    if ~exist(DirOtherOnsets,'dir'), mkdir(DirOtherOnsets); end
    %Onset file, Movement
    load(f{i}.fOnset{1});
    names2 = names;
    f{i}.GroupOnsets = 0;
    for j=2:size(f{i}.fEPI,2)
        load(f{i}.fOnset{j});
        if size(names2,2) == size(names,2)
            for m=1:size(names,2)
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
            write_log(flog,strvcat(['Subject ' int2str(i) ': Different number of onsets in session ' int2str(j)],...
                'Onsets will be grouped into the same type'));
            f{i}.GroupOnsets = 1;
            %Don't bother trying to permute onsets if number of onsets
            %doesn't match
            permute_onsets = 0;
        end
    end
    
    if permute_onsets && ~GroupOnsets
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
            [~,fil ext] = fileparts(f{i}.fOnset{j});
            save(fullfile(DirOtherOnsets,[fil '_grouped' ext]),'names','onsets','durations','TR');
        end      
    end

    %Move onset files that will not be used
    if f{i}.GroupOnsets
        DirNotUsed = fullfile(DirAnalysis,f{i}.DirSubj,'EEG','NotUsed');
        if ~exist(DirNotUsed,'dir'), mkdir(DirNotUsed); end
        for j=1:size(f{i}.fEPI,2)
            [~,fil ext] = fileparts(f{i}.fOnset{j});
            movefile(f{i}.fOnset{j},[DirNotUsed '\' fil ext]);
            f{i}.fOnset{j} = fullfile(DirOtherOnsets,[fil '_grouped' ext]);
        end
    end
catch
    write_log(flog,'Could not group onset files');
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
        end
    catch
    end
else
    Clabel = '';
end
    
%Run with Derivs 
noDerivs = 0; 
%Standard statistics
run_stats_and_results(f{i},dStats{i},['Stats' Clabel],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);

%Statistics on normalized w files
%if normalizeOn, run_stats_and_results(f{i},dStats{i},['StatsW' Clabel],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs); end

%Volterra
run_stats_and_results_Volterra(f{i},dStats{i},['StatsVolt' Clabel],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs)

%Volterra on normalized files
%run_stats_and_results_Volterra(f{i},dStats{i},['StatsVoltW' Clabel],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs)

%
try
    if McGilldelaysOn
        %delays 
        for d=1:size(delay,2)
            g = f{i};
            for j=1:size(f{i}.fEPI,2)
                load(f{i}.fOnset{j});
                for m=1:size(names,2)
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
                [~,fil ext] = fileparts(f{i}.fOnset{j});                
                g.fOnset{j} = fullfile(DirOtherOnsets,[ fil '_' int2str(delay(d)) 's' ext]);
                save(g.fOnset{j},'names','onsets','durations','TR');
            end

            %Statistics with onset delays, to move peak of HRF
            %run_stats_and_results(g,dStats{i},['StatsW_' int2str(delay(d)) 's delay'],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog);
            run_stats_and_results(g,dStats{i},['Stats'  Clabel '_' int2str(delay(d)) 's delay'],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);
            %run_stats_and_results(g,dStats{i},['Stats_' int2str(delay(d))
            %'s delay'],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog);
        end %end for d=1:size(delays,2)
    end
catch
    write_log(flog,'McGill delays stats failed to run');
end

%Movement
try
    if MovementOn
        g = f{i};
        for j=1:size(f{i}.fMVT1,2)
            g.fMVT1{j} = f{i}.fMVT2{j};
        end
        %run_stats_and_results(g,dStats{i},'StatsW_DerMVT','swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog);
        run_stats_and_results(g,dStats{i},['StatsW' Clabel '_DerMVT'],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);

        g = f{i};
        for j=1:size(f{i}.fMVT1,2)
            g.fMVT1{j} = f{i}.fMVT22{j};
        end
        %run_stats_and_results(g,dStats{i},'StatsW_sqrMVT','swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog);
        run_stats_and_results(g,dStats{i},['StatsW' Clabel '_sqrMVT'],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);
    end
catch
    write_log(flog,'Movement stats failed to run');
end

% %Without derivatives
noDerivs = 1;
try
    %Statistics on normalized w files
    %if normalizeOn, run_stats_and_results(f{i},dStats{i},['StatsWnoD' Clabel],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs); end
    run_stats_and_results(f{i},dStats{i},['StatsnoD' Clabel],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);
    %run_stats_and_results_Volterra(f{i},dStats{i},['StatsVoltWnoD' Clabel],'swa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);
    %run_stats_and_results_Volterra(f{i},dStats{i},['StatsVoltnoD' Clabel],'sa',DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs);
catch
    write_log(flog,'Stats without derivs failed to run');
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

function [TR,fOnset] = readOnsets(file,rem_onsets,add_pulse_regressor)
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
    if ~strcmpi(tag,'Scan Start') && ~strcmpi(tag,'mvt') && ~strncmp(tag,'R',1) ...
             && ~strcmpi(tag,'Bad Interval') && ~strncmp(tag,'B',1) ...
             && (~strncmp(tag,'S',1) || strncmpi(tag,'Spk',3) || strncmpi(tag,'SWD',3)) && ~strncmp(tag,'T',1) ...
             && (~strcmp(tag,rem_onsets{1})) ...
             && (~strcmp(tag,rem_onsets{2})) ...
             && (~strcmp(tag,rem_onsets{3}))
        %Possible new onset
        not_found = 1; %true
        for i=1:NTonset
            %check whether onset is one of those already found
            if strcmpi(names{i},tag) %case insensitive
                not_found = 0; %found
                k{i} = k{i}+1; %#ok<*AGROW>                
                onsets{i}(k{i}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                durations{i}(k{i}) = dt*(str2double(temp(indx(3):indx(4)))-1); %subtract 1 so point onsets have zero duration
                break
            end
        end
        if not_found
            %new type of onset
            NTonset = NTonset + 1;
            k{NTonset} = 1;  
            names{NTonset} = tag;
            onsets{NTonset}(k{NTonset}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
            durations{NTonset}(k{NTonset}) = dt*(str2double(temp(indx(3):indx(4)))-1); %subtract 1 so point onsets have zero duration                   
        end
    else
        %ignore pulses, etc., before first Scan Start
        %if TRdone
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
        %end
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
    [~,idx] = find(dvR > medPR*(1+var_allowed)); %??
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
catch 
end
    
%Save as .mat
fOnset = [file '.mat'];
if add_pulse_regressor
    try
        save(fOnset,'names','onsets','durations','TR','vRi')
    catch
        save(fOnset,'names','onsets','durations','TR')
    end
else
    save(fOnset,'names','onsets','durations','TR')
end
end

function run_stats_and_results(f,dStats,DirOut,typI,DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs)
%Inputs f, dStats, output directory, type of files, analysis directory,
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
                catch
                    write_log(flog,strvcat('No TR value found in first onsets file, or no onsets file found',...
                        'Using TR value of 3.013 for Stats'));
                    TR = 3.013;
                end


                clear matlabbatch
                matlabbatch{1}.spm.stats.fmri_spec.dir = {DirStats};
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

                for j=1:size(f.fEPI,2)
                    try
                        %Generate file names required 
                        data = {};        
                        %for j=1:size(f.fEPI,2)
                        [~,fil ext] = fileparts(f.fEPI{j}(1).fname);
                        if spm_existfile([DirAnalysis f.DirSubj sep 'EPI' sep typI 'u' fil ext(1:end-2)])
                            %run on 'au' files 
                            typ = [typI 'u'];
                        else
                            %run on 'ar' files
                            typ = [typI 'r'];
                        end           
                        for k=1:size(f.fEPI{j},1)   
                            [~,fil ext] = fileparts(f.fEPI{j}(k).fname);
                            data = [data; [DirAnalysis f.DirSubj sep 'EPI' sep typ fil ext]];               
                        end
                    catch
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
                if noDerivs
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                else
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
                end
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
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
                matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
                matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                load(f.fOnset{1});
                if size(f.fEPI,2) > 1
                    %if there is more than one session  
                    %Get 3 N + 4 contrasts for each onset type
                    for m=1:size(names,2)
                        if noDerivs
                             %T contrast for activations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{1+2*(m-1)}.tcon.name = [names{m} 'Activation'];
                            matlabbatch{3}.spm.stats.con.consess{1+2*(m-1)}.tcon.convec = [zeros(1,(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{1+2*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{2+2*(m-1)}.tcon.name = [names{m} 'Deactivation'];
                            matlabbatch{3}.spm.stats.con.consess{2+2*(m-1)}.tcon.convec = [zeros(1,(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{2+2*(m-1)}.tcon.sessrep = 'both';
                        else
                            %F contrast for canonical + derivatives: one contrast 
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.name = [names{m} 'No Avg']; 
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.convec = {[zeros(3,3*(m-1)) eye(3)]}';
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.sessrep = 'replna';
                            %F contrast for canonical + derivatives: (N+1) contrasts for N sessions
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.name = names{m};
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.convec = {[zeros(3,3*(m-1)) eye(3)]}';
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.sessrep = 'both';
                            %T contrast for activations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.name = [names{m} 'Activation'];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.name = [names{m} 'Deactivation'];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.sessrep = 'both';
                        end
                    end
                    %Gives one contrast for movement, the last contrast
                    if noDerivs
                        matlabbatch{3}.spm.stats.con.consess{1+2*size(names,2)}.fcon.name = 'MVT';
                        matlabbatch{3}.spm.stats.con.consess{1+2*size(names,2)}.fcon.convec = {[zeros(6,size(names,2)) eye(6)]}';
                        matlabbatch{3}.spm.stats.con.consess{1+2*size(names,2)}.fcon.sessrep = 'both';
                    else
                        matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.name = 'MVT';
                        matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.convec = {[zeros(6,3*size(names,2)) eye(6)]}';
                        matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.sessrep = 'both';
                        
                    end
                    
                    %Don't delete contrasts!
                    matlabbatch{3}.spm.stats.con.delete = 0;

                    %Contrasts for combining all spike types into one 
                    if size(names,2) > 1
                        convec_temp = []; convec_tempA = []; convec_tempD = [];
                        for m=1:size(names,2)
                            if noDerivs
                                convec_temp = [convec_temp eye(1)];
                                convec_tempA = [convec_tempA 1 ];
                                convec_tempD = [convec_tempD -1 ];
                            else
                                convec_temp = [convec_temp eye(3)];
                                convec_tempA = [convec_tempA 1 0 0];
                                convec_tempD = [convec_tempD -1 0 0];
                            end
                        end
                        if noDerivs
                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{2+2*size(names,2)}.tcon.name = 'All spikes combined Activation';
                            matlabbatch{3}.spm.stats.con.consess{2+2*size(names,2)}.tcon.convec = convec_tempA;
                            matlabbatch{3}.spm.stats.con.consess{2+2*size(names,2)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{3+2*size(names,2)}.tcon.name = 'All spikes combined Deactivation';
                            matlabbatch{3}.spm.stats.con.consess{3+2*size(names,2)}.tcon.convec = convec_tempD;
                            matlabbatch{3}.spm.stats.con.consess{3+2*size(names,2)}.tcon.sessrep = 'both'; 
                            
                        else
                            %F contrast on all spike types [eye(3) eye(3) ...]
                            matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.name = 'All spikes combined No Avg';                 
                            matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.convec = {convec_temp}';
                            matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.sessrep = 'replna';

                            matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.name = 'All spikes combined';                 
                            matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.convec = {convec_temp}';
                            matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.sessrep = 'both';

                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.name = 'All spikes combined Activation';
                            matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.convec = convec_tempA;
                            matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.name = 'All spikes combined Deactivation';
                            matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.convec = convec_tempD;
                            matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.sessrep = 'both'; 
                        end
                    end %end if size(names,2) > 1

                else
                    %if there is only one session, the second group of
                    %contrasts disappears
                    %Get 3 contrasts for each onset type
                    for m=1:size(names,2)
                        if noDerivs
                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{-1+2*m}.tcon.name = [names{m} 'Activation'];
                            matlabbatch{3}.spm.stats.con.consess{-1+2*m}.tcon.convec = [zeros(1,(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{-1+2*m}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{2*m}.tcon.name = [names{m} 'Deactivation'];
                            matlabbatch{3}.spm.stats.con.consess{2*m}.tcon.convec = [zeros(1,(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{2*m}.tcon.sessrep = 'both';
                        else                           
                            %F contrast for canonical + derivatives: Get one
                            %contrast
                            matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.name = [names{m}];
                            matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.convec = {[zeros(3,3*(m-1)) eye(3)]}';
                            matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.sessrep = 'replna';
                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.name = [names{m} 'Activation'];
                            matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.name = [names{m} 'Deactivation'];
                            matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.sessrep = 'both';
                        end
                    end
                    if noDerivs
                        matlabbatch{3}.spm.stats.con.consess{1+2*(size(names,2)-1)}.fcon.name = 'MVT';
                        matlabbatch{3}.spm.stats.con.consess{1+2*(size(names,2)-1)}.fcon.convec = {[zeros(6,3*size(names,2)) eye(6)]}';
                        matlabbatch{3}.spm.stats.con.consess{1+2*(size(names,2)-1)}.fcon.sessrep = 'both';

                    else
                        matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.name = 'MVT';
                        matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.convec = {[zeros(6,3*size(names,2)) eye(6)]}';
                        matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.sessrep = 'both';
                    end
                    %Contrasts for combining all spike types into one 
                    if size(names,2) > 1
                        convec_temp = []; convec_tempA = []; convec_tempD = [];
                        for m=1:size(names,2)
                            if noDerivs
                                convec_temp = [convec_temp eye(1)];
                                convec_tempA = [convec_tempA 1 ];
                                convec_tempD = [convec_tempD -1];
                            else
                                convec_temp = [convec_temp eye(3)];
                                convec_tempA = [convec_tempA 1 0 0];
                                convec_tempD = [convec_tempD -1 0 0];
                            end
                        end
                        if noDerivs
                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{2+2*(size(names,2)-1)}.tcon.name = 'All spikes combined Activation';
                            matlabbatch{3}.spm.stats.con.consess{2+2*(size(names,2)-1)}.tcon.convec = convec_tempA;
                            matlabbatch{3}.spm.stats.con.consess{2+2*(size(names,2)-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{3+2*(size(names,2)-1)}.tcon.name = 'All spikes combined Deactivation';
                            matlabbatch{3}.spm.stats.con.consess{3+2*(size(names,2)-1)}.tcon.convec = convec_tempD;
                            matlabbatch{3}.spm.stats.con.consess{3+2*(size(names,2)-1)}.tcon.sessrep = 'both';

                        else
                            %F contrast on all spike types [eye(3) eye(3) ...]
                            matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.name = 'All spikes combined No Avg';                 
                            matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.convec = {convec_temp}';
                            matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.sessrep = 'replna';
                            %T contrast for activations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.name = 'All spikes combined Activation';
                            matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.convec = convec_tempA;
                            matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: 1 contrast
                            matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.name = 'All spikes combined Deactivation';
                            matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.convec = convec_tempD;
                            matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.sessrep = 'both';
                        end
                    end %end if size(names,2) > 1
                end %if size(f.fEPI,2) > 1

                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,'Standard Stats failed to run');
            end 
        end %end if exist('SPM.mat','file')
    end %end if dStats  
    if dStats
    cd(DirStats);
    if run_report
        if exist('SPM.mat','file')
            spmmat = spm_select('List',pwd,'^SPM.mat'); 
            try
                clear matlabbatch
                    %Results report
                    matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
                    load(f.fOnset{1});
                    if noDerivs
                        %No need to separate case with one or more than one
                        %session
                        mm = 0;
                        for m=1:size(names,2)
                            
                            cna = 1 + 2* (size(f.fEPI,2)+1)*(m-1);
                            cnd = 2 + size(f.fEPI,2)+ 1 + 2* (size(f.fEPI,2)+1)*(m-1);
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. T-Activation'];;
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                       
                        
                        end
                        
                    else    
                        mm = 0;
                        for m=1:size(names,2)
                            mm = mm + 1;
                            if size(f.fEPI,2) > 1
                                %contrast number for F-stat of type No Averaging
                                cn = 1 + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation of type Averaging
                                cna = 3 + 2* size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation of type Averaging
                                cnd = 4 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast for movement                  
                            else
                                cn  = 1 + 3 * (m-1);
                                cna = 2 + 3 * (m-1);
                                cnd = 3 + 3 * (m-1);
                            end
                            %Masked activations
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            %Masked deactivations
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                            mm = mm + 1;
                            matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation'];
                            matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                            matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                            matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                            matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                            matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        end %end for m=1:size(names,2)
                        if reportContrastsBySession && size(f.fEPI,2) > 1
                            %Masked activations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked activations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked activations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked activations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j

                            %Masked deactivations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked deactivations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked deactivations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j
                            %Masked deactivations
                            for j=1:size(f.fEPI,2)
                                for m=1:size(names,2)
                                    %contrast number for F-stat for session j and contrast m
                                    cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T activation 
                                    cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                    %contrast number for T deactivation 
                                    cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                    mm = mm + 1;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation Session ' int2str(j)];
                                    matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                    matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                    matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                                end %end for m
                            end %end for j

                        end %end if reportContrastsBySession && size(f.fEPI,2) > 1
                        %Movement
                        if size(f.fEPI,2) > 1                   
                            %contrast for movement
                            cnm = 6 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (size(names,2)-1);
                        else

                            cnm = 4 + 3 * (size(names,2)-1); 
                        end
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = 'FWE Movement';
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnm;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = 'unc. Movement';
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnm;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                        %Add various unmasked stats
                    end
                matlabbatch{1}.spm.stats.results.units = 1;
                matlabbatch{1}.spm.stats.results.print = true;
                %end specification of Results Report

                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,'Results Report failed to run');
            end   

            %end of stats - copy Results Report
            try
                temp_file = spm_select('List',pwd,'^spm_');      
                for nf=1:size(temp_file,1)
                    [~, fil ext] = fileparts(temp_file(nf,:));
                    if exist(temp_file(nf,:),'file')
                        copyfile(temp_file(nf,:),fullfile(DirResultsAll,[f.DirSubj '_' DirOut '_' fil ext]));
                    end
                 end
            catch 
                    write_log(flog,['Could not copy ' DirOut ' Results Report']);
            end
        end %end if exist('SPM.mat','file')
    end
    end
end

function run_stats_and_results_Volterra(f,dStats,DirOut,typI,DirAnalysis,DirResultsAll,regenerate_reports,reportContrastsBySession,flog,noDerivs)
%Inputs f, dStats, output directory, type of files, analysis directory,
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
                catch
                    write_log(flog,strvcat('No TR value found in first onsets file, or no onsets file found',...
                        'Using TR value of 3.013 for Stats'));
                    TR = 3.013;
                end
                temp_mat = [ zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
                mat_Volt = [zeros(3,9); eye(3) zeros(3,6); zeros(3,3) temp_mat]; %a 9x9 matrix
                    
                clear matlabbatch
                matlabbatch{1}.spm.stats.fmri_spec.dir = {DirStats};
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

                for j=1:size(f.fEPI,2)
                    try
                        %Generate file names required 
                        data = {};        
                        %for j=1:size(f.fEPI,2)
                        [~,fil ext] = fileparts(f.fEPI{j}(1).fname);
                        if spm_existfile([DirAnalysis f.DirSubj sep 'EPI' sep typI 'u' fil ext(1:end-2)])
                            %run on 'au' files 
                            typ = [typI 'u'];
                        else
                            %run on 'ar' files
                            typ = [typI 'r'];
                        end           
                        for k=1:size(f.fEPI{j},1)   
                            [~,fil ext] = fileparts(f.fEPI{j}(k).fname);
                            data = [data; [DirAnalysis f.DirSubj sep 'EPI' sep typ fil ext]];               
                        end
                    catch
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
                if noDerivs
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                else
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
                end
                matlabbatch{1}.spm.stats.fmri_spec.volt = 2;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
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
                matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
                matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                load(f.fOnset{1});
                if size(f.fEPI,2) > 1
                    %if there is more than one session  
                    %Get 3 N + 4 contrasts for each onset type
                    for m=1:size(names,2)
                        if noDerivs
                             %T contrast for activations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.tcon.name = [names{m} 'p1'];
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.tcon.convec = [zeros(1,(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.tcon.name = [names{m} 'n1'];
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.tcon.convec = [zeros(1,(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.tcon.sessrep = 'both';
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.name = [names{m} 'p2'];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.convec = [zeros(1,(m-1)) 0 1];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.name = [names{m} 'n2'];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.convec = [zeros(1,(m-1)) 0 -1];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.sessrep = 'both';
                            
                        else
                            %F contrast for canonical + derivaties: one contrast 
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.name = [names{m} 'No Avg']; 
                            temp_mat1 = [zeros(3,3*(m-1)) eye(3) zeros(3,3*(size(names,2)-m)); ...
                                zeros(6,3*size(names,2))]; %this is 9 x (3*size(names,2))
                            %this will not work if size(names,2) > 2: need to
                            %locate in the design matrix the interactions of stimulus 2 with
                            %stimulus 3
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.convec = ...
                                {[temp_mat1 zeros(9,18*(m-1)) mat_Volt]}';
                            %{[zeros(9,12*(m-1)) mat_Volt]}';
                            matlabbatch{3}.spm.stats.con.consess{1+4*(m-1)}.fcon.sessrep = 'replna';
                            %F contrast for canonical + derivaties: (N+1) contrasts for N sessions
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.name = names{m};
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.convec = {[temp_mat1 zeros(9,18*(m-1)) mat_Volt]}';
                            %{[zeros(9,12*(m-1)) mat_Volt]}';
                            matlabbatch{3}.spm.stats.con.consess{2+4*(m-1)}.fcon.sessrep = 'both';
                            %T contrast for activations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.name = [names{m} 'Activation'];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) 1];
                            matlabbatch{3}.spm.stats.con.consess{3+4*(m-1)}.tcon.sessrep = 'both';
                            %T contrast for deactivations: (N+1) contrasts
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.name = [names{m} 'Deactivation'];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) -1];
                            matlabbatch{3}.spm.stats.con.consess{4+4*(m-1)}.tcon.sessrep = 'both';
                        end
                    end
                    %Gives one contrast for movement, the last contrast
                    %NOT CLEAR WHY this contrast is not working
                    matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.name = 'MVT';
                    matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.convec = {[zeros(6,3*size(names,2) + ...
                        9*size(names,2)*(size(names,2)+1)/2) eye(6)]}';
                    matlabbatch{3}.spm.stats.con.consess{1+4*size(names,2)}.fcon.sessrep = 'both';
                    matlabbatch{3}.spm.stats.con.delete = 0;

                    %Contrasts for combining all spike types into one 
                    if size(names,2) > 1
                        convec_temp = []; convec_tempA = []; convec_tempD = [];
                        for m=1:size(names,2)
                            convec_temp = [convec_temp mat_Volt];
                            convec_tempA = [convec_tempA 1 zeros(1,11)];
                            convec_tempD = [convec_tempD -1 zeros(1,11)];
                        end
                        %F contrast on all spike types [eye(3) eye(3) ...]
                        matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.name = 'All spikes combined No Avg';                 
                        matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.convec = {convec_temp}';
                        matlabbatch{3}.spm.stats.con.consess{2+4*size(names,2)}.fcon.sessrep = 'replna';

                        matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.name = 'All spikes combined';                 
                        matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.convec = {convec_temp}';
                        matlabbatch{3}.spm.stats.con.consess{3+4*size(names,2)}.fcon.sessrep = 'both';

                        %T contrast for activations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.name = 'All spikes combined Activation';
                        matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.convec = convec_tempA;
                        matlabbatch{3}.spm.stats.con.consess{4+4*size(names,2)}.tcon.sessrep = 'both';
                        %T contrast for deactivations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.name = 'All spikes combined Deactivation';
                        matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.convec = convec_tempD;
                        matlabbatch{3}.spm.stats.con.consess{5+4*size(names,2)}.tcon.sessrep = 'both';                                
                    end %end if size(names,2) > 1
                    if size(names,2) == 2
                        %cross-stimulus interactions
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.fcon.name = 'Cross-interactions';                 
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.fcon.convec = {[zeros(9,15) eye(9)]}';
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.fcon.sessrep = 'replna';
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.fcon.name = 'Cross-interactions';                 
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.fcon.convec = {[zeros(9,15) eye(9)]}';
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.fcon.sessrep = 'both';

                    end
                else
                    %if there is only one session, the second group of
                    %contrasts disappears
                    %Get 3 contrasts for each onset type
                    for m=1:size(names,2)
                        %F contrast for canonical + derivaties: Get one
                        %contrast
                        matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.name = [names{m}];
                        temp_mat1 = [zeros(3,3*(m-1)) eye(3) zeros(3,3*(size(names,2)-m)); ...
                            zeros(6,3*size(names,2))]; %this is 9 x (3*size(names,2))
                        %this will not work if size(names,2) > 2: need to
                        %locate in the design matrix the interactions of stimulus 2 with
                        %stimulus 3
                        matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.convec = ...
                            {[temp_mat1 zeros(9,18*(m-1)) mat_Volt]}';                 
                        matlabbatch{3}.spm.stats.con.consess{1+3*(m-1)}.fcon.sessrep = 'replna';
                        %T contrast for activations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.name = [names{m} 'Activation'];
                        matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) 1];
                        matlabbatch{3}.spm.stats.con.consess{2+3*(m-1)}.tcon.sessrep = 'both';
                        %T contrast for deactivations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.name = [names{m} 'Deactivation'];
                        matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.convec = [zeros(1,3*(m-1)) -1];
                        matlabbatch{3}.spm.stats.con.consess{3+3*(m-1)}.tcon.sessrep = 'both';
                    end
                    matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.name = 'MVT';
                    matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.convec = {[zeros(6,3*size(names,2) + ...
                        9*size(names,2)*(size(names,2)+1)/2) eye(6)]}';
                    matlabbatch{3}.spm.stats.con.consess{4+3*(size(names,2)-1)}.fcon.sessrep = 'both';

                    %Contrasts for combining all spike types into one 
                    if size(names,2) > 1
                        convec_temp = []; convec_tempA = []; convec_tempD = [];
                        for m=1:size(names,2)                          
                            convec_temp = [convec_temp eye(3)];
                            convec_tempA = [convec_tempA 1 zeros(1,3)];
                            convec_tempD = [convec_tempD -1 zeros(1,3)];
                        end
                        convec_temp = [convec_temp; zeros(6,3*size(names,2))];
                        convec_temp2 = [];
                        switch size(names,2)
                            case 1
                                convec_temp2 = [convec_temp2 mat_Volt];
                            case 2
                                convec_temp2 = [convec_temp2 mat_Volt zeros(9,9) mat_Volt];
                        end
                                                             
                        %F contrast on all spike types [eye(3) eye(3) ...]
                        matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.name = 'All spikes combined No Avg';                 
                        matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.convec = {[convec_temp convec_temp2]}';
                        matlabbatch{3}.spm.stats.con.consess{5+3*(size(names,2)-1)}.fcon.sessrep = 'replna';
                        %T contrast for activations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.name = 'All spikes combined Activation';
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.convec = convec_tempA;
                        matlabbatch{3}.spm.stats.con.consess{6+3*(size(names,2)-1)}.tcon.sessrep = 'both';
                        %T contrast for deactivations: 1 contrast
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.name = 'All spikes combined Deactivation';
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.convec = convec_tempD;
                        matlabbatch{3}.spm.stats.con.consess{7+3*(size(names,2)-1)}.tcon.sessrep = 'both';                                                                          
                    end %end if size(names,2) > 1
                    if size(names,2) == 2
                        %cross-stimulus interactions
                        matlabbatch{3}.spm.stats.con.consess{8+3*(size(names,2)-1)}.fcon.name = 'Cross-interactions';                 
                        matlabbatch{3}.spm.stats.con.consess{8+3*(size(names,2)-1)}.fcon.convec = {[zeros(9,15) eye(9)]}';
                        matlabbatch{3}.spm.stats.con.consess{8+3*(size(names,2)-1)}.fcon.sessrep = 'replna';
                    end

                end %if size(f.fEPI,2) > 1

                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,['Volterra Stats failed to run for ' DirOut]);
            end 
        end %end if exist('SPM.mat','file')
    end %end if dStats  
    
    if dStats
    if run_report
        cd(DirStats);
        if exist('SPM.mat','file')
            spmmat = spm_select('List',pwd,'^SPM.mat'); 
            try
                clear matlabbatch
                    %Results report
                    matlabbatch{1}.spm.stats.results.spmmat = {spmmat};
                    load(f.fOnset{1});

                    mm = 0;
                    for m=1:size(names,2)
                        mm = mm + 1;
                        if size(f.fEPI,2) > 1
                            %contrast number for F-stat of type No Averaging
                            cn = 1 + (3 * size(f.fEPI,2) + 4) * (m-1);
                            %contrast number for T activation of type Averaging
                            cna = 3 + 2*size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
                            %contrast number for T deactivation of type Averaging
                            cnd = 4 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
                            %contrast for movement  
                            %cnm = ;
                            %contrast for cross interactions
                            %ci = ;
                        else
                            cn  = 1 + 3 * (m-1);
                            cna = 2 + 3 * (m-1);
                            cnd = 3 + 3 * (m-1);
                            %cnm = 4 + 3 * (m-1);
                            %ci = 
                        end
                        %Masked activations
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        %Masked deactivations
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                        mm = mm + 1;
                        matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation'];
                        matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                        matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                        matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                        matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                        matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
                    end %end for m=1:size(names,2)
    %                 
    %                 if size(f.fEPI,2) > 1
    %                         %contrast number for F-stat of type No Averaging
    %                         %cn = 1 + (3 * size(f.fEPI,2) + 4) * (m-1);
    %                         %contrast number for T activation of type Averaging
    %                         %cna = 3 + size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
    %                         %contrast number for T deactivation of type Averaging
    %                         %cnd = 4 + 2 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (m-1);
    %                         %contrast for movement  
    %                         cnm = 5 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * size(names,2);
    %                         %contrast for cross interactions
    %                         ci = 6 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * size(names,2);
    %                     else
    %                         %cn  = 1 + 3 * (m-1);
    %                         %cna = 2 + 3 * (m-1);
    %                         %cnd = 3 + 3 * (m-1);
    %                         cnm = 4 + 3 * size(names,2);
    %                         ci = 5 + 3 * size(names,2);
    %                 end
    %                     
    %                 %Cross interactions
    %                 if size(names,2) == 2
    %                     mm = mm + 1;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = ['Cross-interactions FWE Activation'];
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = ci;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = ci;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
    %                     mm = mm + 1;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = ['Cross-interactions unc.'];
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = ci;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
    %                     mm = mm + 1;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = ['Cross-interactions FWE'];
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = ci;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
    %                     mm = mm + 1;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = ['Cross-interactions unc.'];
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = ci;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
    %                     matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;
    %                 end
    %                 if reportContrastsBySession
    %                     %to do... 
    %                 end
                %Movement
    %             if size(f.fEPI,2) > 1                   
    %                 %contrast for movement
    %                 cnm = 6 + 3 * size(f.fEPI,2) + (3 * size(f.fEPI,2) + 4) * (size(names,2)-1);
    %             else
    %                 cnm = 4 + 3 * (size(names,2)-1); 
    %             end
    %             mm = mm + 1;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = 'FWE Movement';
    %             matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnm;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
    %             matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
    %             mm = mm + 1;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = 'unc. Movement';
    %             matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cnm;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
    %             matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
    %             matlabbatch{1}.spm.stats.results.conspec(mm).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
    % 
    %             matlabbatch{1}.spm.stats.results.units = 1;
    %             matlabbatch{1}.spm.stats.results.print = true;
    %             %end specification of Results Report
    % 
                        if reportContrastsBySession && size(f.fEPI,2) > 1
                        %Masked activations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked activations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked activations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Activation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked activations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Activation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cna;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        
                        %Masked deactivations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked deactivations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 10;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked deactivations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' FWE F-Deactivation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'FWE';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        %Masked deactivations
                        for j=1:size(f.fEPI,2)
                            for m=1:size(names,2)
                                %contrast number for F-stat for session j and contrast m
                                cn = 1 + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T activation 
                                cna = 2 + size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);
                                %contrast number for T deactivation 
                                cnd = 3 + 2 * size(f.fEPI,2) + j + (3 * size(f.fEPI,2) + 4) * (m-1);                                                                 
                                mm = mm + 1;
                                matlabbatch{1}.spm.stats.results.conspec(mm).titlestr = [names{m} ' unc. F-Deactivation Session ' int2str(j)];
                                matlabbatch{1}.spm.stats.results.conspec(mm).contrasts = cn;
                                matlabbatch{1}.spm.stats.results.conspec(mm).threshdesc = 'none';
                                matlabbatch{1}.spm.stats.results.conspec(mm).thresh = 0.001;
                                matlabbatch{1}.spm.stats.results.conspec(mm).extent = 0;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.contrasts = cnd;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.thresh = 0.05;
                                matlabbatch{1}.spm.stats.results.conspec(mm).mask.mtype = 0;                               
                            end %end for m
                        end %end for j
                        
                    end %end if reportContrastsBySession && size(f.fEPI,2) > 1
                spm_jobman('run',matlabbatch);
            catch
                write_log(flog,'Results Report failed to run');
            end   

            %end of stats - copy Results Report
            try
                temp_file = spm_select('List',pwd,'^spm_');      
                for nf=1:size(temp_file,1)
                    [~, fil ext] = fileparts(temp_file(nf,:));
                    if exist(temp_file(nf,:),'file')
                        copyfile(temp_file(nf,:),fullfile(DirResultsAll,[f.DirSubj '_' DirOut '_' fil ext]));
                    end
                 end
            catch 
                    write_log(flog,['Could not copy ' DirOut ' Results Report']);
            end
        end %end if exist('SPM.mat','file')
    end
    end
end

