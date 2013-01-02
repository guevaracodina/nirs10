function out = nirs_run_AnalyzerOnsets(job)
%Reads Markers file exported by Analyzer, to extract 'multiple conditions'
%data required by SPM: names, onsets, durations
%Returns TR, the EPI repetition time, and the name of the new onset file
outNIRSmat = {};
try
    outNIRSmat = job.NIRSmat;
    if ~isempty(outNIRSmat{1})
        nsubj = length(outNIRSmat);
        NIRSok = 1;
    else
        nsubj = 1;
        NIRSok = 0;
    end
catch
    NIRSok = 0;
    nsubj = 1;
end
try
    t = job.raw_onset_files; %spm_select([1 Inf],'any','Select files of onsets');
    if size(t,1) == 1
        if isempty(t{1,1})
            new_onset_files = 0;
            dp_NIRS1 = [];
            freq_NIRS1 = [];
        end
    else
        new_onset_files = 1;
        try
            dp_NIRS1 = job.dp_NIRS1{1,1};
            freq_NIRS1 = job.freq_NIRS1{1,1};
        catch
            dp_NIRS1 = [];
            freq_NIRS1 = 19.5312;
        end
    end
catch
    new_onset_files = 0;
    dp_NIRS1 = [];
    freq_NIRS1 = [];
end
if isfield(job,'onset_to_keep')
    onset_to_keep = job.onset_to_keep;
else
    onset_to_keep = '';
end
%res_factor = 10; %factor to control precision of interpolation for confound regressors
%Code allows up to 3 onset types to be removed. If want to keep all onsets,
%make sure all rem_onsets are set to ''
rem_onsets{1} = 'S250'; %'spkLFT';
rem_onsets{2} = 'S252';
rem_onsets{3} = 'S254';
rem_onsets{4} = 'F 0';
rem_onsets{5} = 'sCS';
rem_onsets{6} = 'eCS';
rem_onsets{7} = 'GSWD L>R';
rem_onsets{8} = 'burst start';
rem_onsets{9} = '250'; %'spkLFT';
rem_onsets{10} = '252';
rem_onsets{11} = '254';
%Exit if nothing available
if ~new_onset_files && ~NIRSok, return; end
%Exit if more than one subject and new onset files were specified
if nsubj > 1 && new_onset_files, return; end

%**********************************************************************************************
%Ke Peng, introduce heart rate repair process
%17/04/2012

if isfield(job.cardiac_repair,'cardiac_repair_on')
    cardiac_repair_on = 1;
    avg_number = job.cardiac_repair.cardiac_repair_on.avg_number;
    gap_def = job.cardiac_repair.cardiac_repair_on.gap_def;
else
    cardiac_repair_on = 0;
end
%***********************************************************************************************

%list of NIRS.mat locations, one per subject
for Idx=1:nsubj
    %Load NIRS.mat information
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'A2OnsetsOK') || job.force_redo)
            
            if NIRSok
                %load(job.NIRSmat{Idx,1});
                if isempty(freq_NIRS1)
                    freq_NIRS1 = NIRS.Cf.dev.fs;
                end
                NC = NIRS.Cf.H.C.N;
                %use last step of preprocessing
                lst = length(NIRS.Dt.fir.pp);
                rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
            end
            if ~new_onset_files
                t = NIRS.Dt.fir.rons;
            end
            for i3=1:size(t,1)
                clear names onsets durations
                [dir1 fil1 ext1] = fileparts(t{i3});
                if strcmp(ext1,'.mat')
                    load(t{i3});
                else
                    if NIRSok
                        d = fopen_NIR(rDtp{i3},NC);
                        dp_NIRS1 = size(d,2);
                        clear d
                    end
                    
                    
                    file = fullfile(dir1, [fil1 ext1]);
                    fp = fopen(file,'r');
                    
                    %Extract frequency information from first line
                    temp = fgetl(fp);
                    i1 = findstr(':',temp); i2 = findstr('ms',temp);
                    %Delta t between points in seconds
                    dt = str2double(temp(i1(2)+1:i2-1))/1000;
                    
                    %find first 255 marker -- ignore onsets that would arrive before the
                    %first Scan Start
                    while ~feof(fp)
                        temp = fgetl(fp);
                        indx = findstr('Scan Start',temp);
                        if isempty(indx)
                            indx = findstr('S255',temp);
                        end
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
                                && ~strcmpi(tag,'Bad Interval') && ~strncmp(tag,'B',1) && ~strncmpi(tag,'CS',2) ...
                                && (~strncmp(tag,'S',1) || strncmpi(tag,'Spk',3) || strncmpi(tag,'SWD',3)) && ~strncmp(tag,'T',1)...
                                && (~strcmp(tag,rem_onsets{1})) ...
                                && (~strcmp(tag,rem_onsets{2})) ...
                                && (~strcmp(tag,rem_onsets{3})) ...
                                && (~strcmp(tag,rem_onsets{4})) ...
                                && (~strcmp(tag,rem_onsets{5})) ...
                                && (~strcmp(tag,rem_onsets{6})) ...
                                && (~strcmp(tag,rem_onsets{7})) ...
                                && (~strcmp(tag,rem_onsets{8})) ...
                                && (~strcmp(tag,rem_onsets{9})) ...
                                && (~strcmp(tag,rem_onsets{10})) ...
                                && (~strcmp(tag,rem_onsets{11}))
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
                                if isempty(onset_to_keep) || strcmp(tag,onset_to_keep)
                                    %new type of onset
                                    NTonset = NTonset + 1;
                                    k{NTonset} = 1;
                                    names{NTonset} = tag;
                                    onsets{NTonset}(k{NTonset}) = str2double(temp(indx(2):indx(3)))*dt-start_scan;
                                    durations{NTonset}(k{NTonset}) = dt*(str2double(temp(indx(3):indx(4)))-1); %subtract 1 so point onsets have zero duration
                                end
                            end
                        else
                            %Look for movement or pulse markers
                            if strncmp(tag,'R',1) %pulse
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
                        end
                    end
                    
                    fclose(fp);
                    
                    %Save as .mat
                    %fOnset = fullfile(dir1,[fil1 '.mat']);
                    %No longer required to save onsets separately - better
                    %accounted for in NIRS.mat
                    %save(fOnset,'names','onsets','durations')
                    %write heart R peaks
                    %fRpeaks = fullfile(dir1,[fil1 '_Rpeaks.mat']);
                    %save(fRpeaks,'vR')
                    try
                        %Calculate pulse regressor
                        %if ~isempty(dp_NIRS1) && ~isempty(freq_NIRS1)
                        lpi = linspace(0,dp_NIRS1/freq_NIRS1,dp_NIRS1); %*res_factor);
                        
                        dvR = diff(vR);
                        
                        %**********************************************************************************************
                        %Ke Peng, introduce heart rate repair process
                        %17/04/2012
                        
                        if cardiac_repair_on
                            try
                                [dvR,vR,Flag_Repair,TotalAdded_Repair] = HeartRate_Repair(dvR,vR,avg_number,gap_def);
                                disp(['Flag_Repair = ' num2str(Flag_Repair) ' for session: ' num2str(i3)]);
                                disp(['TotalAdded_Repair = ' num2str(TotalAdded_Repair) ' for session: ' num2str(i3)]);
                                iR = iR + TotalAdded_Repair;
                            catch
                                disp('Cardiac repair failed.');
                            end
                        end
                        %**********************************************************************************************
                        
                        %vRi = interp1(vR(2:end),dvR,lpi,'linear');
                        vRi = interp1(vR(2:end),dvR,lpi,'pchip');
                        
                        %derivative of pulse rate - too noisy
                        %ddvR = diff(dvR);
                        %vRi2 = interp1(vR(2:end-1),ddvR,lpi,'spline');
                        
                        %Calculate movement regressors - to be done, not clear precisely what
                        %to do
                        
                        %Save
                        %fRegress = fullfile(dir1,['rp_' fil1 '.txt']);
                        %subtract the mean
                        %vRi = vRi - mean(vRi); %PP not required: will be done by spm_detrend
                        fR = 60./vRi'; % vRi2']; %cardiac frequency in heart beats per minute, as a column vector
                        h = figure; plot(lpi,fR); 
                        title([NIRS.Dt.s.p ': Heart rate (BPM), Sess' int2str(i3)],'interpreter','none');
                        [dirNIRS fil0] = fileparts(newNIRSlocation);
                        fname = fullfile(dirNIRS, ['HeartRate_Sess' int2str(i3)]); 
                        %warning('off')
                        print(h,'-dpng',[fname '.png'],'-r300');
                        saveas(h,[fname '.fig']);
                        %warning('on')
                        close(h);
                        %save(fRegress,'fR','-ascii');
                        %*****************************************************
                        %Ke Peng, 12/03/2012
                        %*****************************************************
                        
                        %[fR,Flag_Repair,TotalAdded_Repair] = HeartRate_Repair(fR,5,2.5);
                        %disp(['Flag_Repair = ' num2str(Flag_Repair) ' for session: ' num2str(i3)]);
                        %disp(['TotalAdded_Repair = ' num2str(TotalAdded_Repair) ' for session: ' num2str(i3)]);
                        %iR = iR + TotalAdded_Repair;
                        
                    catch
                        disp('Analyze onsets error');
                    end
                end
                if NIRSok
                    U = [];
                    P.name = 'none';
                    P.h    = 0;
                    for kk=1:length(names)
                        U(kk).name = names(kk);
                        U(kk).ons = onsets{kk};
                        U(kk).dur = durations{kk};
                        U(kk).P = P;
                    end
                    NIRS.Dt.fir.Sess(i3).U = U;
                    try
                        NIRS.Dt.fir.Sess(i3).vR = {vR};
                        NIRS.Dt.fir.Sess(i3).fR = {fR};
                    end
                    NIRS.flags.A2OnsetsOK = 1;
                    save(job.NIRSmat{Idx,1},'NIRS');
                else
                    outfile = fullfile(dir1, [fil1 '.mat']);
                    try
                        save(outfile, 'names' , 'onsets', 'durations', 'fR', 'vR');
                    catch
                        save(outfile, 'names' , 'onsets', 'durations');
                    end
                end
            end
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not process onsets for subj' int2str(Idx)]);
    end
    
end
out.NIRSmat = outNIRSmat;