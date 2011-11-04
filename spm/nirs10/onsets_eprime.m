function out = onsets_eprime(job)
NIRSmat = job.NIRSmat;
excelp = job.excelp;
eprime = job.eprime;

% To find the total number of sheets/rouns per excel file
[typ, desc] = xlsfinfo(excelp);

% Read columns already chosen from the readEprime UI
% loop over sessions
for iR=1:length(desc)
    % Retrieving necessary data from excelEprime sheet
    [data, header] = xlsread(excelp, iR);

    eprime(iR).col.trig.clp = xlsColNum2Str(find(strcmp(eprime(1).col.trig.n, header(1,:))));
    eprime(iR).col.stim.clp = xlsColNum2Str(find(strcmp(eprime(1).col.stim.n, header(1,:))));
    eprime(iR).col.cond.clp = xlsColNum2Str(find(strcmp(eprime(1).col.cond.n, header(1,:))));
%     save(fullfile(sDtp,'NIRS.mat'),'eprime');
    
    TrigOnsetTime = [];
    StimOnsetTime =[];
    StimTags = [];
    
    Trigc = eprime(iR).col.trig.clp{1,1};
    ttl = [Trigc '2'];
    TrigOnsetTime = xlsread(excelp,iR,ttl);
    
    Stimc = eprime(iR).col.stim.clp{1,1};
    s = Stimc;
    Scl = [s ':' s];
    StimOnsetTime = xlsread(excelp,iR,Scl);
    
    Condc = eprime(iR).col.cond.clp{1,1};
    c = Condc;
    Ccl = [c ':' c];
    StimTags = xlsread(excelp,iR,Ccl);
    
    % Initialize timing by removing the delay time before start point sent
    % by trigger
    OnsetTime = StimOnsetTime - TrigOnsetTime;
    
    % Number of conditions (stimulus tag) of the protocol
    durations={};
    nC = max(StimTags) + 1;
    D = 1.36;  % The presentation time for each stimulus set in e-prime
    for d=1:nC
        durations{1,d} = D;
    end
    names={};
    names{1,1}='Control';
    names{1,2}='Concret';
    names{1,3}='Pseudo_Concret';
    names{1,4}='Abstract';
    names{1,5}='Pseudo_Abstract';
    
    % Getting the TTL trigger, the beggining of onsets with their tag
    %     ControlOnsetTime = [];
    if iR == 1
        j=20;
    ControlOnsetTime = OnsetTime(1:j);
    else
        j=0;
        ControlOnsetTime = [];
    end
    ConcretStimOnset=[];
    P_ConcretStimOnset=[];
    AbsStimOnset=[];
    P_AbsStimOnset=[];
    for iT=1:size(StimTags,1)
        switch StimTags(iT)
            case 1
                ConcretStimOnset = [ConcretStimOnset OnsetTime(iT+j)];
            case 2
                P_ConcretStimOnset = [P_ConcretStimOnset OnsetTime(iT+j)];
            case 3
                AbsStimOnset = [AbsStimOnset OnsetTime(iT+j)];
            case 4
                P_AbsStimOnset = [P_AbsStimOnset OnsetTime(iT+j)];
                %             otherwise
                %             ControlOnsetTime = [ControlOnsetTime OnsetTime(iT)];
        end
    end
    
    % Calculating offset
    % offset in sec. is the delay between lunching e-prime and cw6
    %             eprimep = fileparts(excelp);
    %             subjectp = fileparts(eprimep);
    %             jobO.nirsp = fullfile(subjectp, 'NIRS');
    jobO.NIRSmat = NIRSmat;
    jobO.iR =iR;
    offset = offset_nirs(jobO); % by default, we assume the delay between eprime and cw6 is zero
    %         offset = 0; %-9.6;
    onsets={};
    onsets{1,1}=ControlOnsetTime'./1000 + offset;
    onsets{1,2}=ConcretStimOnset./1000 + offset;
    onsets{1,3}=P_ConcretStimOnset./1000 + offset;
    onsets{1,4}=AbsStimOnset./1000 + offset;
    onsets{1,5}=P_AbsStimOnset./1000 + offset;
    
    if iR == 1;
        durations = durations;
        onsets = onsets;
        names = names;
    else
        durations = durations(1,2:5);
        onsets = onsets(1,2:5);
        names = names(1,2:5);
        ControlOnsetTime = [];
    end
    %     [nSubj, ext] = fileparts(FileName);
    %     Onsets{1,iR} = sprintf('Onsets_%s_%s', nSubj, desc{1,iR});
    %     Sess{1,iR}.name = names;
    %     Sess{1,iR}.ons = onsets;
    %     Sess{1,iR}.dur = durations;
    
    %Ignore parametric modulations - cf spm_run_fmri_design.m
    P.name = 'none';
    P.h    = 0;
    for kk = 1:size(names, 2)
        Sess(iR).U(kk).name = names(kk);
        Sess(iR).U(kk).ons = onsets{1,kk};
        Sess(iR).U(kk).dur = durations{kk};
        Sess(iR).U(kk).P = P;
    end
end

out = Sess;


