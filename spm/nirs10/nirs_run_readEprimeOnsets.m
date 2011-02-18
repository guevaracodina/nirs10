function out = nirs_run_readEprimeOnsets(job)
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% 02/2011
% Retrieving data from e-prime (already transfered to excel format)
% Read the timings of important blocs

dir = job.dir;

% Getting the MRI sending trigger 5, the beggining of onsets with their
% tags
[Trig5MRI, headerTrig] = xlsread('904-17111.xlsx',1,'AD1:AD90');
[StimOnset, headerStim] = xlsread('904-17111.xlsx',1,'BU1:BU90');
[StimTags, headerTag] = xlsread('904-17111.xlsx',1,'BY1:BY90');
OnsetTime = StimOnset - Trig5MRI;
ControlOnsetTime = OnsetTime(1:10);
ConcretStimOnset=[];
P_ConcretStimOnset=[];
AbsStimOnset=[];
P_AbsStimOnset=[];
for i=1:length(StimTags)
    if StimTags(i) == 1
       ConcretStimOnset = [ConcretStimOnset OnsetTime(i+10)]; 
    elseif StimTags(i) == 2
        P_ConcretStimOnset = [P_ConcretStimOnset OnsetTime(i+10)];
    elseif StimTags(i) == 3
        AbsStimOnset = [AbsStimOnset OnsetTime(i+10)];
    elseif StimTags(i) == 4
        P_AbsStimOnset = [P_AbsStimOnset OnsetTime(i+10)];
    end
end
onsets={};
onsets{1,1}=ControlOnsetTime./1000;
onsets{1,2}=ConcretStimOnset./1000;
onsets{1,3}=P_ConcretStimOnset./1000;
onsets{1,4}=AbsStimOnset./1000;
onsets{1,5}=P_AbsStimOnset./1000;

durations={};
for d=1:5
durations{1,d}=1.36;
end

names={};
names{1,1}='Control';
names{1,2}='Concret';
names{1,3}='Pseudo_Concret';
names{1,4}='Abstract';
names{1,5}='Pseudo_Abstract';

save Onsets_S1 names onsets durations

%% Constructing the onset.mat for the 2nd round
% Read the timings of important blocs
clear all
[Trig5MRI, headerTrig] = xlsread('904-17222.xlsx',1,'AD1:AD90');
[StimOnset, headerStim] = xlsread('904-17222.xlsx',1,'BU1:BU90');
[StimTags, headerTag] = xlsread('904-17222.xlsx',1,'BY1:BY90');
OnsetTime = StimOnset - Trig5MRI;
ControlOnsetTime = OnsetTime(1:10);
ConcretStimOnset=[];
P_ConcretStimOnset=[];
AbsStimOnset=[];
P_AbsStimOnset=[];
for i=1:length(StimTags)
    if StimTags(i) == 1
       ConcretStimOnset = [ConcretStimOnset OnsetTime(i+10)]; 
    elseif StimTags(i) == 2
        P_ConcretStimOnset = [P_ConcretStimOnset OnsetTime(i+10)];
    elseif StimTags(i) == 3
        AbsStimOnset = [AbsStimOnset OnsetTime(i+10)];
    elseif StimTags(i) == 4
        P_AbsStimOnset = [P_AbsStimOnset OnsetTime(i+10)];
    end
end

onsets={};
onsets{1,1}=ControlOnsetTime./1000;
onsets{1,2}=ConcretStimOnset./1000;
onsets{1,3}=P_ConcretStimOnset./1000;
onsets{1,4}=AbsStimOnset./1000;
onsets{1,5}=P_AbsStimOnset./1000;

durations={};
for d=1:5
durations{1,d}=1.36;
end

names={};
names{1,1}='Control';
names{1,2}='Concret';
names{1,3}='Pseudo_Concret';
names{1,4}='Abstract';
names{1,5}='Pseudo_Abstract';

save Onsets_S2 names onsets durations

%% Round 3 of the same scan
% Read the timings of important blocs
clear all
[Trig5MRI, headerTrig] = xlsread('904-17333.xlsx',1,'AD1:AD90');
[StimOnset, headerStim] = xlsread('904-17333.xlsx',1,'BU1:BU90');
StimTags = xlsread('904-17333.xlsx',1,'BY1:BY90');
OnsetTime = StimOnset - Trig5MRI;
ControlOnsetTime = OnsetTime(1:10);
ConcretStimOnset=[];
P_ConcretStimOnset=[];
AbsStimOnset=[];
P_AbsStimOnset=[];
for i=1:length(StimTags)
    if StimTags(i) == 1
       ConcretStimOnset = [ConcretStimOnset OnsetTime(i)]; 
    elseif StimTags(i) == 2
        P_ConcretStimOnset = [P_ConcretStimOnset OnsetTime(i)];
    elseif StimTags(i) == 3
        AbsStimOnset = [AbsStimOnset OnsetTime(i)];
    elseif StimTags(i) == 4
        P_AbsStimOnset = [P_AbsStimOnset OnsetTime(i)];
    end
end

onsets={};
onsets{1,1}=ControlOnsetTime./1000;
onsets{1,2}=ConcretStimOnset./1000;
onsets{1,3}=P_ConcretStimOnset./1000;
onsets{1,4}=AbsStimOnset./1000;
onsets{1,5}=P_AbsStimOnset./1000;

durations={};
for d=1:5
durations{1,d}=1.36;
end

names={};
names{1,1}='Control';
names{1,2}='Concret';
names{1,3}='Pseudo_Concret';
names{1,4}='Abstract';
names{1,5}='Pseudo_Abstract';

save Onsets_S3 names onsets durations

end