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

% 
% % % to calculate RT and accuracy from excel ePrime files
% % to get the respons culmns 
% clear all;
% cd('\\recherche2\documents$\amirim\Documents\My Project\Pilote; functional localizer\Excel Files')
% RT_R1 = xlsread('904-17111.xlsx',1,'BX12:BX90');
% lateRT_R1 = xlsread('904-17111.xlsx',1,'CG12:CG90');
% StimOnsetTimeR1 = xlsread('904-17111.xlsx',1,'BU12:BU90');
% StimTagsR1 = xlsread('904-17111.xlsx',1,'BY12:BY90');
% lateRT_R1 = lateRT_R1(2:end);
% lateRT_R1 = [lateRT_R1; 0];
% RTs_R1 = RT_R1 + lateRT_R1;
% RTTime_R1 = RTs_R1 - StimOnsetTimeR1;
% % RT1_Round2 = xlsread('903-16222.xlsx',-1);
% % RT1_Round3 = xlsread('903-16333.xlsx',-1);
% RT_R2 = xlsread('904-17222.xlsx',1,'BX12:BX90');
% lateRT_R2 = xlsread('904-17222.xlsx',1,'CG12:CG90');
% StimOnsetTimeR2 = xlsread('904-17222.xlsx',1,'BU12:BU90');
% StimTagsR2 = xlsread('904-17222.xlsx',1,'BY12:BY90');
% lateRT_R2 = lateRT_R2(2:end);
% lateRT_R2 = [lateRT_R2; 0];
% RTs_R2 = RT_R2 + lateRT_R2;
% RTTime_R2 = RTs_R2 - StimOnsetTimeR2;
% % RTs for round3
% % RT_R3 = xlsread('903-16333.xlsx',1,'BI2:BI80');
% % lateRT_R3 = xlsread('903-16333.xlsx',1,'BR2:BR80');
% % StimOnsetTimeR3 = xlsread('903-16333.xlsx',1,'BF2:BF80');
% % StimTagsR3 = xlsread('903-16333.xlsx',1,'BJ2:BJ80');
% RT_R3 = xlsread('904-17333.xlsx',1,'BX12:BX90');
% lateRT_R3 = xlsread('904-17333.xlsx',1,'CG12:CG90');
% StimOnsetTimeR3 = xlsread('904-17333.xlsx',1,'BU12:BU90');
% StimTagsR3 = xlsread('904-17333.xlsx',1,'BY12:BY90');
% lateRT_R3 = lateRT_R3(2:end);
% lateRT_R3 = [lateRT_R3; 0];
% RTs_R3 = RT_R3 + lateRT_R3;
% RTTime_R3 = RTs_R3 - StimOnsetTimeR3;
% 
% TAGs = [StimTagsR1; StimTagsR2; StimTagsR3];
% RTs = [RTTime_R1; RTTime_R2; RTTime_R3];
% save Pilote904 TAGs RTs