function nirs_testStimuli_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS - test stimuli

global nirs10
  
nirs10.readOnsets.addTestStimuli.testStimulusName = 'Test';
nirs10.readOnsets.addTestStimuli.keepChannels = 1:2;
nirs10.readOnsets.addTestStimuli.testChannels = 1;
nirs10.readOnsets.addTestStimuli.testParadigmType =  1; 
nirs10.readOnsets.addTestStimuli.testStimuliNumber = 100;
nirs10.readOnsets.addTestStimuli.testSessionNumber = 1;
nirs10.readOnsets.addTestStimuli.testWavelength = 2;
nirs10.readOnsets.addTestStimuli.testAmplitudeTarget.testAmplitude = 1;
nirs10.readOnsets.addTestStimuli.testAmplitudeTarget.testSNR = -15;
nirs10.readOnsets.addTestStimuli.testAmplitude2 = 50;
%NoFrequentSpikes
nirs10.readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testExpSlowSpike1 = 20;
nirs10.readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testRescaleOn1 = 1;
nirs10.readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testSeed1 = 0;
%With FrequentSpikes
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpSlowSpike2 = 20;
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpFastSpike = [0.5 1.5]; %min and max
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumSlowSpikes_perGroup = 3;
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumFastSpikes_perGroup = 5;
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testRescaleOn2 = 2;
nirs10.readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testSeed2 = 0;
