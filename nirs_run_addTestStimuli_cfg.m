function addTestStimuli = nirs_run_addTestStimuli_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'addStim');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests: add stimuli and corresponding HRFs to raw NIRS signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testDupChannels         = cfg_menu;
testDupChannels.tag     = 'testDupChannels';
testDupChannels.name    = 'Duplicate number of channels';
testDupChannels.help    = {'For each specified channel, add stimuli.'
    'and also add a copy of that channel, without stimuli.'}';
testDupChannels.labels = {
    'True'
    'False'
    }';
testDupChannels.values = {1, 0};
testDupChannels.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testDupChannels', val{:});

testStimulusName         = cfg_entry;
testStimulusName.name    = 'Name of the stimulus';
testStimulusName.tag     = 'testStimulusName';
testStimulusName.strtype = 's';
testStimulusName.num     = [1 Inf];
testStimulusName.def  = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testStimulusName', val{:});
testStimulusName.help    = {'Name of the stimulus.'};

keepChannels         = cfg_entry;
keepChannels.name    = 'List of channels to keep';
keepChannels.tag     = 'keepChannels';
keepChannels.strtype = 'r';
keepChannels.num     = [1 Inf];
keepChannels.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.keepChannels', val{:});
keepChannels.help    = {'Enter channel numbers to keep.'};

AllChannels           = cfg_branch;
AllChannels.name      = 'Keep ALL channels';
AllChannels.tag       = 'AllChannels';
AllChannels.help      = {'Keep all channels'};

keepAllChannels        = cfg_choice;
keepAllChannels.name   = 'Specify channels to keep';
keepAllChannels.tag    = 'keepAllChannels';
keepAllChannels.values = {AllChannels keepChannels};
%Do not know how to specify a default value by a call using .def for a
%cfg_choice object
keepAllChannels.val    = {keepChannels};
keepAllChannels.help   = {'Choose whether to keep all channels or select a subset.'};

%Which channels to test
testChannels         = cfg_entry;
testChannels.name    = 'Test channel numbers';
testChannels.tag     = 'testChannels';
testChannels.strtype = 'r';
testChannels.num     = [1 Inf];
testChannels.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.testChannels', val{:});
testChannels.help    = {'Enter channel numbers where the stimuli will be added.'
    'It is sufficient to give only the channel numbers for the first '
    'wavelength, and the code will add the mirror channels for the other wavelengths'}';

%Test stimuli number
testStimuliNumber         = cfg_entry;
testStimuliNumber.name    = 'Test stimuli number';
testStimuliNumber.tag     = 'testStimuliNumber';
testStimuliNumber.strtype = 'r';
testStimuliNumber.num     = [1 1];
testStimuliNumber.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testStimuliNumber', val{:});
testStimuliNumber.help    = {'Enter number of stimuli to be added.'};

testBP           = cfg_branch; %empty branch
testBP.name      = 'Block Paradigm';
testBP.tag       = 'testBP';
testBP.val       = {};
testBP.help      = {'Block Paradigm'};

testSeed1         = cfg_entry;
testSeed1.name    = 'Random seed';
testSeed1.tag     = 'testSeed1';
testSeed1.strtype = 'r';
testSeed1.num     = [1 1];
testSeed1.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testSeed1', val{:});
testSeed1.help    = {'Enter a seed for the random number generator.'};

testSeed2         = cfg_entry;
testSeed2.name    = 'Random seed';
testSeed2.tag     = 'testSeed2';
testSeed2.strtype = 'r';
testSeed2.num     = [1 1];
testSeed2.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testSeed2', val{:});
testSeed2.help    = {'Enter a seed for the random number generator.'};

testExpFastSpike         = cfg_entry;
testExpFastSpike.name    = 'Time interval for frequent spikes';
testExpFastSpike.tag     = 'testExpFastSpike';
testExpFastSpike.strtype = 'r';
testExpFastSpike.num     = [1 2];
testExpFastSpike.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpFastSpike', val{:});
testExpFastSpike.help    = {'Enter boundary min and max in seconds,'
    'of the uniform distribution of time intervals,'
    'between frequent spikes.'}';

testExpSlowSpike1         = cfg_entry;
testExpSlowSpike1.name    = 'Parameter for infrequent spikes interval';
testExpSlowSpike1.tag     = 'testExpSlowSpike1';
testExpSlowSpike1.strtype = 'r';
testExpSlowSpike1.num     = [1 1];
testExpSlowSpike1.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.NoFrequentSpikes.testExpSlowSpike1', val{:});
testExpSlowSpike1.help    = {'Enter characteristic scale of infrequent spikes,'
    'in seconds'}';

testExpSlowSpike2         = cfg_entry;
testExpSlowSpike2.name    = 'Parameter for infrequent spikes interval';
testExpSlowSpike2.tag     = 'testExpSlowSpike2';
testExpSlowSpike2.strtype = 'r';
testExpSlowSpike2.num     = [1 1];
testExpSlowSpike2.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testExpSlowSpike2', val{:});
testExpSlowSpike2.help    = {'Enter characteristic scale of infrequent spikes,'
    'in seconds'}';

testAvgNumFastSpikes_perGroup         = cfg_entry;
testAvgNumFastSpikes_perGroup.name    = 'Avg number of frequent spikes per group';
testAvgNumFastSpikes_perGroup.tag     = 'testAvgNumFastSpikes_perGroup';
testAvgNumFastSpikes_perGroup.strtype = 'r';
testAvgNumFastSpikes_perGroup.num     = [1 1];
testAvgNumFastSpikes_perGroup.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumFastSpikes_perGroup', val{:});
testAvgNumFastSpikes_perGroup.help    = {
    'Follows a Poisson process: Enter typical number.'
    'of fast spikes in a discharge'}';

testAvgNumSlowSpikes_perGroup         = cfg_entry;
testAvgNumSlowSpikes_perGroup.name    = 'Average number of infrequent spikes per group';
testAvgNumSlowSpikes_perGroup.tag     = 'testAvgNumSlowSpikes_perGroup';
testAvgNumSlowSpikes_perGroup.strtype = 'r';
testAvgNumSlowSpikes_perGroup.num     = [1 1];
testAvgNumSlowSpikes_perGroup.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testPType.testEP.FrequentSpikes.testAvgNumSlowSpikes_perGroup', val{:});
testAvgNumSlowSpikes_perGroup.help    = {
    'Follows a Poisson process: Enter typical number.'
    'of infrequent spikes in a discharge'}';

testRescaleOn1         = cfg_menu;
testRescaleOn1.tag     = 'testRescaleOn1';
testRescaleOn1.name    = 'Rescale time to get specified number of spikes';
testRescaleOn1.help    = {'When rescaling, guaranteed to get the specified '
    'number of spikes. When not rescaling, the program will generate spikes '
    'until the specified number is reached, or the end of the dataset is reached, '
    'which could lead to a having fewer spikes than the number specified.'}';
testRescaleOn1.labels = {
    'True'
    'False'
    }';
testRescaleOn1.values = {1, 0};
testRescaleOn1.val = {1};

testRescaleOn2         = cfg_menu;
testRescaleOn2.tag     = 'testRescaleOn2';
testRescaleOn2.name    = 'Rescale time to get specified number of spikes';
testRescaleOn2.help    = {'When rescaling, guaranteed to get the specified '
    'number of spikes. When not rescaling, the program will generate spikes '
    'until the specified number is reached, or the end of the dataset is reached, '
    'which could lead to a having fewer spikes than the number specified.'}';
testRescaleOn2.labels = {
    'True'
    'False'
    }';
testRescaleOn2.values = {1, 0};
testRescaleOn2.val = {0};

NoFrequentSpikes           = cfg_branch;
NoFrequentSpikes.name      = 'No frequent spikes';
NoFrequentSpikes.tag       = 'NoFrequentSpikes';
NoFrequentSpikes.val       = {testExpSlowSpike1 testRescaleOn1 testSeed1};
NoFrequentSpikes.help      = {'No frequent spikes, only one type of spikes'};

FrequentSpikes           = cfg_branch;
FrequentSpikes.name      = 'Frequent spikes included';
FrequentSpikes.tag       = 'FrequentSpikes';
FrequentSpikes.val       = {testExpSlowSpike2 testExpFastSpike ...
    testAvgNumFastSpikes_perGroup testAvgNumSlowSpikes_perGroup testRescaleOn2 testSeed2};
FrequentSpikes.help      = {'Frequent spikes following a uniform distribution '
    'interspersed with slow spikes.'}';

%test Event paradigm
testEP           = cfg_choice;
testEP.name      = 'Event Paradigm';
testEP.tag       = 'testEP';
testEP.values    = {NoFrequentSpikes FrequentSpikes};
testEP.val       = {FrequentSpikes};
%testEP.def     = @(val)nirs_get_defaults('readOnsets.addTestStimuli.testPType', val{:});
testEP.help      = {'Event Paradigm'
    'Choose whether to have only one type of spikes,'
    'all following an exponential distribution, '
    'or to also include fast spikes, which follow '
    'a uniform distribution.'
    'With two types of spikes, model works as follows,'
    'alternating between infrequent and frequent'
    'bunches of spikes (with numbers specified by separate Poisson distributions) '
    'putting the intervals one after the other. '
    'a number of infrequent spikes is '
    'obtained from a Poisson distribution; for each, a time interval to the next'
    'spike is generated from an exponential distribution '
    'This is then repeated for fast spikes, until the total desired number of'
    'spikes is obtained. '}';

%Paradigm type
testPType      = cfg_choice;
testPType.tag  = 'testPType';
testPType.name = 'Paradigm: block or event';
testPType.values = {testBP testEP};
testPType.val = {testEP};
testPType.help = {'Choose test paradigm type: block (the number of stimuli '
    'will be evenly spread in time in the data file; or event (random stimuli '
    'generation, following a combination of uniform, exponential and Poisson distributions).'}';

%Test session number
testSessionNumber         = cfg_entry;
testSessionNumber.name    = 'Test session number';
testSessionNumber.tag     = 'testSessionNumber';
testSessionNumber.strtype = 'r';
testSessionNumber.num     = [1 1];
testSessionNumber.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testSessionNumber', val{:});
testSessionNumber.help    = {'Enter NIRS data session number where stimuli are to be added.'};

%Test wavelength
testWavelength         = cfg_entry;
testWavelength.name    = 'Test wavelength number(s)';
testWavelength.tag     = 'testWavelength';
testWavelength.strtype = 'r';
testWavelength.num     = [1 Inf];
testWavelength.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testWavelength', val{:});
testWavelength.help    = {'Enter wavelength number(s) of the data where the '
    'stimuli will be added. For example, enter 2 if OD at 690 nm is '
    'desired for the test and is the second wavelength.'}';

testAmplitude         = cfg_entry;
testAmplitude.name    = 'Test response amplitude';
testAmplitude.tag     = 'testAmplitude';
testAmplitude.strtype = 'r';
testAmplitude.num     = [1 Inf];
testAmplitude.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testAmplitudeTarget.testAmplitude', val{:});
testAmplitude.help    = {'Enter desired response amplitude as a percentage. '
    'This will be understood as a percentage of the median of the raw intensity '
    'signal. If several channels are tested, a vector of percentages can '
    'be entered, that will be applied channelwise.'
    'effective SNR will be calculated for each channel and protocole.'}';

testSNR         = cfg_entry;
testSNR.name    = 'Test SNR value';
testSNR.tag     = 'testSNR';
testSNR.strtype = 'r';
testSNR.num     = [1 Inf];
testSNR.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testAmplitudeTarget.testSNR', val{:});
testSNR.help    = {'Enter desired SNR in dB. '
    'The amplitude to be applied to each channel will be calculated '
    'Based on the formula SNR = 10 log10 (A * Ep/Eb),'
    'Where A is the amplitude, Eb is the power of the baseline (the sum'
    'of squares of amplitudes at each time point) and Ep is the power of'
    'the protocole after convolution with the HRFs.'
    'This option will result in a constant SNR value for each channel '
    'and each protocole, but the amplitude will vary.'
    'The calculated value for the amplitude will be stored.'}';

testAmplitudeTarget        = cfg_choice;
testAmplitudeTarget.name   = 'Target Amplitude or SNR';
testAmplitudeTarget.tag    = 'testAmplitudeTarget';
testAmplitudeTarget.values = {testAmplitude testSNR};
%Do not know how to specify a default value by a call using .def for a
%cfg_choice object
testAmplitudeTarget.val    = {testSNR};
testAmplitudeTarget.help   = {'Choose whether target an amplitude level'
    'Or a signal-to-noise ratio (SNR).'}';


testAmplitude2         = cfg_entry;
testAmplitude2.name    = '2nd Volterra response amplitude';
testAmplitude2.tag     = 'testAmplitude2';
testAmplitude2.strtype = 'r';
testAmplitude2.num     = [1 Inf];
testAmplitude2.def     = @(val)nirs_get_defaults(...
    'readOnsets.addTestStimuli.testAmplitude2', val{:});
testAmplitude2.help    = {'Enter desired response amplitude as a percentage '
    'of the 1st Volterra kernel amplitude.'
    'Only used when 2nd Volterra is turned on.'}';

% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
voltAddStim         = cfg_menu;
voltAddStim.tag     = 'voltAddStim';
voltAddStim.name    = 'Model Interactions (Volterra)';
voltAddStim.help    = {
    'Generalized convolution of inputs (U) with basis set (bf).'
    ''
    'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
    'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
    }';
voltAddStim.labels = {
    'Do not model Interactions'
    'Model Interactions'
    }';
voltAddStim.values = {1, 2};
voltAddStim.val = {1};

testGamma         = cfg_menu;
testGamma.tag     = 'testGamma';
testGamma.name    = 'Use Gamma function';
testGamma.help    = {
    'Use gamma function or canonical HRF for test'
    'Note that canonical HRF is leads to lower estimability due to the undershoot'}';
testGamma.labels = {
    'Gamma Function'
    'Canonical HRF'
    }';
testGamma.values = {1, 2};
testGamma.val = {1};

testFilterX         = cfg_menu;
testFilterX.tag     = 'testFilterX';
testFilterX.name    = 'Filter design matrix';
testFilterX.help    = {
    'Filter design matrix prior to calculating power of protocol'}';
testFilterX.labels = {
    'Filter On'
    'Filter Off'
    }';
testFilterX.values = {1, 0};
testFilterX.val = {1};

testFilterData         = cfg_menu;
testFilterData.tag     = 'testFilterData';
testFilterData.name    = 'Filter data';
testFilterData.help    = {
    'Filter data channels prior to calculating power of baseline'}';
testFilterData.labels = {
    'Filter On'
    'Filter Off'
    }';
testFilterData.values = {1, 0};
testFilterData.val = {1};

testStdvsPower        = cfg_menu;
testStdvsPower.tag     = 'testStdvsPower';
testStdvsPower.name    = 'Normalization choice';
testStdvsPower.help    = {
    'Normalize added protocol based on channel standard deviation or power'
    'Either option should give the same result...'}';
testStdvsPower.labels = {
    'Standard deviation'
    '(square root of) Power'
    }';
testStdvsPower.values = {1, 0};
testStdvsPower.val = {1};


testBfNorm        = cfg_menu;
testBfNorm.tag     = 'testBfNorm';
testBfNorm.name    = 'Further Normalization choice';
testBfNorm.help    = {
    'Normalize further by rescaling both 1st and 2nd Volterra by the max value of the 1st Volterra'}';
testBfNorm.labels = {
    'Normalize'
    'Do not Normalize'
    }';
testBfNorm.values = {1, 0};
testBfNorm.val = {0};

testHPFButterOn         = cfg_menu;
testHPFButterOn.tag     = 'testHPFButterOn';
testHPFButterOn.name    = 'Butter HPF';
testHPFButterOn.help    = {
    'Preferred option: ON, at 0.004'}';
testHPFButterOn.labels = {
    'Filter On'
    'Filter Off'
    }';
testHPFButterOn.values = {1, 0};
testHPFButterOn.val = {1};

testHPFbutterCutoff         = cfg_entry;
testHPFbutterCutoff.name    = 'HPF Butter cutoff';
testHPFbutterCutoff.tag     = 'testHPFbutterCutoff';
testHPFbutterCutoff.strtype = 'r';
testHPFbutterCutoff.num     = [1 1];
testHPFbutterCutoff.val = {0.004};
testHPFbutterCutoff.help    = {'HPF cutoff in Hz. (0.004 Hz recommended)'}';

testHPFbutterOrder         = cfg_entry;
testHPFbutterOrder.name    = 'HPF Butter order';
testHPFbutterOrder.tag     = 'testHPFbutterOrder';
testHPFbutterOrder.strtype = 'r';
testHPFbutterOrder.num     = [1 1];
testHPFbutterOrder.val = {5};
testHPFbutterOrder.help    = {'HPF order (3 recommended)'}';

testLPFGaussianOn         = cfg_menu;
testLPFGaussianOn.tag     = 'testLPFGaussianOn';
testLPFGaussianOn.name    = 'Gaussian LPF';
testLPFGaussianOn.help    = {
    'Preferred option: ON, at 1.5 s'}';
testLPFGaussianOn.labels = {
    'Filter On'
    'Filter Off'
    }';
testLPFGaussianOn.values = {1, 0};
testLPFGaussianOn.val = {1};

testLPFGaussianFWHM         = cfg_entry;
testLPFGaussianFWHM.name    = 'LPF Gaussian FWHM';
testLPFGaussianFWHM.tag     = 'testLPFGaussianFWHM';
testLPFGaussianFWHM.strtype = 'r';
testLPFGaussianFWHM.num     = [1 1];
testLPFGaussianFWHM.val = {1.5};
testLPFGaussianFWHM.help    = {'LPF Gaussian FWHM in seconds ( 1.5 s recommended)'}';

testWaveletMDLOn         = cfg_menu;
testWaveletMDLOn.tag     = 'testWaveletMDLOn';
testWaveletMDLOn.name    = 'Wavelet MDL HPF';
testWaveletMDLOn.help    = {
    'Preferred option: Off'}';
testWaveletMDLOn.labels = {
    'Filter On'
    'Filter Off'
    }';
testWaveletMDLOn.values = {1, 0};
testWaveletMDLOn.val = {0};


% Executable Branch
addTestStimuli      = cfg_exbranch;
addTestStimuli.name = 'Add Stimuli with HRFs for testing';
addTestStimuli.tag  = 'addTestStimuli';
addTestStimuli.val  = {NIRSmat redo1 NIRSmatCopyChoice testStimulusName testStimuliNumber ...
    testSessionNumber testWavelength testAmplitudeTarget ...
    voltAddStim testAmplitude2 keepAllChannels testChannels testDupChannels testPType ...
    testGamma testFilterX testFilterData testStdvsPower testBfNorm ...
    testHPFButterOn testHPFbutterCutoff testHPFbutterOrder ...
    testLPFGaussianOn testLPFGaussianFWHM testWaveletMDLOn};
addTestStimuli.prog = @nirs_run_addTestStimuli;
addTestStimuli.vout = @nirs_cfg_vout_addTestStimuli;
addTestStimuli.help = {'Module to add stimuli convoluted with HRFs for'
    'simulation and testing purposes'}';

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_addTestStimuli(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});