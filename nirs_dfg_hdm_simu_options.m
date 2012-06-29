function simuOn = nirs_dfg_hdm_simu_options(BOLDon)
% % % Simulations % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simuA         = cfg_entry;
simuA.tag     = 'simuA';
simuA.name    = '% Signal amplitude';
simuA.help    = {'Enter signal amplitude, as a percentage of the signal (e.g. enter 1 for a 1% amplitude)'};
simuA.strtype = 'e';
simuA.num     = [1 1];
simuA.val     = {100};

simuS         = cfg_entry;
simuS.tag     = 'simuS';
simuS.name    = 'which_condition to simulate';
simuS.help    = {'Enter array of stimuli types to simulated.'
    'Enter 0 to include all stimuli types.'}';
simuS.strtype = 'e';
simuS.num     = [1 Inf];
simuS.val     = {0};

simuP         = cfg_entry;
simuP.tag     = 'simuP';
simuP.name    = 'Parameters to randomize';
simuP.help    = {'Enter array of parameters to be sampled.'
    'Enter 0 to randomize all parameters.'}';
simuP.strtype = 'e';
simuP.num     = [1 Inf];
simuP.val     = {1};

simuIt         = cfg_entry;
simuIt.tag     = 'simuIt';
simuIt.name    = 'Number of simulations';
simuIt.help    = {'Enter number of simulations'};
simuIt.strtype = 'e';
simuIt.num     = [1 1];
simuIt.val     = {1};

noiseNo             = cfg_branch;
noiseNo.tag         = 'noiseNo';
noiseNo.name        = 'No noise in simulated data';
noiseNo.val         = {};
noiseNo.help        = {'Simulated data will be noiseless (0 baseline).'};


% % % Experimental noise % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restscans         = cfg_files;
restscans.tag     = 'restscans';
restscans.name    = 'Resting state scans';
restscans.help    = {'Select resting state scans on which the simulated data will be added.'
    'The code was developed to work specifically on the 4D.nii scans in folder'
    '11-ep2d_bold_ax_4x4x4_TR_1010 of Michele subject 28'
    'Assumptions are that these scans are  compatible with the SPM and xSPM structures'
    '(same number number of scans in particular)'}';
restscans.filter = 'image';
restscans.ufilter = '.*';
restscans.val     = {''};
restscans.num     = [0 Inf];

restscans_BOLD         = cfg_files;
restscans_BOLD.tag     = 'restscans_BOLD';
restscans_BOLD.name    = 'Resting state scans for BOLD';
restscans_BOLD.help    = {'ONLY for simulating BOLD+ASL'}';
restscans_BOLD.filter = 'image';
restscans_BOLD.ufilter = '.*';
restscans_BOLD.val     = {''};
restscans_BOLD.num     = [0 Inf];

restscans_ASL         = cfg_files;
restscans_ASL.tag     = 'restscans_ASL';
restscans_ASL.name    = 'Resting state scans for ASL';
restscans_ASL.help    = {'ONLY for simulating BOLD+ASL'}';
restscans_ASL.filter = 'image';
restscans_ASL.ufilter = '.*';
restscans_ASL.val     = {''};
restscans_ASL.num     = [0 Inf];

noiseYes             = cfg_branch;
noiseYes.tag         = 'noiseYes';
noiseYes.name        = 'Add noise to simulated data';
if BOLDon
    noiseYes.val         = {restscans restscans_BOLD restscans_ASL};
else
    noiseYes.val         =  {};
end
noiseYes.help        = {['Noise will be added to simulated data before inversion.'...
    ' This noise is added in the form of resting state data measured experimentally.' ...
    ' Specify those experimental scans below.']};

simuNoise           = cfg_choice;
simuNoise.name      = 'Include baseline noise';
simuNoise.tag       = 'simuNoise';
simuNoise.values    = {noiseYes noiseNo};
simuNoise.val       = {noiseYes};
simuNoise.help      = {'Include noise background scans; if No, the simulated data will be noiseless, i.e. on 0 background.'}';

simuUpsample         = cfg_entry;
simuUpsample.tag     = 'simuUpsample';
simuUpsample.name    = 'Data upsampling factor';
simuUpsample.help    = {'Enter an upsampling factor (max = 16)'};
simuUpsample.strtype = 'e';
simuUpsample.num     = [1 1];
simuUpsample.val     = {1};

simuInterp         = cfg_entry;
simuInterp.tag     = 'simuInterp';
simuInterp.name    = 'Data interpolation factor';
simuInterp.help    = {'Enter an interpolation factor for simulated data. Data will first be '...
    ' simulated at 16 time bins per TR, then decimated to TR, or between 1-1/16 of TR '...
    ' if an upsampling factor is specified. Then, noise will be added if specified. ' ...
    ' Only then, data will be interpolated by this factor. Large interpolation factors will increase '...
    ' computation time significantly.' };
simuInterp.strtype = 'e';
simuInterp.num     = [1 1];
simuInterp.val     = {1};

simuR         = cfg_entry;
simuR.tag     = 'simuR';
simuR.name    = 'Parameter range';
simuR.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled randomly uniformly between (100-range) and (100+range) % of the prior value.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR.strtype = 'e';
simuR.num     = [1 Inf];
simuR.val     = {25};

simuPrior         = cfg_entry;
simuPrior.tag     = 'simuPrior';
simuPrior.name    = 'Prior values of parameters';
simuPrior.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
    'If nothing is entered, the default prior values will be used.'}';
simuPrior.strtype = 'e';
simuPrior.num     = [0 Inf];
simuPrior.val     = {''};

distr_uniform           = cfg_branch;
distr_uniform.tag       = 'distr_uniform';
distr_uniform.name      = 'Uniform';
distr_uniform.val       = {simuPrior simuR};
distr_uniform.help      = {'Simulated priors drawn pseudo-randomly from a uniform distribution.'};

simuMean1         = cfg_entry;
simuMean1.tag     = 'simuMean1';
simuMean1.name    = 'Mean for Gaussian #1';
simuMean1.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
    'If nothing is entered, default prior values will be used.'};
simuMean1.strtype = 'e';
simuMean1.num     = [0 Inf];
simuMean1.val     = {''};

% SimuMean2         = cfg_entry;
% SimuMean2.tag     = 'SimuMean2';
% SimuMean2.name    = 'Mean for Gaussian #2';
% SimuMean2.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
%     'If nothing is entered, default prior values will be used.'};
% SimuMean2.strtype = 'e';
% SimuMean2.num     = [0 Inf];
% SimuMean2.val     = {''};

simuMean21         = cfg_entry;
simuMean21.tag     = 'simuMean21';
simuMean21.name    = '% Difference between Gaussian means';
simuMean21.help    = {['Enter array, for all parameters, of % values that ' ...
    ' specify the difference between the means of the 2 Gaussians distributions.' ...
    ' For example, if you specify "25" for one parameter, the first half of parameters are drawn ' ...
    'from a Gaussian distribution with mean M specified above, and the second half ' ...
    'from a Gaussian with 1.25*M.' ...
    'If nothing is entered, all parameters will be drawn from 1 Gaussian distribution.']};
simuMean21.strtype = 'e';
simuMean21.num     = [0 Inf];
simuMean21.val     = {''};

simuR1         = cfg_entry;
simuR1.tag     = 'simuR1';
simuR1.name    = 'Parameter range #1';
simuR1.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled from a gaussian distribution with sigma = this range.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR1.strtype = 'e';
simuR1.num     = [1 Inf];
simuR1.val     = {25};

simuR2         = cfg_entry;
simuR2.tag     = 'simuR2';
simuR2.name    = 'Parameter range #2';
simuR2.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled from a gaussian distribution with sigma = this range.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR2.strtype = 'e';
simuR2.num     = [1 Inf];
simuR2.val     = {25};

distr_bimodal           = cfg_branch;
distr_bimodal.tag       = 'distr_bimodal';
distr_bimodal.name      = '2 Gaussians';
distr_bimodal.val       = {simuMean1 simuMean21 simuR1 simuR2};
distr_bimodal.help      = {['Half of simulated priors drawn pseudo-randomly '...
    'each from a gaussian distribution.']};

simuParamDistr         = cfg_choice;
simuParamDistr.tag     = 'simuParamDistr';
simuParamDistr.name    = 'Parameter distribution';
simuParamDistr.values  = {distr_uniform distr_bimodal};
simuParamDistr.val     = {distr_uniform};
simuParamDistr.help  = {['Type of distribution from which are pseudo-randomly' ...
    'drawn the parameters used to simulate data (priors).']};

simuYes         = cfg_branch;
simuYes.tag     = 'simuYes';
simuYes.name    = 'HDM on simulated data';
simuYes.val     = {simuIt simuA simuS simuP simuParamDistr simuUpsample simuInterp simuNoise};
simuYes.help    = {'Perform HDM on simulated data'}';

simuNo         = cfg_branch;
simuNo.tag     = 'simuNo';
simuNo.name    = 'HDM on real data';
simuNo.val     = {};
simuNo.help    = {'No simulations; perform HDM on real data'}';

simuOn         = cfg_choice;
simuOn.tag     = 'simuOn';
simuOn.name    = 'Perform simulations';
simuOn.values  = {simuNo simuYes};
simuOn.val     = {simuNo};
simuOn.help    = {'Perform simulations'}';
