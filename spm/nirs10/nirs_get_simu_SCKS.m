function S = nirs_get_simu_SCKS(job)
if isfield(job.simuOn,'simuYes')
    S.simuOn = 1;
    S.simuS     = job.simuOn.simuYes.simuS; %Stimuli types to include
    S.simuIt    = job.simuOn.simuYes.simuIt; %Number of random iterations
    S.simuA     = job.simuOn.simuYes.simuA; %Signal amplitude, as % of BOLD signal
    S.simuP     = job.simuOn.simuYes.simuP; %Parameters to vary
    if isfield(job.simuOn.simuYes.simuParamDistr,'distr_uniform')
        S.simuPriorDistr = 1;
        S.simuPrior = job.simuOn.simuYes.simuParamDistr.distr_uniform.simuPrior; %Priors to use
        S.simuR     = job.simuOn.simuYes.simuParamDistr.distr_uniform.simuR; %Range to sample
    end
    if isfield(job.simuOn.simuYes.simuParamDistr,'distr_bimodal')
        S.simuPriorDistr = 2;
        S.simuPrior1 = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuMean1; %Priors to use
        S.simuR1     = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuR1; %Range to sample
        S.simuDiffPrior = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuMean21; %Priors to use
        S.simuR2     = job.simuOn.simuYes.simuParamDistr.distr_bimodal.simuR2; %Range to sample
    end

    S.simuUpsample = job.simuOn.simuYes.simuUpsample; %Upsampling factor on data
    try % for back-compatibility
        S.simuInterp = job.simuOn.simuYes.simuInterp; %Interpolation on simulated data (after decimating)
    catch
        S.simuInterp = 1;
    end
    simuNoise1 = job.simuOn.simuYes.simuNoise; %Yes to include background noise based on restscans
    if isfield(simuNoise1,'noiseYes')
        S.simuNoise = 1;
        switch modal
            case 1
                restscans = job.simuOn.simuYes.simuNoise.noiseYes.restscans; %Rest scans to add signal to
        end
    else
        S.simuNoise = 0;
    end
else
    S.simuOn = 0; 
end