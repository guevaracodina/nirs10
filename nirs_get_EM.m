function EM = nirs_get_EM(job)
%EM parameters
EM.spm_integrator = job.EM_parameters.spm_integrator;
EM.Niterations = job.EM_parameters.Niterations;
EM.dFcriterion = job.EM_parameters.dFcriterion;
EM.LogAscentRate = job.EM_parameters.LogAscentRate;
EM.MaxLogAscentRate = job.EM_parameters.MaxLogAscentRate;
EM.Mstep_iterations = job.EM_parameters.Mstep_iterations;
