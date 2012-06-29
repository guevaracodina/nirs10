function LPF = nirs_get_LPF(job)
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
else
    LPF.lpf_gauss_On = 0;
end