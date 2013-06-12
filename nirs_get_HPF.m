function HPF = nirs_get_HPF(job)
if isfield(job,'hpf_filter')
    %NEW, SCKS
    if isfield(job.hpf_filter,'SPM_cosine_filter')
        HPF.hpf_filter = 1;
        HPF.cosine_freq = job.hpf_filter.SPM_cosine_filter.cosine_freq;
        if isfield(job.hpf_filter.SPM_cosine_filter.band_pass_filter,'band_pass_filter_on')
            HPF.band_pass_filter = 1;
            HPF.cosine_freq_band = job.hpf_filter.SPM_cosine_filter.band_pass_filter.band_pass_filter_on.cosine_freq_band;
        else
            HPF.band_pass_filter = 0;
        end
    else
        if isfield(job.hpf_filter,'hpf_butter_On')
            HPF.hpf_filter = 2;
            HPF.hpf_butter_freq = job.hpf_filter.hpf_butter_On.hpf_butter_freq;
            HPF.hpf_butter_order = job.hpf_filter.hpf_butter_On.hpf_butter_order;
        else
            if isfield(job.hpf_filter,'hpf_filter_Off')
                HPF.hpf_filter = 0;
            end
        end
    end
else %OLD, HDM
    if isfield(job.hpf_butter,'hpf_butter_On')
        HPF.hpf_butter_On = 1;
        HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
        HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
    else
        HPF.hpf_butter_On = 0;
    end
end