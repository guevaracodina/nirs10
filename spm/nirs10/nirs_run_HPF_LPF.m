function out = nirs_run_HPF_LPF(job)
%filename prefix
prefix = 'f'; %for "filter" 
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
                
        
        fs = NIRS.Cf.dev.fs;
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        Cgp = NIRS.Cf.H.C.gp;
        Cwl = NIRS.Cf.H.C.wl;
        NC = NIRS.Cf.H.C.N;
        wl = NIRS.Cf.dev.wl;
        
        %LPF
        try
            fwhm = job.nirs_lpf2.lpf_gauss2.fwhm1;
            df = job.nirs_lpf2.lpf_gauss2.downsamplingFactor;
            
            %Boolean for applying low-pass filtering
            LPFon = 1;
            %Boolean for applying downsampling
            if df > 1
                DSon = 1;
            else
                DSon = 0;
            end
        catch
            LPFon = 0;
            DSon = 0;
        end
        
        %HPF
        KH = [];
        try
            HPF = ['DCT, ' int2str(job.nirs_hpf.hpf_dct.hpf_dct_cutoff)];
            HPFon = 1;
            KH.HParam.type = 'DCT';
        catch
            try 
                job.nirs_hpf.hpf_wavelet;
                HPF = ['wavelet,' int2str(job.nirs_hpf.hpf_wavelet.hpf_wavelet_iter)];
                HPFon = 1;
                KH.HParam.type = 'Wavelet-MDL';
            catch
                try 
                    job.nirs_hpf.hpf_none;
                    HPF = 'none';
                    HPFon = 0;
                catch
                    disp('Unrecognized high pass filter');
                    HPFon = 0;
                end
            end
        end
     
        try
            %loop over data files
            for f=1:size(rDtp,1)
                [dir1 fil1 ext1] = fileparts(rDtp{f});          
                d = fopen_NIR(rDtp{f},NC);                 
                %bring in markers - to write out to the next NIRS.Dt.fir.pp
                %link and importantly for filtering over segments
                try 
                    bpi = NIRS.Dt.fir.pp(lst).bpi{f,1}; %bad point indices
                    bpd = NIRS.Dt.fir.pp(lst).bpd{f,1}; %bad point durations
                    markers_available = 1;
                catch
                    bpi = [];
                    bpd = [];
                    markers_available = 0;
                end
                
                %Define filter
                if LPFon
                    K = [];
                    K.LParam.type = 'Gaussian';
                    K.LParam.FWHM = fwhm;
                    K.RT = 1/fs;
                    K.row = 1:size(d,2);
                    K.HParam.type = 'none';
                    K = spm_filter_HPF_LPF_WMDL(K); %get filter structure
                end

                if LPFon
                    [NIRS d  bpi bpd] = NIRS_SPM_filter(NIRS,K,d,DSon,...
                        df,bpi,bpd,markers_available);                   
                end
                
                if HPFon
                    K = KH;
                    %taken from precoloring code of NIRS_SPM
                    switch K.HParam.type
                        case 'Wavelet-MDL'
                            K.X = X; %Need design matrix
                        case 'DCT'
                            SPM.xX.xKXs = spm_sp('Set', ...
                                spm_filter_HPF_LPF_WMDL(K, X)); % KX
                    end

                    %filtering of the data
                    KY = spm_filter_HPF_LPF_WMDL(K, d');
                    d = KY';
                end
              
                outfile = fullfile(dir1,[prefix fil1 ext1]);
                fwrite_NIR(outfile,d);
                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'HPF_LPF';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
                NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
                NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
            end
        catch
            disp('Could not high pass and low pass filter');
        end
        
        save(job.NIRSmat{Idx,1},'NIRS');
    catch
        disp(['Conversion of optical intensities to hemoglobin ',...
            'concentrations failed for subject ' int2str(Idx)]);
    end
    
end
out.NIRSmat = job.NIRSmat;