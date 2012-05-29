function out = nirs_run_ODtoHbOHbR(job)

% Filename prefix for data file containing dOD as "d" (data)
prefix = 'h'; % for "hemoglobin" concentrations
% Loop over subjects

%**************************************************************************
%Introduce nirs data jump filling process
%Ke Peng
%**************************************************************************

if isfield(job.nirs_filling_jumps, 'nirs_fill_jumps_on')
    fill_jump_on = 1;
    num_standard_deviation = jobs.nirs_filling_jumps.nirs_fill_jumps_on.num_standard_deviation;
    num_points = jobs.nirs_filling_jumps.nirs_fill_jumps_on.num_points;
    size_gap = jobs.nirs_filling_jumps.nirs_fill_jumps_on.size_gap;
    
    if isfield(job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable, 'HPF_enable_on')
        HPF_enable_on = 1;
        hpf_butter_order = job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable.HPF_enable_on.hpf_butter_order;
        hpf_butter_freq = job.nirs_filling_jumps.nirs_fill_jumps_on.HPF_enable.HPF_enable_on.hpf_butter_freq;
    else
        HPF_enable_on = 0;
    end
else
    fill_jump_on = 0;
end

%**************************************************************************




for Idx=1:size(job.NIRSmat,1)
    % Load NIRS.mat information
    try
        clear PVF PVF2 DPF EPF
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'concOK') || job.force_redo)
            
            %bl_m = job.Normalize_OD;% method to calculate baseline
            
            age = NIRS.Dt.s.age;
            % Perform computation on output of the last step of preprocessing
            % that has been performed
            lst = length(NIRS.Dt.fir.pp);
            rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
            Cgp = NIRS.Cf.H.C.gp; % source-detector distances
            Cwl = NIRS.Cf.H.C.wl; % channels wavelength indexes (e.g. wl #1 or #2)
            NC = NIRS.Cf.H.C.N; % number of channels
            wl = NIRS.Cf.dev.wl; % device wavenlengths
            
            % Partial volume correction factor
            if isfield(job.PVF,'PVFval') && ~isempty(job.PVF.PVFval)
                PVF = job.PVF.PVFval; % 1 x nLambda
                PVF = repmat(PVF',[1 NC]); % nLambda x NC
            elseif isfield(job.PVF,'PVFsim') && ~isempty(job.PVF.PVFsim)
                simPVF = load(job.PVF.PVFsim{:});
                PVF = simPVF.PVF; % 1 x NC
                nLambda = length(wl);
                PVF2 = zeros(nLambda,NC);
                for iwl = 1:nLambda
                    for Ci = 1:NC
                        PVF2(Cwl(1,Ci),Ci) = PVF(1,Ci);
                    end
                end
                PVF = PVF2; % nLambda x NC
            else % for compatibility with older version
                PVF = job.PVF; % 1 x nLambda
                PVF = repmat(PVF',[1 NC]); % nLambda x NC
            end
            
            try
                %exs(:,1): HbO for each wavelength; exs(:,2): HbR for each wavelength
                %Alexis's choice of extinction coefficients corresponds to case 1 in
                %Homer, which appears to be their preferred choice too
                [exs,NIRS.Dt.pro.extcoeff_ref] = GetExtinctions(wl,1);
            catch
                disp('Problem loading extinction coefficients in nirs_run_ODtoHbOHbR.');
            end
            
            % Differential pathlength factor
            if isfield(job,'DPF') && isfield(job.DPF,'DPFval') && ~isempty(job.DPF.DPFval)
                DPF = job.DPF.DPFval; % 1 x nLambda
                DPF = repmat(DPF',[1 NC]); % nLambda x NC
            elseif isfield(job,'DPF') && isfield(job.DPF,'DPFsim') && ~isempty(job.DPF.DPFsim)
                load(job.DPF.DPFsim{:}); % nTissues x NC
                DPF = sum(PDPF_Vi,1); % 1 x NC
                nLambda = length(wl);
                DPF2 = zeros(nLambda,NC);
                for iwl = 1:nLambda
                    for Ci = 1:NC
                        DPF2(Cwl(1,Ci),Ci) = DPF(1,Ci);
                    end
                end
                DPF = DPF2; % nLambda x NC
            else % Use literature value
                % From Alexis Machado
                % Differential path length factor --
                % From 1996 Duncan et al. Measurement of cranial optical path length...
                table = [690 ; 830 ; 744 ; 807];
                a = [5.38 ; 4.67 ; 5.11 ; 4.99];
                b = [0.049 ; 0.062 ; 0.106 ; 0.067];
                c = [0.877 ; 0.819 ; 0.723 ; 0.814];
                for iwl =1:size(wl,2)
                    [raw,col] = find(wl(iwl)==table);
                    DPF(iwl,1) = a(raw,col)+b(raw,col)*age^c(raw,col);
                end
                % nLambda x 1
                DPF = repmat(DPF,[1 NC]); % nLambda x NC
            end
            
            EPF = zeros(1,NC); %EPF = L * PPF = L * DPF/PVF at each wavelength
            for Ci = 1:NC
                %EPF(1,Ci) = Cgp(1,Ci)*DPF(Cwl(1,Ci),1)./PVF(1,Cwl(1,Ci));
                EPF(1,Ci) = Cgp(1,Ci)*DPF(Cwl(1,Ci),Ci)./PVF(Cwl(1,Ci),Ci);
                %             EPF(1,Ci) = Cgp(Ci,1)*DPF(Cwl(1,Ci),1)./PVF(1,Cwl(1,Ci)); %PP??? why not ./???
            end
            
            inv_exs = pinv(exs(:,1:2)); % size 2 x #wl
            inv_exs2 = kron(inv_exs,eye(NC/size(wl,2))); % size #pairs*2 x #pairs*#wl
            % (number of channels NC = number of pairs x number of wavelengths)
            
            %try
            %loop over data files
            for f=1:size(rDtp,1)
                %[dir1 fil1 ext1] = fileparts(rDtp{f});
                %if ~strcmp(fil1(1:3),prefix)       %???
                d = fopen_NIR(rDtp{f},NC);
                %bring in markers - to write out to the next NIRS.Dt.fir.pp
                %link and importantly for filtering over segments
                
                
%**************************************************************************
%Filling jumps in nirs data
%Ke Peng
%**************************************************************************
                
                if fill_jump_on
                    
                    OP.Sb = num_standard_deviation;
                    OP.Nr = num_points;
                    OP.Mp = size_gap;
                    
                    if HPF_enable_on
                        OP.ubf = 1;
                        OP.fs = NIRS.Cf.dev.fs;
                        OP.bf = hpf_butter_freq;
                        OP.bo = hpf_butter_order;
                    else
                        OP.ubf = 0;
                    end
                    
                    d = nirs_remove_jumps(d,OP);
                end


%**************************************************************************
                try
                    bpi = NIRS.Dt.fir.pp(lst).bpi{f,1}; %bad point indices
                    bpd = NIRS.Dt.fir.pp(lst).bpd{f,1}; %bad point durations
                    si = NIRS.Dt.fir.pp(lst).si{f,1};
                    ei = NIRS.Dt.fir.pp(lst).ei{f,1};
                    if ~isempty(bpi)
                        markers_available = 1;
                    else
                        markers_available = 0;
                    end
                catch
                    bpi = [];
                    bpd = [];
                    markers_available = 0;
                end
                
                %Multiply by 1e6 to get micromolar units
                %negative sign so that an increase in chromophore
                %concentration corresponds to a decrease in intensity due
                %to light absorption
                if markers_available
                    for i1 =1:length(si)
                        d(:,si(i1):ei(i1)) = -1e6 * real(log(d(:,si(i1):ei(i1))));
                    end
                else
                    d = -1e6 * real(log(d));
                end
                
                %Effective path length
                EPF2 = ones(size(d,2),1)*EPF; %PP used to be EPF2 = EPF *ones(1,size(d,2));
                d = d ./ EPF2'; %PPused to be d = d ./ EPF2;
                %MBLL - c consists now of HbO and HbR, even if we had more
                %than two wavelengths to begin
                c = inv_exs2 * d;
                % size (#pairs*2 x nt) = (#pairs*2 x #pairs*#wl) x (#pairs*#wl x nt)
                % 2 here is for 2 Hb types, HbO and HbR
                
                % c has structure:
                %  [ conc_HbO_pair1(time1) conc_HbO_pair1(time2) ...;
                %    conc_HbO_pair2(time1) conc_HbO_pair2(time2) ...;
                %    ...
                %    conc_HbR_pair1(time1) conc_HbR_pair1(time2) ...;
                %    ... ];
                
                
                % %             CB: corrige dans les lignes qui suivent
                % %             if isnan(c(:))
                % %                 disp(['Some elements are NaN for file ' int2str(f)]);
                % %                 NIRS.WARNING = 'Some NaN elements in data - GLMs will fail!';
                % %                 %CB: puisqu'on corrige est ce qu'on garde le WARNING ????
                % %             end
                
                % Loop over channels and Hb types
                for iC=1:size(c,1)
                    try
                        tNaN =(1:size(c,2)).*isnan(c(iC,:));
                        tNaNo = tNaN(tNaN~=0);
                        % Replace by previous time point
                        c(iC, tNaNo) = c(iC, tNaNo-1);
                        if iC==size(c,1)
                            NIRS.WARNING = 'Some NaN elements were corrected';
                        end
                    catch
                        NIRS.WARNING = 'Some NaN elements in data - GLMs will fail!';
                    end
                end
                [dir1,fil1,ext1] = fileparts(rDtp{f});
           
                outfile = fullfile(dir1,[prefix fil1 ext1]);

                fwrite_NIR(outfile,c);
                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'ODtoHbOHbR';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
                try
                    NIRS.Dt.fir.pp(lst+1).bpi{f,1} = bpi; %bad point indices
                    NIRS.Dt.fir.pp(lst+1).bpd{f,1} = bpd; %bad point durations
                    NIRS.Dt.fir.pp(lst+1).si{f,1} = si;
                    NIRS.Dt.fir.pp(lst+1).ei{f,1} = ei;
                catch
                end
            end
            NIRS.flags.concOK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Conversion of optical intensities to hemoglobin ',...
            'concentrations failed for subject ' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;
end

