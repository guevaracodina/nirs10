function out = nirs_run_generate_vhdr_vmrk(job)
group_consecutive_markers = 1; %boolean
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
           
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        SamplingInterval =floor(1000000/fs); %in microseconds
        ChannelLabels = NIRS.Cf.H.C.n;
        
        for f=1:size(rDtp,1)
            %find markers
            switch NIRS.Dt.fir.pp(lst).pre 
                case 'remove_chn_stdev'
                    mtype = 'Bstd';
                case 'mark_negative'
                    mtype = 'Bneg';
                case 'mark_movement'
                    mtype = 'Bmvt';
                otherwise
                    mtype = 'Bmvt';
            end
            %marker positions, in data points
            bpi = NIRS.Dt.fir.pp(lst).bpi{f,1};
            bpd = NIRS.Dt.fir.pp(lst).bpd{f,1};
           
%             if group_consecutive_markers
%                 [bpi bpd] = find_marker_durations(bpi);
%             end
            %reset and generate cell array
            Markers = {};
            Markers{length(bpi)+1,1} = {};
            %create markers structure
            Markers{1}.Type='New Segment';
            Markers{1}.Description='';
            Markers{1}.Position=1;
            Markers{1}.Size=1;
            Markers{1}.ChNumber=0; %all channels
            %loop over markers
            for i=1:length(bpi)
                Markers{i+1}.Type='Comment';
                Markers{i+1}.Description=mtype;
                Markers{i+1}.Position=bpi(i);
                Markers{i+1}.Size=bpd(i);
                Markers{i+1}.ChNumber=0; %all channels
            end
            
            [dir1 fil1 ext1] = fileparts(rDtp{f,1});
                    
            shortfileOutRoot = fil1;
            fileOutRoot = fullfile(dir1,shortfileOutRoot);
            %fileOut=[fileOutRoot '.nir'];
            fileOutRoot_vhdr = [fileOutRoot '.vhdr'];
            fileOutRoot_vmrk = [fileOutRoot '.vmrk'];
            fileOut_nir = [shortfileOutRoot '.nir'];
            fileOut_vmrk = [shortfileOutRoot '.vmrk'];
        
            %Header file
            nirs_boxy_write_header(fileOutRoot_vhdr,... %Output file
                    fileOut_nir,... %DataFile
                    fileOut_vmrk,... %MarkerFile,...
                    'generate_vhdr_vmrk',... %Function that created the header
                    '',... %Channel Resolution
                    '',... %Channel Units
                    ChannelLabels,... %names given as a column of cells 
                    SamplingInterval); %SamplingInterval in microseconds
            %Marker file
            nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                    fileOut_nir,... %DataFile
                    Markers);      
            
            
        end                   
    catch
        disp(['Could not generate .vhdr and .vmrk Analyzer ',...
            'files for subject' int2str(Idx)]);
    end
        
    %save(job.NIRSmat{Idx,1},'NIRS');    
end
out.NIRSmat = job.NIRSmat;