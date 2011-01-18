function out = nirs_run_boxy(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: generate NIRS.mat info and convert Boxy files to .nir binary files, 

%Input files
%prjname: .prj file for project
%fileIn: all raw Boxy data
%fileOut: .txt, .mat, .vhdr, .vmkr

%mat_Mtg: 2-dim array with the x, y, z coordinates as the columns and, in
%the rows, each source first, then the detectors, then the electrodes

jobP = job.config_path;
N_subj = size(job.subj,2);
%Big loop over all subjects
for Idx_subj=1:N_subj 
    %Reinitialize NIRS matrix for each subject
    tNIRS = job.configuration; %temporary NIRS structure, will not be saved
    %System Configuration: NIRS.Cf
    NIRS = [];
    NIRS.Cf.dev.n = 'ISS Imagent';
    NIRS.Cf.dev.wl = job.configuration.Lambda;
    %NIRS.Cf.dev.gn = 1; %gain?
    NIRS.Cf.dev.fs = job.configuration.freq/job.configuration.resample;
    %if job.configuration.resample > 1
    
    
    %end
    %Extract info from chosen files
    %Number of BOXY data files for this subject
    tNIRS.N_dfiles = size(job.subj(1,Idx_subj).fnames,1);
    %NIRS.N_dfiles = size(job.fnames,1);
    
    %path and expname
    temp_fname = job.subj(1,Idx_subj).fnames(1,:);
    tNIRS.BOXYfiles = temp_fname;
    %temp_fname = job.fnames(1,:);
    [subj_BOXY_path expname] = fileparts(temp_fname{1,1});
    %find subject root directory
    temp_idx = strfind(subj_BOXY_path,'\');
    % = 
    
    %Data information
    NIRS.Dt.s.p = fullfile(subj_BOXY_path(1:temp_idx(end)-1), filesep); %was NIRS.subj_path
    
    tNIRS.prjname = [NIRS.subj_path jobP.prj_path '\' expname '.prj']; 
    %path for raw data
    NIRS.pathboxy = subj_BOXY_path;
    %path for output of NIRS.mat and of .nir data files
    NIRS.pathoutput = [NIRS.subj_path jobP.output_path '\'];
    %create this output directory if it doesn't already exist
    if ~isdir(NIRS.pathoutput), mkdir(NIRS.pathoutput); end    
    %Path for T1
    NIRS.pathoutput_T1 = [NIRS.subj_path jobP.T1_path '\'];
    %create T1 output directory if needed
    if ~isdir(NIRS.pathoutput_T1), mkdir(NIRS.pathoutput_T1); end 
    
    %load project montage
    try
        LS = load(NIRS.prjname,'-mat');    
        %assign shorter names to project structures
        mtg = LS.SaveStruct.m_Helmet.Mtg_Data;
        holes = LS.SaveStruct.m_Helmet.v_Holes;
        %coordinates of fiducials, converted to centimeters
        mat_Fid = 100*LS.SaveStruct.m_DigSubjFids.matFiducials; 
    catch
        disp('Problem loading HSJ montage file. Check .prj file is in correct place. Aborting.');
        return
    end
    clear LS %get rid of this big project structure

    %find real world x,y,z coordinates of sources, detectors & electrodes; 
    %sources followed by detectors are in mat_Mtg
    [mat_Mtg mat_Ele] = nirs_boxy_get_montage_coordinates(mtg.v_pSrc,...
                            mtg.v_pDet,mtg.v_pEle,holes,mtg.v_HolesMtg);

    %Add to NIRS file describing montage, to be used by boxy_convert
    %Complex code very specific to BOXY ISS Imagent
    %convert distmin and distmax to meters, as they will be converted back 
    %to centimeters later on.
    NIRS = nirs_boxy_associate_sources_detectors(NIRS,NIRS.distmax/100,NIRS.distmin/100,mat_Mtg,...
            NIRS.nb_Det,mat_Ele,mtg.v_pEle,mtg.v_HolesEle,NIRS.MaxSources);

    %no need anymore for mtg, mat_Mtg, mat_Ele, holes
    clear mtg mat_Mtg mat_Ele holes

    %if electrodes are not specified in the montage, load a standard reference,
    %scaled to the fiducials... not straightforward

    %Continue building the channels into structure NIRS 
    %Approximate positions of channels by 10-20 or 10-10 electrode system
    NIRS = nirs_boxy_approx10_20or10_10(NIRS,NIRS.use10_10system);

    %Extend channel names and lengths to channels of both wavelengths
    NIRS = nirs_boxy_extend_channel_names_lengths(NIRS);
    
    %add useful info to NIRS.mat: frequency, sampling interval
    %NIRS.SamplingFrequency = NIRS.freq; %in Hertz
    NIRS.SamplingInterval =floor(1000000/NIRS.freq); %in microseconds

    %Big Loop over each of the BOXY data files 
    for Idx_file = 1:NIRS.N_dfiles
        %read, resample and write data -- also some reordering of channels: DS.ml
        %     bN.fileIn = [];
    %     for int_bloc = 1:length(bN.mat_Bloc(:,1)) %Make block of input file names 
    %         bN.fileIn = [bN.fileIn; sprintf('%s%s%s%03.0f', ...
    %             bN.pathboxy,bN.expname,'.', bN.mat_Bloc(int_bloc))];
    %     end
        %Output file name
        temp_fName = job.subj(1,Idx_subj).fnames(Idx_file,:);
        %temp_fName = job.fnames(Idx_file,:);
        fName = temp_fName{1,1};
        [path1 expname1 ext1] = fileparts(fName); %#ok<ASGLU>
        ext1 = ext1(2:end); %skip the initial dot "."
        shortfileOutRoot = [expname1 '_' ext1];
        fileOutRoot = [NIRS.pathoutput shortfileOutRoot];
        %expname_cond = %sprintf('%s%s',bN.expname,char(bN.cell_Cond));
        %fileOutRoot = [bN.pathoutput bN.expname_cond];
        fileOut=[fileOutRoot '.nir'];
        fileOutRoot_vhdr = [fileOutRoot '.vhdr'];
        fileOutRoot_vmrk = [fileOutRoot '.vmrk'];
        fileOut_nir = [shortfileOutRoot '.nir'];
        fileOut_vmrk = [shortfileOutRoot '.vmrk'];
         
        NIRS = nirs_boxy_convert(NIRS,fName,fileOut,Idx_file); 
    
        NIRS.NIRfile{Idx_file} = fileOut;
        
        %Header file
        %if job.write_Analyzer_header
        nirs_boxy_write_header(fileOutRoot_vhdr,... %Output file
                fileOut_nir,... %DataFile
                fileOut_vmrk,... %MarkerFile,...
                'nirs_convert_boxy',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                NIRS.ChannelLabels,... %names given as a column of cells 
                NIRS.SamplingInterval); %SamplingInterval in microseconds
        %end
        %Marker file
        %if job.write_Analyzer_markers
        temp_markers{1,1} = NIRS.Markers{Idx_file,1};
        nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);     
    %end
    
    end %for Files_idx -- Big loop over BOXY files for this subject  

    NIRS.FidPos = mat_Fid;
    NIRS.n_Fid = size(mat_Fid,1);
    NIRS.n_Src = size(NIRS.SrcPos,1);
    NIRS.n_Det = size(NIRS.DetPos,1);
    NIRS.n_Opt = NIRS.n_Src+ NIRS.n_Det;
    NIRS.n_Chn = size(NIRS.ChnPos,1);    
    
    %for NIRS_SPM and coregistration in MNI coordinates between T1 MRI and
    %NIRS: write in T1 folder text file of fiducial, source and detector
    %positions, putting number of fiducials, and sum of number of sources and
    %of detectors in the name of the text file as they are required inputs in NIRS_SPM
    fid = fopen([NIRS.pathoutput_T1 'FidSrcDet_' int2str(NIRS.n_Fid) '_' int2str(NIRS.n_Opt) '.txt'],'wt');
    for Idx=1:NIRS.n_Fid
        fprintf(fid,'%s%s%s%s%s\n',num2str(mat_Fid(Idx,1)),' ',...
            num2str(mat_Fid(Idx,2)),' ',num2str(mat_Fid(Idx,3)));
    end
    for Idx=1:NIRS.n_Src
        fprintf(fid,'%s%s%s%s%s\n',num2str(NIRS.SrcPos(Idx,1)),' ',...
            num2str(NIRS.SrcPos(Idx,2)),' ',num2str(NIRS.SrcPos(Idx,3)));
    end
    for Idx=1:NIRS.n_Det
        fprintf(fid,'%s%s%s%s%s\n',num2str(NIRS.DetPos(Idx,1)),' ',...
            num2str(NIRS.DetPos(Idx,2)),' ',num2str(NIRS.DetPos(Idx,3)));
    end
    fclose(fid);

    %also for NIRS_SPM and coregistration, create text file of channel pairings
    %between sources and detectors
    fid = fopen([NIRS.pathoutput_T1 'Source Detector Pairs.txt'], 'wt');
    fprintf(fid,'%s\n','HSJ');
    fprintf(fid,'%s%s%s%s%s%s\n',int2str(NIRS.n_Src),'x',int2str(NIRS.n_Det),'_',int2str(NIRS.n_Chn),'ch');
    fprintf(fid,'%s\n\n','1set');
    for Idx=1:NIRS.n_Chn    %convention: add total number of sources n_Src to detector number
        fprintf(fid,'%s%s%s\n',int2str(NIRS.ml(Idx,1)),'  ',int2str(NIRS.ml(Idx,2)+NIRS.n_Src));
    end
    fclose(fid);

    %write out NIRS in .mat file
    save([NIRS.subj_path 'NIRS'],'NIRS'); 
    %make NIRS.mat available as virtual output for the first subject only
    if Idx_subj == 1
        out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
    end
end %end for  Big loop over subjects