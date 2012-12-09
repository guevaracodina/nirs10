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
%BUT WHAT IS THE 4th column of mat_Mtg?
%mat_Chn: same as mat_Mtg, but for the channels

jobP = job.config_path;
N_subj = size(job.subj,2);
outNIRSmat = {};
%Big loop over all subjects
for Idx_subj=1:N_subj
    tmpNIRS = [job.subj2.Apath{1} jobP.output_path filesep 'NIRS.mat'];
    if spm_existfile(tmpNIRS) && ~job.force_redo
        load(tmpNIRS);
    else
        NIRS.flags = [];
    end
    if ~isfield(NIRS.flags,'NIRS_OK') || job.force_redo
        %Reinitialize NIRS matrix for each subject
        NIRS = job.cf1;
        if NIRS.resample > 1
            NIRS.freq = NIRS.freq/NIRS.resample;
        end
        %Extract info from chosen files
        %Number of BOXY data files for this subject
        NIRS.N_dfiles = size(job.subj(1,Idx_subj).fnames,1);
        %NIRS.N_dfiles = size(job.fnames,1);
        
        %path and expname
        temp_fname = job.subj(1,Idx_subj).fnames(1,:);
        NIRS.BOXYfiles = temp_fname;
        %temp_fname = job.fnames(1,:);
        [subj_BOXY_path expname] = fileparts(temp_fname{1,1});
        %find subject root directory
        temp_idx = strfind(subj_BOXY_path,filesep);
        NIRS.subj_path = [subj_BOXY_path(1:temp_idx(end)-1) filesep];
        NIRS.prjname = [NIRS.subj_path jobP.prj_path filesep expname '.prj'];
        %path for raw data
        NIRS.pathboxy = subj_BOXY_path;
        %path for output of NIRS.mat and of .nir data files
        NIRS.pathoutput = [NIRS.subj_path jobP.output_path filesep];
        %create this output directory if it doesn't already exist
        if ~isdir(NIRS.pathoutput), mkdir(NIRS.pathoutput); end
        %Path for T1
        NIRS.pathoutput_T1 = [NIRS.subj_path jobP.T1_path filesep];
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
            %Output file name
            temp_fName = job.subj(1,Idx_subj).fnames(Idx_file,:);
            fName = temp_fName{1,1};
            [path1 expname1 ext1] = fileparts(fName); %#ok<ASGLU>
            ext1 = ext1(2:end); %skip the initial dot "."
            shortfileOutRoot = [expname1 '_' ext1];
            fileOutRoot = [NIRS.pathoutput shortfileOutRoot];
            fileOut=[fileOutRoot '.nir'];
            fileOutRoot_vhdr = [fileOutRoot '.vhdr'];
            fileOutRoot_vmrk = [fileOutRoot '.vmrk'];
            fileOut_nir = [shortfileOutRoot '.nir'];
            fileOut_vmrk = [shortfileOutRoot '.vmrk'];
            
            NIRS = nirs_boxy_convert(NIRS,fName,fileOut,Idx_file);
            
            NIRS.NIRfile{Idx_file} = fileOut;
            %Header file
            nirs_boxy_write_header(fileOutRoot_vhdr,... %Output file
                fileOut_nir,... %DataFile
                fileOut_vmrk,... %MarkerFile,...
                'nirs_convert_boxy',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                NIRS.ChannelLabels,... %names given as a column of cells
                NIRS.SamplingInterval); %SamplingInterval in microseconds
            %Marker file
            temp_markers{1,1} = NIRS.Markers{Idx_file,1};
            nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);
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
        
        
        oldNIRS = NIRS;
        
        %Fill out standardized NIRS structure
        NIRS = [];
        NIRS.flags.NIRS_OK = 1;
        NIRS.Cf.dev.n = 'ISS Imagent';
        NIRS.Cf.dev.wl = job.cf1.Lambda;
        %NIRS.Cf.dev.gn = 1; %gain?
        NIRS.Cf.dev.fs = job.cf1.freq/job.cf1.resample;
        NIRS.Dt.s.p = oldNIRS.subj_path;
        %Try adding the anatomical image
        if ~isempty(job.subj(1,Idx_subj).anatT1)
            NIRS.Dt.ana.T1 = job.subj(1,Idx_subj).anatT1{1};
        end
        NIRS.Dt.s.age = job.subj(1,Idx_subj).age1;
        %Data
        NIRS.Dt.fir.pp(1).p = oldNIRS.NIRfile'; %?
        NIRS.Dt.fir.pp(1).pre = 'readBOXY';
        NIRS.Dt.fir.pp(1).job = job;
        %Information on Helmet
        NIRS.Cf.H.n = 'Sainte-Justine';
        
        %Information on fiducials
        NIRS.Cf.H.F.r.o.mm.p = mat_Fid';
        %sources
        %NIRS.Cf.H.S.n
        NIRS.Cf.H.S.N = oldNIRS.n_Src;
        NIRS.Cf.H.S.r.o.mm.p = oldNIRS.SrcPos';
        %detectors
        %NIRS.Cf.H.D.n
        NIRS.Cf.H.D.N = oldNIRS.n_Det;
        NIRS.Cf.H.D.r.o.mm.p = oldNIRS.DetPos';
        %channels
        NIRS.Cf.H.C.n = [oldNIRS.ChnNames' oldNIRS.ChnNames'];
        NIRS.Cf.H.C.N = 2*oldNIRS.n_Chn;
        NIRS.Cf.H.C.id = [1:size(oldNIRS.ml,1); oldNIRS.ml(:,1:2)'];
        NIRS.Cf.H.C.wl = oldNIRS.ml(:,4)';
        NIRS.Cf.H.C.gp = [oldNIRS.ChnDist' oldNIRS.ChnDist'];
        %try generating onsets - only saving location for now
        try
            NIRS.Dt.fir.rons = job.subj(1,Idx_subj).raw_onset_files;
        catch
            disp('Problem with onset files');
        end
        %write out NIRS in .mat file
        save([oldNIRS.pathoutput 'NIRS'],'NIRS');
        save([oldNIRS.pathoutput 'oldNIRS'],'oldNIRS');
        outNIRSmat = [outNIRSmat; fullfile(oldNIRS.pathoutput,'NIRS.mat')];
    else
        outNIRSmat = [outNIRSmat; tmpNIRS];
    end
end %end for  Big loop over subjects
out.NIRSmat = outNIRSmat;