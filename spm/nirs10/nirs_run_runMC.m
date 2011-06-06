function out = nirs_run_runMC(job)
% si on a la meme origine pour les voxels et les mm, on peut avoir 1 1 1 et
% ca marche. Le code colle la source sur un voxel touchant la tete et il
% lance la simulation
%
% 
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________
%%%%%%%% same code as sensitivity matrix

try
    MCtestOneChannel = job.MCtestOneChannel;
catch
    MCtestOneChannel = 0;
end
try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir0,dummy,dummy2] = fileparts(job.NIRSmat{Idx,1});
        if NewDirCopyNIRS
            dir2 = [dir0 filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
        else
            dir2 = dir0;
        end
        
        f = job.MC_runCUDAchoice.MCX1.MCXconfigFiles;
        if ~isempty(f{1})
            cs_dir =  fileparts(f{1,:});
        else
            cs_dir = dir2;
        end
        %%%%%attention errreur !!!!!
        if strcmp(cs_dir(max(strfind(cs_dir,'\'))+1:max(strfind(cs_dir,'\'))+3),'roi')
            seps = strfind(cs_dir,'\');
            cs_ldir = cs_dir(max(seps(1:end-1))+(length('MC')+1):max(seps)-1);
        else
            cs_ldir = cs_dir(max(strfind(cs_dir,'\'))+(length('MC')+1):end);
        end
        ics =1;
        while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
            ics =ics+1;
        end
        cs = NIRS.Cs.mcs{ics};
        %%%%%%%%

        %Simulate MCX
        if isfield(job.MC_runCUDAchoice,'MCX1')
            
            [t,dummy] = spm_select('FPList',cs_dir,'.inp');
            %run MCX
            if MCtestOneChannel             
                t = t(1,:); %This gives one detector, but also need at least one source, with the same wavelength, for the one channel
            end
            J = job.MC_runCUDAchoice.MCX1;
            for k1=1:size(t,1)
                [dir1,file1,dummy] = fileparts(t(k1,:));
                file2 = [file1 '.inp'];
                cd(dir1);
        
                if k1 == 1
                    copyfile([spm('Dir') '\toolbox\nirs10\mc_exe\mcx.exe'],[dir1 '\mcx.exe']);% le premier etait mcx_dir
                end
                %         res = system(['mcx.exe -t 2048 -T 64 -n 1e5 -f ' file2 ' -s ' file1 ' -r 10 -g 1 -U 1 -d 1 -a 0 -b 0']);
                str_run1 = ['mcx.exe -t ' int2str(J.MCX_t) ' -T ' int2str(J.MCX_T) ' -n ' int2str(cs.par.nphotons)];
                if J.MCX_l
                    str_log = ' -l'; %'MCX_logfile';
                else
                    str_log = '';
                end
                str_run2 = [' -r ' int2str(J.MCX_r) ' -g ' int2str(J.MCX_g) ' -U 1 -S 1 -d 1 -a 0 -b 0'];
                res = system([str_run1 ' -f ' file2 ' -s ' file1 str_run2 str_log]);
            end
            delete([dir1 '\mcx.exe']);
            
        elseif isfield(job.MC_runCUDAchoice,'tMCimg1')
            t = job.MC_runCUDAchoice.tMCimg1.tMCimg_configFiles;
            jobCD.dir{1} = fileparts(t{1,1});
            cfg_run_cd(jobCD);
            for it = 1:size(t,1)
                t_it = t{it,1};
                str_run = [[spm('Dir') 'toolbox\nirs10\mc_exe\tMCimg.exe'] t_it(end-15:end-4)];
                res = system(str_run);
            end
        end
        NIRS.Cs.mcs{ics}.MCX_t = J.MCX_t;
        NIRS.Cs.mcs{ics}.MCX_T = J.MCX_T;
        NIRS.Cs.mcs{ics}.MCX_r = J.MCX_r;
        NIRS.Cs.mcs{ics}.MCX_g = J.MCX_g;
        NIRS.Cs.mcs{ics}.MCX_l = J.MCX_l;
        if NewDirCopyNIRS
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{Idx,1} = newNIRSlocation;
        else
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(['Could not run MonteCarlo simulation for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;