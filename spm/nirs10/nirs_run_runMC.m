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

% le nombre de repetitions a une influence sur le DPF mais les valeurs
% restent de meme ordre de grandeur
% si on garde la meme seed pour plusieurs simulations (tous les autres parametres etant les memes), on obtient la meme chose


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
        
        try
            f = job.MC_runCUDAchoice.MCX1.MCXconfigFiles;
        catch
            f = job.MC_runCUDAchoice.tMCimg1.tMCimg_configFiles;
        end
        if ~isempty(f{1})
            cs_dir =  fileparts(f{1,:});
        else
            cs_dir = dir2;
        end
        
        [dummy,cs_ldir] = fileparts(cs_dir);
        ics =1;
        while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
            ics =ics+1;
        end
        cs = NIRS.Cs.mcs{ics};
        %%%%%%%%
        
        switch cs.alg
            case 1 %MCX
                [t,dummy] = spm_select('FPList',cs_dir,'.inp');
                
                if MCtestOneChannel
                    t = t(1,:); %This gives one detector, but also need at least one source, with the same wavelength, for the one channel
                end
                J = job.MC_runCUDAchoice.MCX1;
                countS =0;
                countD =0;
                for k1=1:size(t,1)
                    [dir1,file1,dummy] = fileparts(t(k1,:));
                    file2 = [file1 '.inp'];
                    cd(dir1);
                    
                    % if his/mch needed
                    if strcmp(file1(1),'S')
                        codeexe = 'mcx_det';
                        if countS ==0
                         copyfile([spm('Dir') '\toolbox\nirs10\mc_exe\' codeexe '.exe'],[dir1 filesep codeexe '.exe']);
                         countS = countS+1;
                        end
                    else
                        codeexe = 'mcx';
                        if countD ==0
                         copyfile([spm('Dir') '\toolbox\nirs10\mc_exe\' codeexe '.exe'],[dir1 filesep codeexe '.exe']);
                         countD = countD+1;
                        end
                    end
                    %else
                    %codeexe = 'mcx.exe';
                    %end
                    
                    %                 str_run1 = ['mcx.exe -t ' int2str(J.MCX_t) ' -T ' int2str(J.MCX_T) ' -n ' int2str(cs.par.nphotons)];
% % % % % % % % % % %                     str_run1 = [codeexe ' -A 2 -n ' int2str(cs.par.nphotons)];
% % % % % % % % % % %                     if J.MCX_l
% % % % % % % % % % %                         str_log = ' -l'; %'MCX_logfile';
% % % % % % % % % % %                     else
% % % % % % % % % % %                         str_log = '';
% % % % % % % % % % %                     end
% % % % % % % % % % %                     str_run2 = [' -r ' int2str(J.MCX_r) ' -g ' int2str(J.MCX_g) ' -U 1 -S 1 -d 1 -a 0 -b 0 -B 1 -z 1'];
% % % % % % % % % % %                     res = system([str_run1 ' -f ' file2 ' -s ' file1 str_run2 str_log]);

% % % % % % % % % % %                     str_run2 = [' -r ' int2str(J.MCX_r) ' -g ' int2str(J.MCX_g) ' -U 1 -S 1 -d 1 -a 0 -b 0 -B 1 -z 1'];
% -z 1 EST ESSENTIEL !!!!!!!!!!!!!!!!!!!
%                     res = system([codeexe ' -A -n ' int2str(cs.par.nphotons) ' -f ' file2 ' -s ' file1 ' -r ' int2str(J.MCX_r) ' -g 1 -b 0 -d 1 -z 1 -l']);
                    res = system([codeexe ' -A -n ' int2str(cs.par.nphotons) ' -f ' file2 ' -s ' file1 ' -r ' int2str(J.MCX_r) ' -U 1 -S 1 -g 1 -b 0 -d 1 -z 1 -l']);
                end
                
                if countD>0, delete([dir1 filesep 'mcx.exe']);end
                if countS>0, delete([dir1 filesep 'mcx_det.exe']);end
                
                NIRS.Cs.mcs{ics}.MCX_t = J.MCX_t;
                NIRS.Cs.mcs{ics}.MCX_T = J.MCX_T;
                NIRS.Cs.mcs{ics}.MCX_r = J.MCX_r;
                NIRS.Cs.mcs{ics}.MCX_g = J.MCX_g;
                NIRS.Cs.mcs{ics}.MCX_l = J.MCX_l;
                
            case 2 %tMCimg
                [t,dummy] = spm_select('FPList',cs_dir,'.cfg');
                
                if MCtestOneChannel
                    t = t(1,:); %This gives one detector, but also need at least one source, with the same wavelength, for the one channel
                end
                J = job.MC_runCUDAchoice.tMCimg1;
                jobCD.dir{1} = fileparts(t(1,:));
                cfg_run_cd(jobCD);
                for k1=1:size(t,1)
                    [dir1,file1,dummy] = fileparts(t(k1,:));
                    str_run = [[spm('Dir') '\toolbox\nirs10\mc_exe\tMCimg.exe'] ' ' file1];
                    res = system(str_run);
                end
        end

        if NewDirCopyNIRS
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{Idx,1} = newNIRSlocation;
        else
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not run MonteCarlo simulation for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;