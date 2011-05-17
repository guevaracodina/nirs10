function out = nirs_run_runMC(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________

%Simulate MCX
if isfield(job.MC_runCUDAchoice,'MCX1')
    t = job.MC_runCUDAchoice.MCX1.MCXconfigFiles;
    %run MCX
    for Idx=1:size(t,1)
        if size(t,1) == 1
            nfile = t;
        else
            nfile=t(Idx,:);
        end
        [dir1 file1] = fileparts(nfile{:});
        file2 = [file1 '.inp'];
        cd(dir1);
        if Idx == 1
            copyfile('D:\Users\Clément\MonteCarlo\mcx.exe',[dir1 '\mcx.exe']);% le premier etait mcx_dir
        end
        tic
%         plus r est petit, plus la simu est rapide et plus ca consomme de 
% ressources
%         res = system(['mcx.exe -t 2048 -T 64 -n 1e5 -f ' file2 ' -s ' file1 ' -r 10 -g 1 -U 1 -d 1 -a 0 -b 0']);
        res = system(['mcx.exe -t 4800 -T 480 -n 1e6 -f ' file2 ' -s ' file1 ' -r 400 -g 1 -U 1 -d 1 -a 0 -b 0']);
        toc
    end
    delete([dir1 '\mcx.exe']);
    
elseif isfield(job.MC_runCUDAchoice,'tMCimg1')
    t = job.MC_runCUDAchoice.tMCimg1.tMCimg_configFiles;
    jobCD.dir{1} = fileparts(t{1,1});
    cfg_run_cd(jobCD);
    for it = 1:size(t,1)
        t_it = t{it,1};
        str_run = ['D:\Users\Clément\MonteCarlo\tMCimg.exe ' t_it(end-15:end-4)];
        res = system(str_run);
    end
end
% on peut charger la matrice NIRS pour y inclure les nouveaux directory
% !!!!
out = fullfile(NIRS.subj_path,'NIRS.mat');