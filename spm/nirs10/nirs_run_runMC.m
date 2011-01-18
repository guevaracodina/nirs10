function out = nirs_run_runMC(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage:

%Simulate MCX
%change working directory as 'system' command cannot take a space in the
%filename path
start_dir = pwd;
algo = job.MC_runCUDAchoice;

if isfield(job.MC_runCUDAchoice,'MCX1')
    t = job.MCXconfigFiles;
    %run MCX
    for Idx=1:size(t,1)
        if size(t,1) == 1
            nfile = t;
        else
            nfile=t(Idx,:);
        end
        [dir1 file1] = fileparts(nfile);
        file2 = [file1 '.inp'];
        cd(dir1);
        if Idx == 1
            copyfile([mcx_dir '\mcx.exe'],[dir1 '\mcx.exe']);
        end
        str_run = ['mcx.exe -t 3584 -T 128 -g 10 -m 100000 -f ' file2 ' -s ' file1 ' -r 10 -a 0 -b 0 -l'];
        res=system(str_run);
    end
    delete([dir1 '\mcx.exe']);
    cd(start_dir);
    
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
out = 1;%.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');