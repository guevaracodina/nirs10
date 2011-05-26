function out = nirs_run_runMC(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________

load(job.NIRSmat{1,1});

%%%%%%%% meme code que generate sensitivity matrix
try
    f = job.MC_runCUDAchoice.MCX1.MCXconfigFiles;
catch
    disp('a coder')
end
cs_dir =  fileparts(f{1,:});
cs_ldir = cs_dir(max(strfind(cs_dir,'\'))+9:end);

ics =1;
while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
    ics =ics+1;
end
cs = NIRS.Cs.mcs{ics};
%%%%%%%%

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
            copyfile([spm('Dir') '\toolbox\nirs10\mc_exe\mcx.exe'],[dir1 '\mcx.exe']);% le premier etait mcx_dir
        end
%         res = system(['mcx.exe -t 2048 -T 64 -n 1e5 -f ' file2 ' -s ' file1 ' -r 10 -g 1 -U 1 -d 1 -a 0 -b 0']);
        res = system(['mcx.exe -t 4800 -T 480 -n ' cs.par.nphotons ' -f ' file2 ' -s ' file1 ' -r 1 -g 1 -U 1 -d 1 -a 0 -b 0']);
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
out = 1;%fullfile(NIRS.subj_path,'NIRS.mat');