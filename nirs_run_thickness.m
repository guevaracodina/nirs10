function out = nirs_run_thickness(job)

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{Idx,1});
    
    % Get "cs" (current simulation) info
    if ~isempty(job.dir_in{1,1})
        cs_dir =  fileparts(job.dir_in{1,1});
        [dummy cs_ldir] = fileparts(cs_dir);
        ics =1;
        while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
            ics =ics+1;
        end
    else
        [cs_dir dummy dummy1] = fileparts(job.NIRSmat{iSubj,1});
        ics = length(NIRS.Cs.n); % Use last simulation run
    end
    cs = NIRS.Cs.mcs{ics};
    if cs.alg==1, Oe='.mch'; elseif cs.alg==2, Oe='.his';end
    
    mu_subj ={};
    mu_subj.litt ={};
    jobW.mu_subj = mu_subj;
    jobW.pve_cfg=0;
    jobW.alg=1;
    jobW.NSkpt = 1;
    jobW.NDkpt = 1;
    jobW.NSinit=NIRS.Cf.H.S.N;
    jobW.Pvoid={};
    jobW.n_b8i =  cs.n_b8i;
    jobW.dir = [cs.dir 'thickness'];
    if ~exist(jobW.dir,'dir')
        mkdir(jobW.dir);
    end
    
    % Get selected P
    esp = [0 strfind(job.SP,' ') length(job.SP)+1];
    
    % Mise a jour des positions c1
    %%%%% calcul des positions
    V = spm_vol(cs.segR);
    inv_mat = spm_imatrix(V.mat);
    
    for i=1:length(esp)-1
        SPn = job.SP(esp(i)+1:esp(i+1)-1);
        switch SPn(1)
            case 'S'
                jobW.Pkpt = str2num(SPn(2:end));
            case 'D'
                jobW.Pkpt = str2num(SPn(2:end))+NIRS.Cf.H.S.N;
            otherwise
                ind = (1:NIRS.Cf.H.Q.N).*strcmp(NIRS.Cf.H.Q,SPn);
                jobW.Pkpt = ind+NIRS.Cf.H.S.N+NIRS.Cf.H.D.N;
        end
        % Mise a jour des positions c1
        Pp_c1_rmm = NIRS.Cs.temp.Pp_roi_c1_rmm(:,jobW.Pkpt);
        % Positions : Transform MNI mm -> MNI isotropic voxels
        Pp_c1_rmv = V.mat\[Pp_c1_rmm;1];
        Pp_c1_rmiv= abs(inv_mat(7:9)').*(Pp_c1_rmv(1:3)/cs.par.voxelSize);
        jobW.P={};
        jobW.P.p = [cs.P.p(:,jobW.Pkpt) Pp_c1_rmiv];
        jobW.P.wd= [cs.P.wd(:,jobW.Pkpt) -cs.P.wd(:,jobW.Pkpt)];
        Sr = cs.par.radiis;
        Dr = cs.par.radiid;
        jobW.P.r = [Sr' Dr' zeros(1,size(cs.Pkpt,1)-2)];
        
        jobW.par.nphotons = cs.par.nphotons;
        jobW.par.seed = cs.par.seed;
        jobW.par.numTimeGates = cs.par.numTimeGates;
        jobW.par.deltaT = cs.par.deltaT;
        jobW.par.voxelSize = cs.par.voxelSize;
        jobW.ROIlimits = cs.ROIlimits;
        jobW.nummed = cs.nummed;
        jobW.wl_dev = NIRS.Cf.dev.wl;
        
        mem = jobW.Pkpt;
        clear jobW.Pkpt;
        jobW.Pkpt = [mem mem];
        
        % Write the files
        nirs_configMC_writeCFGfiles2(jobW);
    end
    copyfile(cs.b8i,fullfile(jobW.dir,cs.n_b8i));
    
    % MCsim
    [t,dummy] = spm_select('FPList',jobW.dir,'.inp');
    countS =1;
    for k1=1:size(t,1)
        [dir1,file1,dummy] = fileparts(t(k1,:));
        file2 = [file1 '.inp'];
        cd(dir1);
        
        if strcmp(file1(1),'S')
            %%% %%% MC simulation %%% %%%
            codeexe = 'mcx_det';
            if isempty(spm_select('FPList',dir1,'mcx_det.exe'))
                copyfile([spm('Dir') '\toolbox\nirs10\mc_exe\' codeexe '.exe'],[dir1 filesep codeexe '.exe']);
            end
            system([codeexe ' -A -n 1e6 -f ' file2 ' -s ' ...
                file1 ' -r 100 -g 1 -b 0 -d 1 -z 0']);
            
            ms=loadmc2([file1 '.mc2'],[V.dim 1],'float');
            nirs_create_vol(fullfile(jobW.dir,[file2 '_Green.nii']),...
                    V.dim, [16,0], V.pinfo, V.mat, log(ms));
            %%% %%% thickness computation %%% %%%
            file3 = [file1 '.mch'];
            [history, header]=loadmch(file3);
            [dummy, tk_file, dummy] = fileparts(file3);
            History{countS,1} = tk_file; % name of file
            History{countS,3} = header; % header info
            History{countS,2} = history; % for each detected photon:
            % det # / number of scattering events / pathlength
            % through each layer (nlayers+2 columns)
            countS = countS+1;
        end
    end
    if k1==size(t,1), delete([dir1 filesep 'mcx_det.exe']);end
    save('History')
end
out.NIRSmat = job.NIRSmat;
%         %segmented T1
%         V = spm_vol(NIRS.Dt.ana.T1seg);
%         Y = spm_read_vols(V);
%
%         max_size = V.dim';
%         min_size = [1;1;1];
%
%         %fitted positions on skin and on cortex
%         Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
%         Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
%         Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
% % SP sampling points

%         SPp_rmm = Pp_rmm(SP);
%         SPfp_rmm = Pfp_rmm(SP);
%         SPp_c1_rmm = Pp_c1_rmm(SP);
%         % - directions : pointing outside of the head
%         SPd_rmm = -(SPp_rmm - SPp_c1_rmm)/norm(SPp_rmm - SPp_c1_rmm);
%
%         T=zeros(5,size(SP));%thicknesses below all selected points
%         for Pi = 1:size(SP)
%             bPi = zeros(3,5);%boundary points below Pi
%             bPi(:,1) = SPp_c1_rmm(:,Pi);%c1
%
%             Kim = bPi(:,1);%current position
%             Kiv = V.mat\[bPi(:,1);1];
%             Kiv = round(Kiv(1:3));
%             %initialization : projection on WM
%             while Y(Kiv(1),Kiv(2),Kiv(3))==1
%                 Kim = Kim + prec*SPd_rmm(:,Pi);
%
%                 Kiv = V.mat\[Kim;1];
%                 Kiv = round(Kiv(1:3));
%                 Kiv = min(max(Kiv,min_size),max_size);
%             end
%             bPi(:,2) = Kim - prec*SPd_rmm(:,Pi);
%
%             %find other boundary points
%             Kim = bPi(:,1);
%             Kiv = V.mat\[bPi(:,1);1];
%             Kiv = round(Kiv(1:3));
%
%             for j=[1 3 4 5] % juste sur les images sans masque
%                 while Y(Ki_v(1),Ki_v(2),Ki_v(3))==j
%                     Ki_m = Ki_m - prec*SPd_rmm(:,Pi);
%
%                     Ki_v = V.mat\[Ki_m;1];
%                     Ki_v = round(Ki_v(1:3));
%                     Ki_v = min(max(Ki_v,min_size),max_size);
%                 end
%                 Ti(:,j) = Ki_m + prec*SPd_rmm(:,Pi);
%
%                 switch j % computing thicknesses
%                     case 1
%                         T(1,Pi) = norm(Ti(:,2) - Ti(:,1));
%                     case 2
%                         T(2,Pi) = 10;
%                     case 3
%                         T(3,Pi) = norm(Ti(:,3) - Ti(:,1));
%                     otherwise
%                         T(j,Pi) = norm(Ti(:,j-1) - Ti(:,j));
%                 end
%             end
%         end