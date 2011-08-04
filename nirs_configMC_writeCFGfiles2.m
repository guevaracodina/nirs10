function nirs_configMC_writeCFGfiles2(job)
% Writes configurqtion files for either MCX or TMCimg
% FORMAT nirs_configMC_writeCFGfiles(P,wl)
% algo       - algorithm used to run the MC simulation
% n_b8i      - name of the bin volume
% mc_dir     - directory of MC simulation
% n          - name of the bin volume
% ROIlimits  - limits of the ROI
% parameters - tissues properties and simulation parameters
% NS         - number of sources
% ND         - number of detectors
% Pvoid      - void points
%_______________________________________________________________________
%
% Do the job that needs to be done
% Write each line in the config files...
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% mise a jour Clement 07/2011


outOP = GetOpt_ppts('wl');
opt_ppts = outOP{1};
opt_ppts_perturb = outOP{2};

pve_cfg = job.pve_cfg;
algo = job.alg;
NSinit = job.NSinit;
NS = job.NSkpt;
ND = job.NDkpt;
Pvoid = job.Pvoid;
Pkpt = job.Pkpt;
Pp_rmiv = job.P.p;
Pwd_rmiv = job.P.wd;
r = job.P.r;
mc_dir = job.dir;

n_b8i = job.n_b8i;

% Create a .cfg or .inp file for each optode and each wavelength
for iwl = 1:size(job.wl_dev,2)
    for iP = 1:NS+ND
        if isempty(Pvoid) || ~Pvoid(1,iP) %skip optodes with no data &&&'nS',iS,wl(iwl)
            if iP<=NS
                if Pkpt(iP)<10, PNo = ['S_No' num2str(0) num2str(Pkpt(iP))]; else PNo = ['S_No' num2str(Pkpt(iP))]; end
            else
                if Pkpt(iP)-NSinit<10, PNo = ['D_No' num2str(0) num2str(Pkpt(iP)-NSinit)]; else PNo = ['D_No' num2str(Pkpt(iP)-NSinit)]; end
            end
            
            if algo==2 % tMCimg
                n_cfg = fullfile(mc_dir,[PNo '_' num2str(job.wl_dev(iwl)) 'nm.cfg']);
                fid = fopen(n_cfg,'w');
                disp('Generating files for sensitivity profiles computation')
                
                % Title
                fprintf(fid, ';; Fichier de configuration\n');
                % ID subject
                fprintf(fid, [';; Sujet # : ' n_b8i '\n\n']);
                % ROI limits relative to whole MRI
                fprintf(fid, ';; ROI limits relative to whole MRI\n');
                fprintf(fid, [';; xmin : ' num2str(job.ROIlimits(1,2)) '\n']);
                fprintf(fid, [';; xmax : ' num2str(job.ROIlimits(2,2)) '\n']);
                fprintf(fid, [';; ymin : ' num2str(job.ROIlimits(1,1)) '\n']);
                fprintf(fid, [';; ymax : ' num2str(job.ROIlimits(2,1)) '\n']);
                fprintf(fid, [';; zmin : ' num2str(job.ROIlimits(1,3)) '\n']);
                fprintf(fid, [';; zmax : ' num2str(job.ROIlimits(2,3)) '\n\n']);
                % Number of photons
                fprintf(fid, ';; Number of photons\n');
                fprintf(fid, 'nphotons = %12.0f\n\n', job.par.nphotons);
                % Seed random number generator
                fprintf(fid, ';; Seed random number generator\n');
                fprintf(fid, 'seed    =  %6.0f\n\n', job.par.seed);
                % Modulation frequency
                fprintf(fid, ';; Modulation frequency\n');
                fprintf(fid, 'freq    = %g \n\n', job.par.modulationFreq);
                % Optode (source)
                fprintf(fid, ';; Source position.  All lengths in mm.\n');
                % Here watch out! The tMCimg code is written in C/C++ thus matrix
                % encoding is row-major; in Matlab it is on the contrary column-major
                % thus for Matlab-generated files one must invert x and y
                % dimensions
                fprintf(fid, 'source { pos = [%4.0f %4.0f %4.0f]\n', Pp_rmiv(2,iP), Pp_rmiv(1,iP), Pp_rmiv(3,iP));%% ici pas besoin d'inverser les colonnes x et y
                fprintf(fid, '	 dir = [%5.4f %5.4f %5.4f] \n', Pwd_rmiv(2,iP), Pwd_rmiv(1,iP), Pwd_rmiv(3,iP));%% ici pas besoin d'inverser les colonnes x et y
                fprintf(fid, ' 	 rad = %1.2f       } \n\n', r(1,iP));
                
                % Time gates
                fprintf(fid, ';; Times are in seconds.\n');
                fprintf(fid, 'start_time = 0.00e-9    \n'); % Always start at t=0...
                fprintf(fid, 'gate_width = %3.2e     \n', job.par.deltaT);
                fprintf(fid, 'ngate      = %2.0f         \n\n', job.par.numTimeGates);
                
                % Name of segmentation file
                fprintf(fid, ';; Segmentation file to load\n');
                %[~,file,ext] = fileparts(job.n);
                fprintf(fid, 'segfile  %s \n\n', n_b8i);%[file ext]);
                
                % Voxel size
                fprintf(fid, ';; Voxel size.  Square voxels only, for now (dz = dy = dx)\n');
                fprintf(fid, 'dx     %1.0f \n\n', job.par.voxelSize);
                
                % System dimensions - see note above about inverting x and y dimensions
                fprintf(fid, ';; Size of system (in voxels) \n');
                fprintf(fid, 'nxvox =  %3.0f       \n', job.ROIlimits(2,2));
                fprintf(fid, 'nyvox =  %3.0f       \n', job.ROIlimits(2,1));
                fprintf(fid, 'nzvox =  %3.0f       \n\n', job.ROIlimits(2,3));
                
                fprintf(fid, ';; Region of interest within system volume (also in voxels) \n');
                fprintf(fid, 'image_x =  0  %3.0f       \n', job.ROIlimits(2,2)-1);
                fprintf(fid, 'image_y =  0  %3.0f       \n', job.ROIlimits(2,1)-1);
                fprintf(fid, 'image_z =  0  %3.0f       \n\n', job.ROIlimits(2,3)-1);
                
                % Optical properties
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 1 = gray matter \n',...
                    opt_ppts{iwl,1});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 2 = white matter \n',...
                    opt_ppts{iwl,2});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 3 = CSF \n',...
                    opt_ppts{iwl,3});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 4 = skull \n',...
                    opt_ppts{iwl,4});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 5 = scalp \n',...
                    opt_ppts{iwl,5});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 6 = perturbation \n',...
                    opt_ppts_perturb{iwl});

                if pve_cfg==1
                    fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 7 = gray matter and bold response \n',...
                    opt_ppts{iwl,1});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 8 = white matter and bold response \n',...
                    opt_ppts{iwl,2});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 9 = CSF and bold response \n',...
                    opt_ppts{iwl,3});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 10 = skull and bold response \n',...
                    opt_ppts{iwl,4});
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 11 = scalp and bold response \n',...
                    opt_ppts{iwl,5});                
                end
                
                % Detectors (all other optodes)
                Rp_rmiv = Pp_rmiv(:,[1:iP-1 iP+1:NS+ND]);
                Rwd_rmiv = Pwd_rmiv(:,[1:iP-1 iP+1:NS+ND]);
                Rp_rmiv = Rp_rmiv(1:3,:);
                
                temp_Rp_rmiv = Rp_rmiv(1,:);
                Rp_rmiv(1,:) = Rp_rmiv(2,:);
                Rp_rmiv(2,:) = temp_Rp_rmiv;
                
                Rwd_rmiv = Rwd_rmiv(1:3,:);
                Rr = r(1,[1:iP-1 iP+1:NS+ND]);
                fprintf(fid, '\ndetector { pos = [%4.0f %4.0f %4.0f]\n dir = [%5.4f %5.4f %5.4f]\n rad = %1.2f }\n',...
                    [Rp_rmiv;Rwd_rmiv;Rr]);
                fclose(fid);
                
                
            elseif algo==1 % MCX
                n_cfg = fullfile(mc_dir,[PNo '_' num2str(job.wl_dev(iwl)) 'nm.inp']);
                fid = fopen(n_cfg,'w');
                
                disp('Generating files for sensitivity profiles computation')
                
                % We create the file src.cfg for the Monte-Carlo simulation by reading
                % template.cfg and adding relevant information
                fprintf(fid,'%12.0f             # total photon (not used)\n', job.par.nphotons);
                fprintf(fid,'%12.0f            # RNG seed, negative to generate\n', job.par.seed);
                fprintf(fid,'%.1f %.1f %.1f              # source position (mm)\n',...
                    Pp_rmiv(1,iP), Pp_rmiv(2,iP), Pp_rmiv(3,iP));
                fprintf(fid,'%.3f %.3f %.3f                # initial directional vector\n',...
                    Pwd_rmiv(1,iP), Pwd_rmiv(2,iP), Pwd_rmiv(3,iP));
                end_time =  job.par.numTimeGates*job.par.deltaT;
                fprintf(fid,'%s %2.2e %1.2e   # time-gates(s): start, end, step\n','0.00e-9', end_time, job.par.deltaT);
                %fprintf(fid,'%s %s %s   # time-gates(s): start, end, step\n', '0','1e-9','1e-9');
                %[~,file,ext] = fileparts(job.n);%%% le .bin
                fprintf(fid,'%s     # volume (''uchar'' format)\n', n_b8i);%[file, ext]);
                
                fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # x: voxel size, dim, start/end indices\n',...
                    job.par.voxelSize,job.ROIlimits(2,1),1,job.ROIlimits(2,1));
                fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # y: voxel size, dim, start/end indices\n',...
                    job.par.voxelSize,job.ROIlimits(2,2),1,job.ROIlimits(2,2));
                fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # z: voxel size, dim, start/end indices\n',...
                    job.par.voxelSize,job.ROIlimits(2,3),1,job.ROIlimits(2,3));
                fprintf(fid, '%g %s \n',job.nummed,' # num of media'); %%%%%% a changer /////////
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,1}([2 3 1 4]),'# GM: scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,2}([2 3 1 4]),'# WM: scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,3}([2 3 1 4]),'# CSF: scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,4}([2 3 1 4]),'# Skull: scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,5}([2 3 1 4]),'# Scalp: scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts_perturb{iwl}([2 3 1 4]),'# Perturbation: scat(1/mm), g, mua (1/mm), n');
                if pve_cfg==1
                        fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,1}([2 3 1 4]),'# GM (BOLD): scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,2}([2 3 1 4]),'# WM (BOLD): scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,3}([2 3 1 4]),'# CSF (BOLD): scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,4}([2 3 1 4]),'# Skull (BOLD): scat(1/mm), g, mua (1/mm), n');
                fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                    opt_ppts{iwl,5}([2 3 1 4]),'# Scalp (BOLD): scat(1/mm), g, mua (1/mm), n');                
                end
                
                % Detectors (all other optodes)
                Rp_rmiv = Pp_rmiv(:,[1:iP-1 iP+1:NS+ND]);
                %             Rwd_rmiv = Pwd_rmiv(:,[1:iP-1 iP+1:NS+ND]);
                Rp_rmiv = Rp_rmiv(1:3,:);
                %             Rwd_rmiv = Rwd_rmiv(1:3,:);
                Rr = r(1,[1:iP-1 iP+1:NS+ND]);
                Rn = (1:size(Rr,2));
                fprintf(fid,'%g %g            # detector number and radius (mm)\n',...
                    size(Rr,2),r(1,NS+1));
                for iR=1:size(Rr,2)
                    fprintf(fid,'%.1f	%.1f	%.1f  # detector %g position (mm)\n',...
                        Rp_rmiv(:,iR),Rn(iR));
                end
                fclose(fid);
            end
        end
    end
end
end