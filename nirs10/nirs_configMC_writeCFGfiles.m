function out = nirs_configMC_writeCFGfiles(job)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   writeCFGfiles.m                   %%%%%
%%%%% Create .cfg files for tMCimg Monte-Carlo simulations %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if CUDA_on
%     file_ext = '.inp';
% else
%     file_ext = '.cfg';
% end


% Write each line to the file...
algo = job.algo;
P = job.P;
NS = job.NS;
ND = job.ND;
Pvoid = job.Pvoid;
Pp_rmiv = job.P.p;
Pwd_rmiv = job.P.wd;
r = job.P.r;
mc_dir = job.mc_dir;
wl = job.wl;

for iP = 1:NS+ND
    if isempty(Pvoid) || ~Pvoid(1,iP) %skip optodes with no data &&&'nS',iS,wl(iwl)
        if iP<=NS
            if iP<10, PNo = ['S_No' num2str(0) num2str(iP)]; else PNo = ['S_No' num2str(iP)]; end
        else
            if iP<NS+10, PNo = ['D_No' num2str(0) num2str(iP-NS)]; else PNo = ['D_No' num2str(iP-NS)]; end
        end
        
        
        
        if algo==2%%%%%%%%%% tMCimg
            n_cfg = fullfile(mc_dir,[PNo '_' num2str(wl) 'nm.cfg']);
            fid = fopen(n_cfg,'w');
            % Title
            fprintf(fid, ';; Fichier de configuration\n');
            % ID subject
            fprintf(fid, [';; Sujet # : ' char(job.n_id) '\n\n']);
            % ROI limits relative to whole MRI
            fprintf(fid, ';; ROI limits relative to whole MRI\n');
            fprintf(fid, [';; xmin : ' num2str(job.ROIlimits(1,1)) '\n']);
            fprintf(fid, [';; xmax : ' num2str(job.ROIlimits(2,1)) '\n']);
            fprintf(fid, [';; ymin : ' num2str(job.ROIlimits(1,2)) '\n']);
            fprintf(fid, [';; ymax : ' num2str(job.ROIlimits(2,2)) '\n']);
            fprintf(fid, [';; zmin : ' num2str(job.ROIlimits(1,3)) '\n']);
            fprintf(fid, [';; zmax : ' num2str(job.ROIlimits(2,3)) '\n\n']);
            % Number of photons
            fprintf(fid, ';; Number of photons\n');
            fprintf(fid, 'nphotons = %12.0f\n\n', job.parameters.nphotons);
            % Seed random number generator
            fprintf(fid, ';; Seed random number generator\n');
            fprintf(fid, 'seed    =  %6.0f\n\n', job.parameters.seed);
            % Modulation frequency
            fprintf(fid, ';; Modulation frequency\n');
            fprintf(fid, 'freq    = %g \n\n', job.parameters.modulationFreq);
            % Optode (source)
            fprintf(fid, ';; Source position.  All lengths in mm.\n');
            % Here watch out! The tMCimg code is written in C/C++ thus matrix
            % encoding is row-major; in Matlab it is on the contrary column-major
            % thus for Matlab-generated files one must invert x and y
            % dimensions
            fprintf(fid, 'source { pos = [%4.0f %4.0f %4.0f]\n', Pp_rmiv(1,iP), Pp_rmiv(2,iP), Pp_rmiv(3,iP));%% ici pas besoin d'inverser les colonnes x et y
            fprintf(fid, '	 dir = [%5.4f %5.4f %5.4f] \n', Pwd_rmiv(1,iP), Pwd_rmiv(2,iP), Pwd_rmiv(3,iP));%% ici pas besoin d'inverser les colonnes x et y
            fprintf(fid, ' 	 rad = %1.2f       } \n\n', r(1,iP));
            
            % Time gates
            fprintf(fid, ';; Times are in seconds.\n');
            fprintf(fid, 'start_time = 0.00e-9    \n'); % Always start at t=0...
            fprintf(fid, 'gate_width = %3.2e     \n', job.parameters.deltaT);
            fprintf(fid, 'ngate      = %2.0f         \n\n', job.parameters.numTimeGates);
            
            % Name of segmentation file
            fprintf(fid, ';; Segmentation file to load\n');
            [~,file,ext] = fileparts(job.n);
            fprintf(fid, 'segfile  %s \n\n', [file ext]);
            
            % Voxel size
            fprintf(fid, ';; Voxel size.  Square voxels only, for now (dz = dy = dx)\n');
            fprintf(fid, 'dx     %1.0f \n\n', job.parameters.voxelSize);
            
            % System dimensions - see note above about inverting x and y dimensions
            fprintf(fid, ';; Size of system (in voxels) \n');
            fprintf(fid, 'nxvox =  %3.0f       \n', job.dim_rmiv(2));
            fprintf(fid, 'nyvox =  %3.0f       \n', job.dim_rmiv(1));
            fprintf(fid, 'nzvox =  %3.0f       \n\n', job.dim_rmiv(3));
            
            fprintf(fid, ';; Region of interest within system volume (also in voxels) \n');
            fprintf(fid, 'image_x =  0  %3.0f       \n', job.dim_rmiv(2)-1);
            fprintf(fid, 'image_y =  0  %3.0f       \n', job.dim_rmiv(1)-1);
            fprintf(fid, 'image_z =  0  %3.0f       \n\n', job.dim_rmiv(3)-1);
            
            % Optical properties
            fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 1 = gray matter \n',...
                job.parameters.gmPpties);
            fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 2 = white matter \n',...
                job.parameters.wmPpties);
            fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 3 = CSF \n',...
                job.parameters.csfPpties);
            fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 4 = skull \n',...
                job.parameters.skullPpties);
            fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 5 = scalp \n',...
                job.parameters.scalpPpties);
            if isfield(job.parameters,'perturbationPpties')
                fprintf(fid, 'tissue { mua = %5.4f mus = %3.2f g = %2.1f n = %3.2f } ; 6 = perturbation \n',...
                    job.parameters.perturbationPpties);
            end
            
            % Detectors (all other optodes)
            Rp_rmiv = Pp_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            Rwd_rmiv = Pwd_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            Rp_rmiv = Rp_rmiv(1:3,:);
            Rwd_rmiv = Rwd_rmiv(1:3,:);
            Rr = r(1,[1:iP-1 iP+1:NS+ND]);
            fprintf(fid, '\ndetector { pos = [%4.0f %4.0f %4.0f]\n dir = [%5.4f %5.4f %5.4f]\n rad = %1.2f }\n',...
                [Rp_rmiv;Rwd_rmiv;Rr]);
            fclose(fid);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif algo==1 %%%% MCX
            n_cfg = fullfile(mc_dir,[PNo '_' num2str(wl) 'nm.inp']);
            fid = fopen(n_cfg,'w');
            
            %%%
            disp('Generating files for sensitivity profiles computation')
            
            % We create the file src.cfg for the Monte-Carlo simulation by reading
            % template.cfg and adding relevant information
            fprintf(fid,'%12.0f             # total photon (not used)\n', job.parameters.nphotons);
            fprintf(fid,'%12.0f            # RNG seed, negative to generate\n', job.parameters.seed);
            fprintf(fid,'%.1f %.1f %.1f              # source position (mm)\n',...
                Pp_rmiv(2,iP), Pp_rmiv(1,iP), Pp_rmiv(3,iP));
            fprintf(fid,'%.1f %.1f %.1f                # initial directional vector\n',...
                Pwd_rmiv(2,iP), Pwd_rmiv(1,iP), Pwd_rmiv(3,iP));
            %end_time =  job.parameters.numTimeGates*job.parameters.deltaT;
            %fprintf(fid,'%s %2.2e %1.2e   # time-gates(s): start, end, step\n','0.00e-9', end_time, job.parameters.deltaT);
            fprintf(fid,'%s %s %s   # time-gates(s): start, end, step\n', '5e-11','5e-9','1e-10');
            [~,file,ext] = fileparts(job.n);%%% le .bin
            fprintf(fid,'%s     # volume (''uchar'' format)\n', [file, ext]);
            
            fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # x: voxel size, dim, start/end indices\n',...
                job.parameters.voxelSize,job.dim_rmiv(1),1,job.dim_rmiv(1));
            fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # y: voxel size, dim, start/end indices\n',...
                job.parameters.voxelSize,job.dim_rmiv(2),1,job.dim_rmiv(2));
            fprintf(fid,'%1.0f %3.0f %3.0f %3.0f            # z: voxel size, dim, start/end indices\n',...
                job.parameters.voxelSize,job.dim_rmiv(3),1,job.dim_rmiv(3));
            fprintf(fid, '%s %s \n','5',' # num of media');
            fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                job.parameters.gmPpties([2 3 1 4]),'# GM: scat(1/mm), g, mua (1/mm), n');
            fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                job.parameters.wmPpties([2 3 1 4]),'# WM: scat(1/mm), g, mua (1/mm), n');
            fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                job.parameters.csfPpties([2 3 1 4]),'# CSF: scat(1/mm), g, mua (1/mm), n');
            fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                job.parameters.skullPpties([2 3 1 4]),'# Skull: scat(1/mm), g, mua (1/mm), n');
            fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
                job.parameters.scalpPpties([2 3 1 4]),'# Scalp: scat(1/mm), g, mua (1/mm), n');
            
            % Detectors (all other optodes)
            Rp_rmiv = Pp_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            %             Rwd_rmiv = Pwd_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            Rp_rmiv = Rp_rmiv(1:3,:);
            %             Rwd_rmiv = Rwd_rmiv(1:3,:);
            Rr = r(1,[1:iP-1 iP+1:NS+ND]);
            Rn = (1:size(Rr,2));
            fprintf(fid,'%g %s            # detector number and radius (mm)\n',...
                size(Rr,2),'0.4');
            for iR=1:size(Rr,2)
                fprintf(fid,'%.1f	%.1f	%.1f  # detector %g position (mm)\n',...
                    Rp_rmiv(:,iR),Rn(iR));
            end
            fclose(fid);
            %%%
            
            %             % System dimensions - we no longer invert x and y dimensions, as MCX
            %             % has option to read with Matlab convention rather than C convention
            %             fprintf(fid, '%12.0f %s\n', job.parameters.nphotons, ' # total photon (not used)');
            %             fprintf(fid, '%12.0f %s\n', job.parameters.seed, ' # RNG seed, negative to generate');
            %             fprintf(fid, '%4.0f %4.0f %4.0f %s\n', Pp_rmiv(2,iP), Pp_rmiv(1,iP), Pp_rmiv(3,iP),' # source position (mm)');
            %             fprintf(fid, '%5.4f %5.4f %5.4f %s\n', Pwd_rmiv(2,iP), Pwd_rmiv(1,iP), Pwd_rmiv(3,iP),' # initial directional vector');
            %             %start_time = 0;
            %             end_time =  job.parameters.numTimeGates*job.parameters.deltaT;
            %             fprintf(fid, '%s %2.2e %1.2e %s\n', '0.00e-9', end_time, job.parameters.deltaT,' # time-gates(s): start, end, step');
            %             [~,file,~] = fileparts(job.n);
            %             fprintf(fid, '%s %s\n', file,' # volume (''uchar'' format)');
            %
            %             % System dimensions - see note about inverting x and y dimensions
            %             fprintf(fid, '%1.0f %3.0f %3.0f %3.0f %s\n', job.parameters.voxelSize, ...
            %                 job.dim_rmiv(1),0,job.dim_rmiv(1)-1,' # x: voxel size, dim, start/end indices');
            %             fprintf(fid, '%1.0f %3.0f %3.0f %3.0f %s\n', job.parameters.voxelSize, ...
            %                 job.dim_rmiv(2),0,job.dim_rmiv(2)-1,' # y: voxel size, dim, start/end indices');
            %             fprintf(fid, '%1.0f %3.0f %3.0f %3.0f %s\n', job.parameters.voxelSize, ...
            %                 job.dim_rmiv(3),0,job.dim_rmiv(3)-1,' # z: voxel size, dim, start/end indices');
            %             fprintf(fid, '%s %s \n','5',' # num of media');
            %             fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
            %                 job.parameters.gmPpties([2 3 1 4]),'# GM: scat(1/mm), g, mua (1/mm), n');
            %             fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
            %                 job.parameters.wmPpties([2 3 1 4]),'# WM: scat(1/mm), g, mua (1/mm), n');
            %             fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
            %                 job.parameters.csfPpties([2 3 1 4]),'# CSF: scat(1/mm), g, mua (1/mm), n');
            %             fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
            %                 job.parameters.skullPpties([2 3 1 4]),'# Skull: scat(1/mm), g, mua (1/mm), n');
            %             fprintf(fid, '%3.2f %2.1f %5.4f %3.2f %s\n', ...
            %                 job.parameters.scalpPpties([2 3 1 4]),'# Scalp: scat(1/mm), g, mua (1/mm), n');
            %
            %             % Detectors (all other optodes)
            %             Rp_rmiv = Pp_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            %             Rwd_rmiv = Pwd_rmiv(:,[1:iP-1 iP+1:NS+ND]);
            %             Rp_rmiv = Rp_rmiv(1:3,:);
            %             Rwd_rmiv = Rwd_rmiv(1:3,:);
            %             Rr = r(1,[1:iP-1 iP+1:NS+ND]);
            %             fprintf(fid, '\ndetector { pos = [%4.0f %4.0f %4.0f]\n dir = [%5.4f %5.4f %5.4f]\n rad = %1.2f }\n',...
            %                 [Rp_rmiv;Rwd_rmiv;Rr]);
            %             fclose(fid);
        end
    end
end
out =1;
end