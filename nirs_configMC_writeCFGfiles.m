function out = nirs_configMC_writeCFGfiles(job)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   writeCFGfiles.m                   %%%%%
%%%%% Create .cfg files for tMCimg Monte-Carlo simulations %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write each line to the file...
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
        % thus for Matlab-generated files one must invert x and y dimensions
        fprintf(fid, 'source { pos = [%4.0f %4.0f %4.0f]\n', Pp_rmiv(2,iP), Pp_rmiv(1,iP), Pp_rmiv(3,iP));
        fprintf(fid, '	 dir = [%5.4f %5.4f %5.4f] \n', Pwd_rmiv(2,iP), Pwd_rmiv(1,iP), Pwd_rmiv(3,iP));
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
     end
end
out =1;
end