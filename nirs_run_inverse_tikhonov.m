function out = nirs_run_inverse_tikhonov(job)
% Clément Bonnéry May 2011

% Le systeme qu'on resout est le suivant :
% on est en un UNIQUE POINT TEMPOREL donc
%              ---------------------
%
% / Y = Xsens*Dmua(vx) + epsilon_channel-noise
% \ Dmua(vx) = [ext]*[DHbO(vx);DHbR(vx)]
%
% on cherche a minimiser :
% argmin(||Y - Xsens*[ext]*beta||^2 + alpha*||beta - beta0||^2)
%
% Y                     = un instant donne dans le fichier .nirs
% Xsens                 = matrice sensitivite avec un point temporel
% [ext]                 = matrice des coefficients d extinctions
% Dmua                  = mu_a dans les voxels
% DHbO\R (beta)         = variation des concentrations dans les voxels
% alpha                 = hyperparametre de beta (DHbR\O)
% epsilon_channel-noise = bruit dans les canaux
if isfield(job.psel_choice,'all_points_downsampled')
    pmethod = 1;
    downfreq = job.psel_choice.all_points_downsampled.downsample_freq;
else
    pmethod = 0;
    %temporal points
    temp_pts = job.psel_choice.specific_points.temp_pts;
end
tikh_constraint = job.tikh_constraint;

for Idx=1:size(job.NIRSmat,1)
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'tikhOK') || job.force_redo)
            %selection of working directory
            if ~isempty(job.dir_in{1})
                f = job.dir_in;
                cs_dir = fileparts(f{1,:});
                [dummy cs_ldir] = fileparts(cs_dir); %name of simulation
                ics = 1;
                while ~strcmp(cs_ldir,NIRS.Cs.n{ics})
                    ics =ics+1;
                end
                %current simulation
                cs = NIRS.Cs.mcs{ics};
            else
                %take last simulation
                cs = NIRS.Cs.mcs{end};
                cs_dir = cs.dir;
            end
            
            %current tomo reconstruction
            ctm = {};
            % b8i is the bin file of the 8bits image segR
            % segR is ROI chosen in the segmented T1 image
            
            % b8i and segR have the same orientations as there is no
            % permutation (C and Matlab conventions in nirs_run_configMC)
            
            %resize segR so that it has same size as sensitivity matrix
            VsegR = spm_vol(cs.segR);
            clear jobR
            jobR.image_in ={cs.segR};
            jobR.out_dir = cs_dir;
            jobR.out_dim = [1 1 1];
            jobR.out_dt = 'same';
            jobR.out_vxsize = job.sens_vxsize;
            jobR.out_autonaming = 0;
            temp.segRR = nirs_resize(jobR);
            %put segRR in memory, calling it YsegRR, and used only for masks
            VsegRR = spm_vol(temp.segRR);
            YsegRR = spm_read_vols(VsegRR);
            YsegRR = reshape(YsegRR,[1 prod(VsegRR.dim)]);
            %         %masks c1 and c5
            %         m_c1 = zeros(size(YsegRR));
            %         m_c1(YsegRR==1)=1;% mask for GM
            %         m_c5 = zeros(size(YsegRR));
            %         m_c5(YsegRR==5)=1;% mask for skin
            %         if tikh_constraint
            %             m_c12 = zeros(size(YsegRR));
            %             m_c12(YsegRR==1 | YsegRR == 2)=1;
            %         end
            %         clear YsegRR;
            
            % Source detector pairs....
            C_cs = cs.C; %former Cmc
            NC_cs = length(C_cs); %Total number of measurements
            %load sensitivity matrix
            load(fullfile(cs_dir,'sens.mat'));
            Xmc = sens;
            clear sens
            % Need to downsample sensitivity matrix
            %         for Ci =1:NC_cs
            %             %inefficient, why write .nii
            %             Xmci = reshape(Xmc(Ci,:),VsegR.dim);
            %             Vmc = nirs_create_vol(fullfile(cs_dir,'Xmci.nii'),...
            %                 VsegR.dim, VsegR.dt, VsegR.pinfo, VsegR.mat, Xmci);
            %
            %             %this resizing has the same parameters as the one of segR to
            %             %segRR (lines 47 to 52)
            %             jobR.image_in ={Vmc.fname};
            %             jobR.out_autonaming = 1;
            %             jobR.out_prefix = 'R';
            %             fname =  nirs_resize(jobR);
            %             clear Xmci Vmc
            %             Vmc = spm_vol(fname);
            %             Ymc = spm_read_vols(Vmc);
            %             %
            %             XmcR(Ci,:) = reshape(Ymc,[1 prod(Vmc.dim)]);
            %         end
            
            dimR = floor(VsegR.dim/jobR.out_vxsize);
            Xmc_cm = zeros(NC_cs,prod(dimR));
            for Ci =1:NC_cs
                %inefficient, why write .nii
                Xmci = reshape(Xmc(Ci,:),VsegR.dim);
                tmp = nirs_resize_no_save(Xmci,dimR);
                Xmc_cm(Ci,:) = tmp(:);
            end
            
            Xmc = Xmc_cm;
            clear Xmc_cm Vmc Ymc
            Nvx = size(Xmc,2);
            NC2mi = NC_cs/2; %only for 2 wavelengths; number of pairs
            
            ext1 = GetExtinctions(NIRS.Cf.dev.wl(1));
            ext2 = GetExtinctions(NIRS.Cf.dev.wl(2));
            Ext = [ext1(1,1) ext1(1,2) ; ext2(1,1) ext2(1,2)];
            %regularization parameter
            alpha =job.alpha;
            
            Xwl{1} = sparse(Xmc(1:NC2mi,:)); %sensitivity matrix for 1st wavelength
            Xwl{2} = sparse(Xmc(NC2mi+1:end,:));
            
            %need to load data files prior to OD->concentration
            pdata = NIRS.Dt.fir.pp;
            wh = length(pdata);
            while wh
                [dummy fil1] = fileparts(pdata(wh).p{1});
                if strcmp(fil1(1),'h') || strcmp(fil1(1:2),'bh')
                    wh = wh-1;
                else
                    pdata = pdata(wh).p;
                    wh = 0;
                end
            end
            
            for ifnirs=1:length(pdata)
                Y = fopen_NIR(pdata{ifnirs},NIRS.Cf.H.C.N)';
                %creates Tmrs : tomographical reconstruction
                if ~isfield(NIRS,'Tm')
                    NIRS.Tm ={};
                end
                if pmethod
                    downstep = round(NIRS.Cf.dev.fs/downfreq);
                    TR = downstep/downfreq/NIRS.Cf.dev.fs;
                    temp_pts = 1:downstep:size(Y,1);
                end
                %mkdir
                if isfield(job.tikh_method,'tikhonov')% Tikhonov regularization
                    disp('Tikhonov inversion');
                    ctm.alg = 'Tikhinv';
                    jobSC.Cp=0;
                elseif isfield(job.tikh_method,'simple_bayes')% Interpretation Bayesienne simple
                    disp('pseudo inverse avec normes ponderes');
                    ctm.alg = 's_Bayes';
                    % IBS.1 : calcul des covariances
                    jobC.cov{1} = 'n';
                    jobC.cov{2} = 'p';
                    Cov = nirs_calculatecovariances(jobC);
                    jobSC.Cn = Cov{1};
                    jobSC.Cn = Cov{2};
                end
                
                daate = strrep(datestr(now),':','-');
                tm_dir = ['tm_' daate '_a' ctm.alg '_' int2str(job.sens_vxsize) 'mm'];
                ctm.p = fullfile(cs_dir,tm_dir);
                if ~exist(ctm.p,'dir'),mkdir(ctm.p);end
                if pmethod, save(fullfile(ctm.p,'TR.mat'),'TR'); end
                if isfield(NIRS.Tm,'tmrs')
                    itm = size(NIRS.Tm.tmrs,2)+1;
                else
                    itm=1;
                end
                ctm.n = tm_dir;
                ctm.hyperparameter = alpha;
                ctm.sens_vxsize = job.sens_vxsize;
                
                SC =0;
                if ~isempty(job.tikh_SC)
                    for i=1:size(job.tikh_SC,2)
                        if isfield(job.tikh_SC(i).tikh_mask,'wgmc')
                            jobSC.Y = YsegRR;
                            jobSC.tikh_mask = 'wgmc';
                            ctm.tikh_mask{1} = 'White and grey matter mask';
                            ctm.tikh_mask{2} = YsegRR;
                            
                        elseif isfield(job.tikh_SC(i).tikh_mask,'timask')
                            jobSC.Y = job.tikh_SC(i).tikh_mask.timask{:};
                            jobSC.tikh_mask = 'image';
                            ctm.tikh_mask{1} = 'User selected mask';
                            ctm.tikh_mask{2} = jobSC.Y;
                            
                        elseif isfield(job.tikh_SC(i).tikh_mask,'samcs')
                            jobSC.Y = cs.PVEmask;
                            jobSC.tikh_mask = 'image';
                            ctm.tikh_mask{1} = 'Same mask as for Monte-Carlo simulation';
                            ctm.tikh_mask{2} = jobSC.Y;
                        end
                        jobSC.alpha = alpha;
                        jobSC.alpha2 = job.tikh_SC(i).alpha2;
                        out = nirs_inverse_tikhonov_SC(jobSC);
                    end
                    SC = SC+out;
                end
                
                [dum,namRR,extRR] = fileparts(temp.segRR);
                ctm.segRR = fullfile(ctm.p,[namRR extRR]);
                copyfile(temp.segRR,ctm.segRR);
                %             delete(temp.segRR);
                NIRS.Tm.tmrs{itm} = ctm;
                NIRS.Tm.n{itm} = ctm.n;
                
                for itp=1:length(temp_pts)
                    %tic
                    Y_t0 = Y(temp_pts(itp),C_cs)';
                    Y_to{1} = Y_t0(1:NC2mi,1);
                    Y_to{2} = Y_t0(NC2mi+1:end,1);
                    
                    % Dmua
                    Dmua{1} = zeros(Nvx,1);
                    Dmua{2} = zeros(Nvx,1);
                    
                    switch ctm.alg
                        case 'Tikhinv' % Tikhonov regularization
                            for iwl=1:2
                                if itp==1
                                    XX{iwl} =Xwl{iwl}'*Xwl{iwl};
                                    sz = size(XX{iwl},2);
                                    %takes 27 GB of memory for 4e4 x 4e4 size
                                    if isempty(job.tikh_SC)
                                        XXLI{iwl} = sparse(XX{iwl} + alpha*eye(sz)); % eq.19 huppert_2010_hierarchical
                                    else
                                        XXLI{iwl} = sparse(XX{iwl} + SC);
                                    end
                                end
                                YY = (Xwl{iwl}'*Y_to{iwl});
                                Dmua{iwl} = XXLI{iwl} \ YY;
                            end
                        case 'sBayes' % Interpretation bayesienne simple
                            for iwl=1:2
                                if itp==1
                                    % IBS.2 : pseudo inversion (si necessaire on fait une svd)
                                    tXCn{iwl} = Xwl{iwl}'*Cn;
                                    tXCnX{iwl} =tXCn{iwl}*Xwl{iwl};
                                    sz = size(XX{iwl},2);
                                    %takes 27 GB of memory for 4e4 x 4e4 size
                                    if isempty(job.tikh_SC) % spatial constraint : mask
                                        XXLI{iwl} = sparse(tXCnX{iwl} + alpha*Cp*eye(sz)); % eq.19 huppert_2010_hierarchical
                                    else
                                        XXLI{iwl} = sparse(tXCnX{iwl} + SC);
                                    end
                                end
                                YY = (tXCn{iwl}'*Y_to{iwl});
                                Dmua{iwl} = XXLI{iwl} \ YY;
                            end
                    end
                    
                    beta = Ext \ [Dmua{1} Dmua{2}]'; %get HbO HbR from Mua
                    betaR_HbO = reshape(full(beta(1,:)),VsegRR.dim);
                    betaR_HbR = reshape(full(beta(2,:)),VsegRR.dim);
                    
                    %write nifti for DHbO DHbR
                    str0 = gen_num_str(itp,4);
                    V_O = nirs_create_vol(fullfile(ctm.p,['O_' str0 '.nii']),...
                        VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbO);
                    V_R = nirs_create_vol(fullfile(ctm.p,['R_' str0 '.nii']),...
                        VsegRR.dim, [16,0], VsegRR.pinfo, VsegRR.mat, betaR_HbR);
                    
                    %toc
                    %disp(int2str(itp));
                end
            end
            NIRS.flags.tikhOK = 1;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        try
            rmdir(ctm.p); %careful, the code might break here, so put breakpoint on line before
        end
        disp(['Could not run MonteCarlo reconstruction for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
%clear out %careful, variable out was used with another meaning previously!
out.NIRSmat = job.NIRSmat;

