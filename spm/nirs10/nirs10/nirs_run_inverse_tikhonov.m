function out = nirs_run_inverse_tikhonov(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: 

%Load NIRS.mat information
load(job.NIRSmat{1,1});

%Compute Tikhonov inversion
%use Delta mu = A^T (AA^T+alpha)^-1 y

%wavelengths
lambda = NIRS.Lambda; %= [830 690];

%Tikhonov regularization parameters
alpha = 1; %to be determined later

%Load measure list - code taken from spm_generate_sensitivity_matrix
%select channel pairs
%[t,sts] = spm_select(1,'.mat','Select SD structure (.mat file in dataSPM\data_all)',[],0);
%if ~sts, return; end
%temp = load(t);
SrcDet = NIRS.ml(:,1:2);
%clear temp
n_chn = size(SrcDet,1)/2;

selectedNIRfile = 1; %Which NIR file to apply the inversion to, if there are several sessions
%Load .nir data (converted from boxy)
t = NIRS.NIRfile{selectedNIRfile};
%[t,sts] = spm_select(1,'.nir','Select NIRS processed data (.nir file in dataSPM\data_all)',[],0);
%if ~sts, return; end

fid = fopen(t,'r');
meas = fread(fid, 'float32',0,'ieee-le'); 
fclose(fid);

t = NIRS.sens_name;

%Load sensitivity matrix for one wavelength
% [t,sts] = spm_select([1 inf],'.mat','Select sensitivity matrix for one wavelength (sens - .mat file in T1\dataMC\)',[],0);
% if ~sts, return; end
load(t); %very long

%See which wavelength the user has selected
[dir1 file1 ext1] = fileparts(t);
for Idx=1:2
   if strfind(file1,int2str(lambda(Idx))), wlength = Idx; end
end

meas = reshape(meas,2*n_chn,[]);
npt = size(meas,2);
%select first or second half of measures based on the chosen wavelength 
meas = meas(1+(wlength-1)*n_chn:wlength*n_chn,:);
%normalize measures
for Idx=1:n_chn
    meas(Idx,:) = meas(Idx,:)/meas(Idx,1) -ones(1,npt);
end

%Compute inverse -- call Edgar's spm_lot_tikh
ATinv_AAT = fsens' / (fsens * fsens' + alpha * eye(size(fsens,1)));
%save the inverse
save([dir1 '\ATinv_AAT_' int2str(lambda(wlength)) '.mat'],'ATinv_AAT','-v7.3'); 
%free up some memory
clear fsens

%or load ATinv_AAT if it's already been calculated - takes about 30s
[t,sts] = spm_select(1,'.mat','Load ATinv_AAT for one wavelength',[],0);
if ~sts, return; end
load(t);

%Need dimensions and details to save volume -- should be a separate function, since it is called for the 3rd time now,
%see also spm_generate_sensitivity_matrix and Generate_bin_cfg_or_inp_files
%select segmented volume 
t = NIRS.anatT1;

%[t,sts] = spm_select(1,'image','Select segmented volume only to get image dimensions (in T1 folder)',[],0);
%if ~sts, return; end
seg_vol = spm_vol(t);

%dimensions of image
dimi = seg_vol.dim;
inv_mat = spm_imatrix(seg_vol.mat); %parameters of inverse transform
scalings = diag(inv_mat(7:9)); %see spm_matrix for definition
%Dimension of volume
dimf = ceil((dimi-1) * scalings);
%Create linear spaces for initial volume (voxel) data and final (in mm)
xi = linspace(0,dimf(1)-1,dimi(1));
yi = linspace(0,dimf(2)-1,dimi(2));
zi = linspace(0,dimf(3)-1,dimi(3));
xf = 0:1:dimf(1)-1;
yf = 0:1:dimf(2)-1;
zf = 0:1:dimf(3)-1;
%Make grids, to use with 3-dimensional interpolation
[xi1,yi1,zi1] = meshgrid(xi,yi,zi);
[xf1,yf1,zf1] = meshgrid(xf,yf,zf);
    
%Apply the inverse to the measures

%save images - in nifti but using format .hdr + .img because of the large
%number of files
%This will take terabytes of data -- but let's just look at a few images to
%see if they have been sensibly reconstructed


for Idx=2400:60:2600 %npt 
    dopt = ATinv_AAT * meas(:,Idx);
    dopt = reshape(dopt,dimf);
    %annoying interpolation to do each time
    dopt_image = interp3(xf1,yf1,zf1,dopt,xi1,yi1,zi1,'nearest'); %takes about 20 sec to run
    save_each_pt = 1;
    if save_each_pt
        if Idx < 10
            nstr = ['00000' int2str(Idx)];
        else if Idx < 100
                nstr = ['0000' int2str(Idx)];
            else if Idx < 1000
                    nstr = ['000' int2str(Idx)];
                else if Idx < 10000
                        nstr = ['00' int2str(Idx)];
                    else if Idx < 100000
                            nstr = ['0' int2str(Idx)];
                        else nstr = int2str(Idx);
                        end
                    end
                end
            end
        end
        
        %Save in .nii format - use structure of segmented file
        V.fname = [dir1 '\nirs_img_' nstr '.nii']; 
        V.dim  = seg_vol.dim; %same as dimi
        V.dt   = [16 0]; %float32 not seg_vol.dt which is int16
        V.mat   = seg_vol.mat;
        V.pinfo = seg_vol.pinfo;
        spm_write_vol(V,dopt_image);
    end
end


%Code de Michèle:
dB = 20;
factor = 10.^(dB/10);
vollog = log(dopt);
mask = vollog > ( max(vollog(:)) / factor );
vol2log = vollog .* ( mask ) .* (dopt>0);
figure, imagesc(vollog(:,:,170)), colorbar, title('log(phi)')
figure, imagesc(mask(:,:,170)), colorbar, title('mask 60 dB')
figure, imagesc(vol2log(:,:,170)), colorbar
title('log(phi) masked at 60 dB')


% Positions projected : 
% ça, c'est bien en ordre x,y,z (il a inversé l'ordre des 2 premières
% coordonnées par rapport au fichier .cfg
% Donc, No10 est à 119(+1) (sur 176), à 69(+33) (sur 256), à 84+114
% No3 à 94, 89+33, 92+114 

% 2e paire (p8), No10 et No6, devrait être ensemble dans la tranche ~67+33.

% phi0 est ~phi6(119,102,198) soit ~10^-2.6. 




save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');