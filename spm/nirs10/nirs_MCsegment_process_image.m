function varargout = nirs_MCsegment_process_image(varargin)
% Runs the processing of an image
% FORMAT V_processed = spm_mc_process_image(file, method)
% file - image to be processed
% method - method used to process the image
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-06

file = varargin{1};
cimethod = varargin{2};

V = spm_vol(file);
Y = spm_read_vols(V);
[dir,name] = fileparts(V.fname);

se_size_pi = varargin{3};
se = strel('disk',se_size_pi);

gaussfilt_size = varargin{4};
gaussfilt_sdev = varargin{5};

switch cimethod
    
    case 0%'median'
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
        end
        
    case 1%'opening'
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
            Y(i,:,:) = imerode(squeeze(Y(i,:,:)),se);
            Y(i,:,:) = imdilate(squeeze(Y(i,:,:)),se);
        end
        
    case 2%'gaussNdilate'
        % Gaussian Filter
        masqueGauss = zeros(gaussfilt_size);
        for i=1:gaussfilt_size
            for j=1:gaussfilt_size
                masqueGauss(i,j) = 1/((2*pi)^(1/2)*gaussfilt_sdev)*exp(-1/20*((floor(gaussfilt_size/2)+1-i)^2+(floor(gaussfilt_size/2)+1-j)^2)/gaussfilt_sdev^2);
            end
        end
        %%%%%% A VERIFIER SI ON MET UNE NORME 2, L'INFLUENCE SUR LES
        %%%%%% TRAITEMENTS
        masqueGauss = masqueGauss/sum(sum(masqueGauss));
        
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
            Y_slice = squeeze(Y(i,:,:));
            Y(i,:,:) = conv2(Y_slice,masqueGauss,'same');
            Y(i,:,:) = imdilate(squeeze(Y(i,:,:)),se);
        end
        
    case 3%'otsu'
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
            level = graythresh(squeeze(Y(i,:,:)));
            Y(i,:,:) = im2bw(squeeze(Y(i,:,:)),level);
        end
end

V_processed = struct('fname',fullfile(dir,['processed_',name,'.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);
V_processed = spm_create_vol(V_processed);
V_processed = spm_write_vol(V_processed, Y);

varargout{1}=V_processed;