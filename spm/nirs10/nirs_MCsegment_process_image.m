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
        
    case 2
        %Preferred method
        
        %Example:
        %This aims to modify the skin (c5) layer so that it more accurately
        %represents the skin -- especially in order to re-attribute to skull
        %many voxels initially attributed to skin by SPM new segment

        %for each brain layer (especially skin and skull),
        % a 3D Gaussian filter is applied (and used as a weight), 
        % then medfilt2 and then
        % a threshold mask on the result
        
        %'gaussNdilate' devenu : gauss poids ;odifie median otsu
        % Gaussian Filter
        masqueGauss = zeros(gaussfilt_size,gaussfilt_size,gaussfilt_size);
        for i=1:gaussfilt_size
            for j=1:gaussfilt_size
                for k=1:gaussfilt_size
                    masqueGauss(i,j,k) = 1/((2*pi)^(1/2)*gaussfilt_sdev)*...
                        exp(-1/20*((floor(gaussfilt_size/2)+1-i)^2+...
                                   (floor(gaussfilt_size/2)+1-j)^2+...
                                   (floor(gaussfilt_size/2)+1-k)^2)/gaussfilt_sdev^2); %possible mistake in number 20?
                end
            end
        end
        %         %%%%%% A VERIFIER SI ON MET UNE NORME 2, L'INFLUENCE SUR LES
        %         %%%%%% TRAITEMENTS
        masqueGauss = masqueGauss/sum(masqueGauss(:));
        %This operation will not conserve probabilities across layers,
        %hence sum of probability over layers at one voxel can differ from unity
        % % % % % %         C = convn(Y,masqueGauss,'same'); %could replace by a convnfft for speed
        % % % % % %         Y2 = Y.*(C.^2);
        % % % % % %         %Median filter in direction x only in voxel space
        % % % % % %         for i=1:size(Y,1)
        % % % % % %             %median filter -- cartoon-like result
        % % % % % %             Y2(i,:,:) = medfilt2(squeeze(Y2(i,:,:)));
        % % % % % %             level = graythresh(squeeze(Y2(i,:,:)));
        % % % % % %             %create a boolean mask instead of probability distribution
        % % % % % %             Y(i,:,:) = im2bw(squeeze(Y2(i,:,:)),level);
        % % % % % %         end
        
        
        C = convn(Y,masqueGauss,'same'); %could replace by a convnfft for speed
        Y2 = Y.*(C.^2);
        %Median filter in direction x only in voxel space
        for i=1:size(Y,1)
            %median filter -- cartoon-like result
            Y2(i,:,:) = medfilt2(squeeze(Y2(i,:,:)));
            level = graythresh(squeeze(Y2(i,:,:)));
            %create a boolean mask instead of probability distribution
            Y(i,:,:) = im2bw(squeeze(Y2(i,:,:)),level);
        end
%         Y = Y.*(C.^2);

    case 3%'otsu'
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
            level = graythresh(squeeze(Y(i,:,:)));
            Y(i,:,:) = im2bw(squeeze(Y(i,:,:)),level);
        end
        
    case 4%median filter then otsu in different orientations...
        % extreme smoothing
        for i=1:size(Y,2)
            Y(:,i,:) = medfilt2(squeeze(Y(:,i,:)));
        end
        
        for i=1:size(Y,1)
            Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
            level = graythresh(squeeze(Y(i,:,:)));
            Y(i,:,:) = im2bw(squeeze(Y(i,:,:)),level);
        end
end
%Write Y
V_processed = nirs_create_vol(fullfile(dir,['processed_',name,'.nii']),...
                    V.dim, V.dt, V.pinfo, V.mat, Y);

varargout{1}=V_processed;