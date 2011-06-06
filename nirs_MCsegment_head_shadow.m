function varargout = nirs_MCsegment_head_shadow(varargin)
% Build the head shadow : it sums all the ci images and smooth the result
% to provide a boolean map indicating the belonging of a voxel to the head
% volume.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-06

fname = varargin{1};
[dir,name] = fileparts(fname);

[dummy,V] = spm_imcalc_ui([strcat(dir,'\','c1',name,'.nii');...
    strcat(dir,'\','c2',name,'.nii');...
    strcat(dir,'\','c3',name,'.nii');...
    strcat(dir,'\','c4',name,'.nii');...
    strcat(dir,'\','c5',name,'.nii')],...
    strcat(dir,'\','head_shadow_ci.nii'),...
    'i1+i2+i3+i4+i5');

Y = spm_read_vols(V);

se_size_hs = varargin{2};
thresh_hs = varargin{3};

% once the sum is calculated an erosion is carried out to smoothen the
% image and to uniformise the mask
se = strel('disk',se_size_hs);

for i=1:size(Y,1)
    Y(i,:,:) = medfilt2(squeeze(Y(i,:,:)));
    Y(i,:,:) = imerode(squeeze(Y(i,:,:)),se);
end
% a threshold is used to achieve the correction of errors
Y(Y >= thresh_hs) = 1;
Y(Y < thresh_hs)  = 0;


V = struct('fname',fullfile(dir,'processed_head_shadow_ci.nii'),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);
V = spm_create_vol(V);
V = spm_write_vol(V, Y);

varargout{1} = strcat(dir,'\','processed_head_shadow_ci','.nii');