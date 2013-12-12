function nirs_load_channels_placement_figs(NIRSmat)
% Loads the figures with channels placement as subplots in a single figure.
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

load(NIRSmat)
if isfield(NIRS, 'jobCoreg')
    [pathName, fileName, fileExt] = fileparts(NIRS.jobCoreg.NIRSmat{1});
    if exist(fullfile(pathName,'extra_coreg'),'dir')
        pathName = fullfile(pathName,'extra_coreg');
    end
    h = figure; set(h,'color', 'w')
    load Split
    colormap(split);
    h = subplot(332);
    regionName = 'dorsal';
    getCoregFig(pathName, regionName, h);
    
    h = subplot(335);
    regionName = 'frontal';
    getCoregFig(pathName, regionName, h);
    
    h = subplot(334);
    regionName = 'left';
    getCoregFig(pathName, regionName, h);
    
    h = subplot(336);
    regionName = 'right';
    getCoregFig(pathName, regionName, h);
    
    h = subplot(337);
    regionName = 'occipital';
    getCoregFig(pathName, regionName, h);
    
    h = subplot(339);
    regionName = 'ventral';
    getCoregFig(pathName, regionName, h);
else
    error('nirs10:nirs_load_channels_placement_figs:noCoregJob', 'No coregistration job done!')
end

function getCoregFig(pathName, regionName, h) 
nirs_importfig(fullfile(pathName, regionName), h);
axis image;
axis off;
title(regionName,'FontSize',12);
% EOF
