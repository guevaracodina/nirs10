function genericROI = nirs_dfg_3D_roi
% % % ROIs % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nameROI         = cfg_entry;
nameROI.name    = 'ROI name';
nameROI.tag     = 'nameROI';
nameROI.strtype = 's';
nameROI.num     = [0 Inf];
nameROI.val     = {''};
nameROI.help    = {'Enter name for ROI. If left blank, ROIs will be enumerated.'}';

radiusROI         = cfg_entry;
radiusROI.name    = 'ROI radius value';
radiusROI.tag     = 'radiusROI';
radiusROI.strtype = 'r';
radiusROI.num     = [1 1];
radiusROI.val     = {5};
radiusROI.help    = {'Radius value in mm'}';

coordinateROI         = cfg_entry;
coordinateROI.name    = 'ROI coordinates';
coordinateROI.tag     = 'coordinateROI';
coordinateROI.strtype = 'r';
coordinateROI.num     = [1 3];
%coordinateROI.val     = {};
coordinateROI.help    = {'Enter MNI coordinates [x y z] in mm for center of ROI'}';

whichROI         = cfg_branch;
whichROI.tag     = 'whichROI';
whichROI.name    = 'Define ROI';
whichROI.val     = {nameROI coordinateROI radiusROI};
whichROI.help    = {'Define ROI'}';

genericROI         = cfg_repeat;
genericROI.tag     = 'genericROI';
genericROI.name    = 'Define ROIs';
genericROI.help    = {'Define here the ROIs to be analyzed'}';
genericROI.values  = {whichROI};
genericROI.num     = [1 Inf];