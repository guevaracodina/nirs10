function varargout = nirs_orthviews(action,varargin)
% This function origins from spm_orthviews
% NIRS version created for nirs10 toolbox
%
% Display orthogonal views of a set of images
% FORMAT H = nirs_orthviews('Image',filename[,position])
% filename - name of image to display
% area     - position of image {relative}
%            -  area(1) - position x
%            -  area(2) - position y
%            -  area(3) - size x
%            -  area(4) - size y
% H        - handle for ortho sections
%
% FORMAT nirs_orthviews('Redraw')
% Redraws the images
%
% FORMAT nirs_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel
%
% FORMAT nirs_orthviews('Space'[,handle[,M,dim]])
% handle   - the view to define the space by, optionally with extra
%            transformation matrix and dimensions (e.g. one of the blobs
%            of a view)
% with no arguments - puts things into mm space
%
% FORMAT nirs_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT nirs_orthviews('MaxBB')
% sets the bounding box big enough display the whole of all images
%
% FORMAT nirs_orthviews('Resolution'[,res])
% res      - resolution (mm)
% sets the sampling resolution for all images. The effective resolution
% will be the minimum of res and the voxel sizes of all images. If no
% resolution is specified, the minimum of 1mm and the voxel sizes of the
% images is used.
%
% FORMAT nirs_orthviews('Zoom'[,fov[,res]])
% fov      - half width of field of view (mm)
% res      - resolution (mm)
% sets the displayed part and sampling resolution for all images. The
% image display will be centered at the current crosshair position. The
% image region [xhairs-fov xhairs+fov] will be shown.
% If no argument is given or fov == Inf, the image display will be reset to
% "Full Volume". If fov == 0, the image will be zoomed to the bounding box
% from spm_get_bbox for the non-zero voxels of the image. If fov is NaN,
% then a threshold can be entered, and spm_get_bbox will be used to derive
% the bounding box of the voxels above this threshold.
% Optionally, the display resolution can be set as well.
%
% FORMAT nirs_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT nirs_orthviews('Reset')
% clears the orthogonal views
%
% FORMAT nirs_orthviews('Pos')
% returns the co-ordinate of the crosshairs in millimetres in the
% standard space.
%
% FORMAT nirs_orthviews('Pos', i)
% returns the voxel co-ordinate of the crosshairs in the image in the
% ith orthogonal section.
%
% FORMAT nirs_orthviews('Xhairs','off') OR nirs_orthviews('Xhairs')
% disables the cross-hairs on the display.
%
% FORMAT nirs_orthviews('Xhairs','on')
% enables the cross-hairs.
%
% FORMAT nirs_orthviews('Interp',hld)
% sets the hold value to hld (see spm_slice_vol).
%
% FORMAT nirs_orthviews('AddBlobs',handle,XYZ,Z,mat,name)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% name     - a name for this blob
% This method only adds one set of blobs, and displays them using a
% split colour table.
%
% FORMAT nirs_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% name     - a name for this blob
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% FORMAT nirs_orthviews('AddColourBar',handle,blobno)
% Adds colourbar for a specified blob set.
% handle   - image number
% blobno   - blob number
%
% FORMAT nirs_orthviews('RemoveBlobs',handle)
% Removes all blobs from the image specified by the handle(s).
%
% FORMAT nirs_orthviews('Addtruecolourimage',handle,filename,colourmap,prop,mx,mn)
% Adds blobs from an image in true colour.
% handle   - image number to add blobs to [default: 1]
% filename - image containing blob data [default: request via GUI]
% colourmap - colormap to display blobs in [default: GUI input]
% prop     - intensity proportion of activation cf grayscale [default: 0.4]
% mx       - maximum intensity to scale to [maximum value in activation image]
% mn       - minimum intensity to scale to [minimum value in activation image]
%
% FORMAT nirs_orthviews('Register',hReg)
% hReg     - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%
% FORMAT nirs_orthviews('AddContext',handle)
% handle   - image number to add context menu to
%
% FORMAT nirs_orthviews('RemoveContext',handle)
% handle   - image number to remove context menu from
%
% FORMAT nirs_orthviews('ZoomMenu',zoom,res)
% FORMAT [zoom, res] = nirs_orthviews('ZoomMenu')
% zoom     - A list of predefined zoom values
% res      - A list of predefined resolutions
% This list is used by spm_image and nirs_orthviews('addcontext',...) to
% create the 'Zoom' menu. The values can be retrieved by calling
% nirs_orthviews('ZoomMenu') with 2 output arguments. Values of 0, NaN and
% Inf are treated specially, see the help for nirs_orthviews('Zoom' ...).
%__________________________________________________________________________
%
% PLUGINS
% The display capabilities of nirs_orthviews can be extended with plugins.
% These are located in the nirs_orthviews subdirectory of the SPM
% distribution.
% The functionality of plugins can be accessed via calls to
% nirs_orthviews('plugin_name', plugin_arguments). For detailed descriptions
% of each plugin see help nirs_orthviews/spm_ov_'plugin_name'.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner et al
% $Id: nirs_orthviews.m 4276 2011-03-31 11:25:34Z spm $


% The basic fields of st are:
%         n        - the number of images currently being displayed
%         vols     - a cell array containing the data on each of the
%                    displayed images.
%         Space    - a mapping between the displayed images and the
%                    mm space of each image.
%         bb       - the bounding box of the displayed images.
%         centre   - the current centre of the orthogonal views
%         callback - a callback to be evaluated on a button-click.
%         xhairs   - crosshairs off/on
%         hld      - the interpolation method
%         fig      - the figure that everything is displayed in
%         mode     - the position/orientation of the sagittal view.
%                    - currently always 1
%
%         st.registry.hReg \_ See spm_XYZreg for documentation
%         st.registry.hMe  /
%
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "nirs_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
%
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which are documented here:
%
%         premul  - a matrix to premultiply the .mat field by.  Useful
%                   for re-orienting images.
%         window  - either 'auto' or an intensity range to display the
%                   image with.
%         mapping - Mapping of image intensities to grey values. Currently
%                   one of 'linear', 'histeq', loghisteq',
%                   'quadhisteq'. Default is 'linear'.
%                   Histogram equalisation depends on the image toolbox
%                   and is only available if there is a license available
%                   for it.
%         ax      - a cell array containing an element for the three
%                   views.  The fields of each element are handles for
%                   the axis, image and crosshairs.
%
%         blobs   - optional.  Is there for using to superimpose blobs.
%                   vol     - 3D array of image data
%                   mat     - a mapping from vox-to-mm (see spm_vol, or
%                             help on image formats).
%                   max     - maximum intensity for scaling to.  If it
%                             does not exist, then images are auto-scaled.
%
%                   There are two colouring modes: full colour, and split
%                   colour.  When using full colour, there should be a
%                   'colour' field for each cell element.  When using
%                   split colourscale, there is a handle for the colorbar
%                   axis.
%
%                   colour  - if it exists it contains the
%                             red,green,blue that the blobs should be
%                             displayed in.
%                   cbar    - handle for colorbar (for split colourscale).
%
% PLUGINS
% The plugin concept has been developed to extend the display capabilities
% of nirs_orthviews without the need to rewrite parts of it. Interaction
% between nirs_orthviews and plugins takes place
% a) at startup: The subfunction 'reset_st' looks for folders
%                'nirs_orthviews' in spm('Dir') and each toolbox
%                folder. Files with a name spm_ov_PLUGINNAME.m in any of
%                these folders will be treated as plugins.
%                For each such file, PLUGINNAME will be added to the list
%                st.plugins{:}.
%                The subfunction 'add_context' calls each plugin with
%                feval(['spm_ov_', st.plugins{k}], ...
%                  'context_menu', i, parent_menu)
%                Each plugin may add its own submenu to the context
%                menu.
% b) at redraw:  After images and blobs of st.vols{i} are drawn, the
%                struct st.vols{i} is checked for field names that occur in
%                the plugin list st.plugins{:}. For each matching entry, the
%                corresponding plugin is called with the command 'redraw':
%                feval(['spm_ov_', st.plugins{k}], ...
%                  'redraw', i, TM0, TD, CM0, CD, SM0, SD);
%                The values of TM0, TD, CM0, CD, SM0, SD are defined in the
%                same way as in the redraw subfunction of nirs_orthviews.
%                It is up to the plugin to do all necessary redraw
%                operations for its display contents. Each displayed item
%                must have set its property 'HitTest' to 'off' to let events
%                go through to the underlying axis, which is responsible for
%                callback handling. The order in which plugins are called is
%                undefined.

global st;
global opt_info;

persistent zoomlist;
persistent reslist;

if isempty(st), reset_st; end

if nargin == 0, action = ''; end

if ~any(strcmpi(action,{'reposition','pos'}))
    spm('Pointer','Watch');
end
    
switch lower(action)
    case 'image',
        H = specify_image(varargin{1});
        if ~isempty(H)
            if numel(varargin)>=2
                st.vols{H}.area = varargin{2};
            else
                st.vols{H}.area = [0 0 1 1];
            end
            if isempty(st.bb), st.bb = maxbb; end
            resolution;
            bbox;
            cm_pos;
        end
        varargout{1} = H;
        mmcentre     = mean(st.Space*[maxbb';1 1],2)';
        st.centre    = mmcentre(1:3);
        redraw_all
        
    case 'bb',
        if ~isempty(varargin) && all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end
        bbox;
        redraw_all;
        
    case 'redraw',
        redraw_all;
        eval(st.callback);
        if isfield(st,'registry'),
            spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        
    case 'reposition',
        if isempty(varargin), tmp = findcent;
        else tmp = varargin{1}; end
        if numel(tmp) == 3
            h = valid_handles(st.snap);
            if ~isempty(h)
                tmp = st.vols{h(1)}.mat * ...
                    round(st.vols{h(1)}.mat\[tmp(:); 1]);
            end
            st.centre = tmp(1:3);
        end
        redraw_all;
        eval(st.callback);
        if isfield(st,'registry')
            spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        cm_pos;
        %Record position
        opt_info.currentC = tmp(1:3);
        %Add text, Ke Peng
        clearErrorMsg;
        try
            if ~isempty(opt_info)
                new_dis = sqrt((opt_info.xSPM.XYZmm(1,:) - tmp(1)).^2 + ...
                        (opt_info.xSPM.XYZmm(2,:) - tmp(2)).^2 + ...
                        (opt_info.xSPM.XYZmm(3,:) - tmp(3)).^2);
                least_dis = min(new_dis);
                if least_dis <= 5
                    Fgraph = spm_figure('GetWin','Graphics');
                    try
                        delete(opt_info.hText.ht);
                    end
                    currentp = find(new_dis == least_dis);
                    currentXYZmm = opt_info.xSPM.XYZmm(:,currentp);
                    hTexthd = axes('Parent',Fgraph,'Position',[0.02 0.95 0.96 0.02],'Visible','off');
                    posstr = ['[' num2str(currentXYZmm(1)) ' ' num2str(currentXYZmm(2)) ' ' num2str(currentXYZmm(3)) ']'];
                    t_str = [opt_info.type ' ' opt_info.xSPM.label(currentp) posstr];
                    text(0.5,0,t_str,'Parent',hTexthd,'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',14);
                    opt_info.currentP = currentp;
                    opt_info.hText.ht = hTexthd;
                else
                    try
                        opt_info.currentP = [];
                        delete(opt_info.hText.ht);
                    end
                end
            end
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'dc_num'
        try
            clearErrorMsg;
            OPnum = str2num(get(findobj('Tag','ByNum'),'String')); %OPnum = get(opt_info.hText.hi,'String');
            DrawTmp = opt_info.xSPM.XYZmm(:,OPnum*2+2);
            nirs_orthviews('Reposition',DrawTmp);
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'dc_num_plus'
        try
            clearErrorMsg;
            OPnum = str2num(get(opt_info.hText.hi,'String'));
            OPnum = OPnum + 1;
            DrawTmp = opt_info.xSPM.XYZmm(:,OPnum*2+2);
            nirs_orthviews('Reposition',DrawTmp);
            set(opt_info.hText.hi,'String',int2str(OPnum));
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'dc_num_minus'
        try
            clearErrorMsg;
            OPnum = str2num(get(opt_info.hText.hi,'String'));
            OPnum = OPnum - 1;
            DrawTmp = opt_info.xSPM.XYZmm(:,OPnum*2+2);
            nirs_orthviews('Reposition',DrawTmp);
            set(opt_info.hText.hi,'String',int2str(OPnum));
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'goleftpa',
        try
            clearErrorMsg;
            DrawTmp = opt_info.xSPM.XYZmm(:,2);
            nirs_orthviews('Reposition',DrawTmp);
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'gonasion',
        try
            clearErrorMsg;
            DrawTmp = opt_info.xSPM.XYZmm(:,1);
            nirs_orthviews('Reposition',DrawTmp);
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'gorightpa',
        try
            clearErrorMsg;
            DrawTmp = opt_info.xSPM.XYZmm(:,3);
            nirs_orthviews('Reposition',DrawTmp);
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'correspondence',
        %Added by Ke Peng, called in nirs_run_liom_orth_coreg
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if ~isempty(opt_info.currentP)
                currentP = opt_info.currentP;
                if currentP ~= 1 && currentP ~= 2 && currentP ~= 3
                    %Find the corresponding point, either left one or right
                    %one
                    p0.XYZmm = opt_info.xSPM.XYZmm(:,currentP);
                    p0.label = opt_info.xSPM.label{currentP};
                    p0str = p0.label(1:(end-4));
                    
                    p1.XYZmm = opt_info.xSPM.XYZmm(:,currentP+1);
                    p1.label = opt_info.xSPM.label{currentP+1};
                    p1str = p1.label(1:(end-4));                    
                    
                    p2.XYZmm = opt_info.xSPM.XYZmm(:,currentP-1);
                    
                    if strcmp(p0str,p1str)
                        DrawTmp = p1.XYZmm(1:3,:);
                    else
                        DrawTmp = p2.XYZmm(1:3,:);
                    end                   
                    nirs_orthviews('Reposition',DrawTmp);
                end
            else
                text(0.5,0,'ERROR: Please choose an optode/a channel first','Parent',hErrTexthd,...
                'HorizontalAlignment','center',...
                'VerticalAlignment','baseline',...
                'FontWeight','Bold','FontSize',14);
                opt_info.hText.he = hErrTexthd;
            end
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregop',
        %Added by Ke Peng, called in nirs_run_liom_orth_coreg
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if ~isempty(opt_info.currentP) && opt_info.currentP ~= 1 && opt_info.currentP ~= 2 && opt_info.currentP ~= 3
                currentP = opt_info.currentP;
                Plabel = opt_info.xSPM.label(currentP);
                SlpOrCtx = Plabel{1};
                if currentP ~= 1 && currentP ~= 2 && currentP ~= 3 && strcmp(SlpOrCtx(end-2 : end),'Slp')
                    set(opt_info.hText.hb.crgop, 'Enable', 'on');
                    text(0.5,0,'Please choose a new position on the scalp for this optode/channel','Parent',hErrTexthd,...
                        'HorizontalAlignment','center','VerticalAlignment','baseline','FontWeight','Bold','FontSize',14);
                    opt_info.RCP = currentP; %Optode to re-coreg
                    set(opt_info.hText.hb.op, 'Enable', 'Off');
                    set(opt_info.hText.hb.crgop, 'Enable', 'On');
                    set(opt_info.hText.hb.fd, 'Enable', 'Off');
                    setallbutton(0);
                else
                    text(0.5,0,'ERROR: Please click on an optode/a channel ON THE SCALP','Parent',hErrTexthd,...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','baseline',...
                        'FontWeight','Bold','FontSize',14);
                end
            else
                text(0.5,0,'ERROR: Please click on an optode/a channel on the scalp','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
            end
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregop_exe'
        %Added by Ke Peng, called in nirs_run_liom_orth_coreg
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info, 'currentC') && isfield(opt_info, 'RCP')
                text(0.5,0,'Co-registration ongoing...','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                Preplace = opt_info.RCP;
                newp = opt_info.currentC;
                Q = opt_info.recoreg.Q;
                newpw = Q * [newp;1];
                newpw_c1 = projection_CS(newpw); %Korean's method
                newp_c1 = Q\newpw_c1;
                newp_c1 = newp_c1(1:3,:);
                opt_info.OPrecoreg.idx = [opt_info.OPrecoreg.idx (Preplace-2)/2]; %the ith Optode or channel
                opt_info.OPrecoreg.SlpXYZmm = opt_info.xSPM.XYZmm(:,Preplace);
                opt_info.OPrecoreg.CtxXYZmm = opt_info.xSPM.XYZmm(:,Preplace + 1);
                opt_info.OPrecoreg.undo = [0 opt_info.OPrecoreg.undo];
                opt_info.xSPM.XYZmm(:,Preplace) = newp; %replace the coordinates on the scalp with the new one
                opt_info.xSPM.XYZmm(:,Preplace + 1) = newp_c1; %replace the coordinates on the cortex with the new one
                nirs_orthdraw(opt_info.xSPM, opt_info.frender);
                delete(hErrTexthd);
                hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');
                text(0.5,0,'Co-registration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                nirs_orthviews('Reposition',newp);
                set(opt_info.hText.hb.crgopu, 'Enable', 'on');             
            else
                text(0.5,0,'ERROR: coordinates not recorded. Please re-do.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);               
            end
            set(opt_info.hText.hb.crgop, 'Enable', 'off');
            set(opt_info.hText.hb.op, 'Enable', 'off');
            set(opt_info.hText.hb.crgopu, 'Enable', 'on');
            set(opt_info.hText.hb.crgopr, 'Enable', 'off');
            set(opt_info.hText.hb.crgopcf, 'Enable', 'on');
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregop_undo',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info, 'OPrecoreg')
                P = opt_info.OPrecoreg.idx(end);
                Cs = opt_info.OPrecoreg.SlpXYZmm;
                Cc = opt_info.OPrecoreg.CtxXYZmm;
                opt_info.OPrecoreg.SlpXYZmm = opt_info.xSPM.XYZmm(:,P*2+2);
                opt_info.OPrecoreg.CtxXYZmm = opt_info.xSPM.XYZmm(:,P*2+3); %Record the old coordinates for undo
                opt_info.xSPM.XYZmm(:,P*2+2) = Cs; %replace the coordinates on the scalp with the new one
                opt_info.xSPM.XYZmm(:,P*2+3) = Cc; %replace the coordinates on the cortex with the new one
                opt_info.OPrecoreg.undo(end) = 1;
                nirs_orthdraw(opt_info.xSPM, opt_info.frender);
                text(0.5,0,'Undo coregistration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                nirs_orthviews('Reposition',Cs);
                set(opt_info.hText.hb.crgopu, 'Enable', 'off');
                set(opt_info.hText.hb.crgopr, 'Enable', 'on');
            else
                text(0.5,0,'ERROR: Undo optode not found.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
            end
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregop_redo'
        try
             [Fgraph hErrTexthd] = clearErrorMsg;
             if isfield(opt_info, 'OPrecoreg')
                P = opt_info.OPrecoreg.idx(end);
                Cs = opt_info.OPrecoreg.SlpXYZmm;
                Cc = opt_info.OPrecoreg.CtxXYZmm;
                opt_info.OPrecoreg.SlpXYZmm = opt_info.xSPM.XYZmm(:,P*2+2);
                opt_info.OPrecoreg.CtxXYZmm = opt_info.xSPM.XYZmm(:,P*2+3); %Record the old coordinates for undo
                opt_info.xSPM.XYZmm(:,P*2+2) = Cs; %replace the coordinates on the scalp with the new one
                opt_info.xSPM.XYZmm(:,P*2+3) = Cc; %replace the coordinates on the cortex with the new one
                opt_info.OPrecoreg.undo(end) = 0;
                nirs_orthdraw(opt_info.xSPM, opt_info.frender);
                text(0.5,0,'Redo coregistration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);                
                nirs_orthviews('Reposition',Cs);                 
             else
                text(0.5,0,'ERROR: Redo optode not found.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
             end
             set(opt_info.hText.hb.crgopu, 'Enable', 'on');
             set(opt_info.hText.hb.crgopr, 'Enable', 'off');
             opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregop_confirm'
        try
             [Fgraph hErrTexthd] = clearErrorMsg;
             if isfield(opt_info, 'OPrecoreg')
                 P = opt_info.OPrecoreg.idx(end);
                 type = opt_info.type;
                 recoreg = opt_info.recoreg;
                 m2v = opt_info.m2v;
                 Q = opt_info.recoreg.Q;
                 CSlp = opt_info.xSPM.XYZmm(:,P*2+2);
                 CCtx = opt_info.xSPM.XYZmm(:,P*2+3);
                 undo = opt_info.OPrecoreg.undo(end);
                 if ~undo
                     recoreg = FillOptodePos(recoreg,Q,m2v,P,type,CSlp,CCtx,opt_info.render_template);
                 end
                 opt_info.recoreg = recoreg;
                 text(0.5,0,'Recoreg op position confirmed.','Parent',hErrTexthd,...
                     'HorizontalAlignment','center',...
                     'VerticalAlignment','baseline',...
                     'FontWeight','Bold','FontSize',14);
             else
                text(0.5,0,'ERROR: Recoreg optode not found.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
             end
             set(opt_info.hText.hb.crgopu, 'Enable', 'off');
             set(opt_info.hText.hb.crgopr, 'Enable', 'off');
             set(opt_info.hText.hb.op, 'Enable', 'on');
             set(opt_info.hText.hb.fd, 'Enable', 'on');
             set(opt_info.hText.hb.crgopcf, 'Enable', 'off');
             setallbutton(1);
             opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'coregfd',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if ~isempty(opt_info.currentP) && (opt_info.currentP == 1 || opt_info.currentP == 2 || opt_info.currentP == 3)
                set(opt_info.hText.hb.crgfd, 'Enable', 'on');
                text(0.5,0,'Coreg will based on the latest saved NIRS. Please choose a new position.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                opt_info.FCP = opt_info.currentP; %Optode to re-coreg
                set(opt_info.hText.hb.crgfd, 'Enable', 'On');
                set(opt_info.hText.hb.fd, 'Enable', 'Off');
                set(opt_info.hText.hb.op, 'Enable', 'Off');
                setallbutton(0);
            else
                text(0.5,0,'ERROR: Please click on a fiducial point','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
            end
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            opt_info.hText.he = hErrTexthd;
        end
        
    case 'coregfd_exe',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info, 'currentC') && isfield(opt_info, 'FCP')
                set(opt_info.hText.hb.crgfd, 'Enable', 'off');
                set(opt_info.hText.hb.crgfdr, 'Enable', 'off');
                set(opt_info.hText.hb.crgfdu, 'Enable', 'off');
                text(0.5,0,'Co-registration ongoing...Please wait...','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                newp = opt_info.currentC;
                fcp = opt_info.FCP;
                opt_info.FDrecoreg.xSPM = opt_info.xSPM; %Record for undo
                opt_info.FDrecoreg.recoreg = opt_info.recoreg; %Record for undo
                opt_info.xSPM.XYZmm(:,fcp) = newp; %Replace with the new location
                Q = opt_info.recoreg.Q;
                fd_omm = opt_info.recoreg.corr.F.r.o.mm.p; % Fiducial coordinates in mm (original, unnormalised)
                Sp_rom = opt_info.recoreg.corr.S.r.o.mm.p;
                Dp_rom = opt_info.recoreg.corr.D.r.o.mm.p;
                Ns = opt_info.recoreg.corr.S.N;
                Pvoid = opt_info.recoreg.corr.P.void;
                Cid = opt_info.recoreg.corr.C.id;
                OPidx = opt_info.OPrecoreg.idx;
                OPundo = opt_info.OPrecoreg.undo;
                if ~isempty(OPidx)
                    OPidx = unique(OPidx(find(OPundo == 0)));%Find optodes/channels that have been aleady recoregistered.
                else
                    OPXYZmm_ed = [];
                end
                type = opt_info.type;
                render_template = opt_info.render_template;
                m2v = opt_info.m2v;
                if render_template
                    fd_wmm = opt_info.xSPM.XYZmm(:,1:3);
                    fd_rmm = Q\[fd_wmm; [1 1 1]];
                    fd_rmm = fd_rmm(1:3,:);
                else
                    fd_rmm = opt_info.xSPM.XYZmm(:,1:3);
                    fd_wmm = Q*[fd_rmm; [1 1 1]]; %Fiducial coordinates in mm
                    fd_wmm = fd_wmm(1:3,:); % Fiducial coordinates in mm (MNI, normalised)
                end
                switch type
                    case 'Channel'
                        sel = 0;
                        if ~isempty(OPidx)
                            OPXYZmm_ed = opt_info.recoreg.corr.C.r.m.mm.fp(:,OPidx);
                        end
                    case 'Source'
                        sel = 1;
                        if ~isempty(OPidx)
                            OPXYZmm_ed = opt_info.recoreg.corr.P.r.m.mm.fp(:,OPidx);
                        end
                    case 'Detector'
                        sel = 2;
                        if ~isempty(OPidx)
                            OPXYZmm_ed = opt_info.recoreg.corr.P.r.m.mm.fp(:,Ns+OPidx);
                        end
                    case 'All Channels/Sources/Detectors'
                        sel = 3;
                    otherwise
                        sel = -1;
                end
                recoreg1 = LIOMrecoreg(fd_omm,fd_rmm,Q,Sp_rom,Dp_rom,Pvoid,Cid,render_template,m2v,OPidx,OPXYZmm_ed,type);
                recoreg1.corr.F.w.m.mm.p = fd_wmm;
                opt_info.recoreg = recoreg1;
                [coreg xSPM opt] = nirs_orthPrep(recoreg1,sel,render_template,0,opt_info.m2v);
                opt_info.xSPM = xSPM;
                nirs_orthdraw(xSPM, opt_info.frender);
                delete(hErrTexthd);
                hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');
                text(0.5,0,'Co-registration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                nirs_orthviews('Reposition',newp);
                set(opt_info.hText.hb.crgfdu, 'Enable', 'on');
                set(opt_info.hText.hb.crgfdcf, 'Enable', 'on');
            else
                try
                    delete(hErrTexthd);
                end
                hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');
                text(0.5,0,'ERROR: coordinates not recorded. Please re-do.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                set(opt_info.hText.hb.crgfd, 'Enable', 'on');
                set(opt_info.hText.hb.crgfdu, 'Enable', 'off');
                set(opt_info.hText.hb.crgfdr, 'Enable', 'off');
                set(opt_info.hText.hb.crgfdcf, 'Enable', 'off');
            end

            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            set(opt_info.hText.hb.crgfd, 'Enable', 'off');
            opt_info.hText.he = hErrTexthd;
        end
        
    case 'coregfd_undo',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info, 'FDrecoreg') && ~isempty(opt_info.FDrecoreg)
                xSPM = opt_info.FDrecoreg.xSPM;
                recoreg = opt_info.FDrecoreg.recoreg;
                opt_info.FDrecoreg.xSPM = opt_info.xSPM;
                opt_info.FDrecoreg.recoreg = opt_info.recoreg;
                opt_info.xSPM = xSPM;
                opt_info.recoreg = recoreg;
                nirs_orthdraw(xSPM, opt_info.frender);
                delete(hErrTexthd);
                hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');
                text(0.5,0,'Undo FD co-registration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                newp = opt_info.xSPM.XYZmm(:,opt_info.FCP);
                nirs_orthviews('Reposition',newp);
                set(opt_info.hText.hb.crgfdr, 'Enable', 'on');
                set(opt_info.hText.hb.crgfdu, 'Enable', 'off');
            else
                text(0.5,0,'ERROR: Undo oldNIRS not found.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
            end
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            opt_info.hText.he = hErrTexthd;
            set(opt_info.hText.hb.fd, 'Enable', 'on');
        end
        
    case 'coregfd_redo',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info, 'FDrecoreg') && ~isempty(opt_info.FDrecoreg)
                xSPM = opt_info.FDrecoreg.xSPM;
                recoreg = opt_info.FDrecoreg.recoreg;
                opt_info.FDrecoreg.xSPM = opt_info.xSPM;
                opt_info.FDrecoreg.recoreg = opt_info.recoreg;
                opt_info.xSPM = xSPM;
                opt_info.recoreg = recoreg;
                nirs_orthdraw(xSPM, opt_info.frender);
                delete(hErrTexthd);
                hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');
                text(0.5,0,'Redo FD co-registration Done','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                newp = opt_info.xSPM.XYZmm(:,opt_info.FCP);
                nirs_orthviews('Reposition',newp);
                set(opt_info.hText.hb.crgfdr, 'Enable', 'off');
                set(opt_info.hText.hb.crgfdu, 'Enable', 'on');
            else
                text(0.5,0,'ERROR: Undo oldNIRS not found.','Parent',hErrTexthd,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
            end
            opt_info.hText.he = hErrTexthd;
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            opt_info.hText.he = hErrTexthd;
            set(opt_info.hText.hb.fd, 'Enable', 'on');
        end
        
    case 'coregfd_confirm',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            if isfield(opt_info.recoreg, 'err')
                opt_info.recoreg.pro.errValofCoreg_mm2 = opt_info.recoreg.err.errValofCoreg_mm2;
                opt_info.recoreg.pro.errValofCoreg_mm2_all = opt_info.recoreg.err.errValofCoreg_mm2_all;
                opt_info.recoreg.pro.errValofCoreg_cm_worst = opt_info.recoreg.err.errValofCoreg_cm_worst;
            end
            opt_info.FDcoreg = [];
            text(0.5,0,'Recoreg fd position confirmed.','Parent',hErrTexthd,...
                'HorizontalAlignment','center',...
                'VerticalAlignment','baseline',...
                'FontWeight','Bold','FontSize',14);
            set(opt_info.hText.hb.fd, 'Enable', 'on');
            set(opt_info.hText.hb.op, 'Enable', 'on');
            set(opt_info.hText.hb.crgfdcf, 'Enable', 'off');   
            set(opt_info.hText.hb.crgfdr, 'Enable', 'off');
            set(opt_info.hText.hb.crgfdu, 'Enable', 'off');
            setallbutton(1);
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
        end
        
    case 'orthsave',
        try
            [Fgraph hErrTexthd] = clearErrorMsg;
            opt_info.recoreg.corr.n = opt_info.NIRS.Cf.H.n;
            opt_info.recoreg.corr.C.n = opt_info.NIRS.Cf.H.C.n;
            opt_info.recoreg.corr.C.wl = opt_info.NIRS.Cf.H.C.wl;
            opt_info.recoreg.corr.C.gp = opt_info.NIRS.Cf.H.C.gp;
            %opt_info.recoreg.corr.P.n = opt_info.NIRS.Cf.H.P.n;
            opt_info.recoreg.pro.extcoeff_ref = opt_info.NIRS.Dt.pro.extcoeff_ref;
            opt_info.NIRS.Cf.H = opt_info.recoreg.corr;
            opt_info.NIRS.Dt.pro = opt_info.recoreg.pro;
            NIRS = opt_info.NIRS;
            save(opt_info.newNIRSlocation,'NIRS');
            [dir_coreg,dummy] = fileparts(opt_info.newNIRSlocation);
            ch_MNI_vx = opt_info.recoreg.corr.C.r.m.vx.c1.p;
            if isfield(opt_info.recoreg.corr.C, 'w')
                ch_MNIw_vx = opt_info.recoreg.corr.C.w.m.vx.c1.p;
            else
                ch_MNIw_vx = opt_info.m2v.w\(opt_info.recoreg.Q*opt_info.recoreg.corr.C.r.m.mm.c1.p);
            end
            use_fSeg = 1;
            if isfield(opt_info.NIRS.Dt.ana,'T1seg') && use_fSeg
                fSeg = opt_info.NIRS.Dt.ana.T1seg;
            else
                fSeg = [];
            end
            rendered_MNI = render_MNI_coordinates_new(ch_MNIw_vx,...
                ch_MNI_vx,opt_info.wT1_info,NIRS.Dt.ana.wT1.VF,opt_info.render_template,fSeg,dir_coreg,0,opt_info.ERad);
            rend_file = fullfile(dir_coreg,'TopoData.mat');
            recoreg_file = fullfile(dir_coreg,'opt_info.mat');
            save(rend_file, 'rendered_MNI');
            save(recoreg_file, 'opt_info');
            nirs_brain_project_2d(NIRS,dir_coreg,rendered_MNI,[],'r','','',[],0); %display
            clear global opt_info
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            opt_info.hText.he = hErrTexthd;
        end
     
        
    case 'setcoords',
        st.centre = varargin{1};
        st.centre = st.centre(:);
        redraw_all;
        eval(st.callback);
        cm_pos;
        
    case 'space',
        if numel(varargin)<1
            st.Space = eye(4);
            st.bb = maxbb;
            resolution;
            bbox;
            redraw_all;
        else
            space(varargin{:});
            resolution;
            bbox;
            redraw_all;
        end
        
    case 'maxbb',
        st.bb = maxbb;
        bbox;
        redraw_all;
        
    case 'resolution',
        resolution(varargin{:});
        bbox;
        redraw_all;
        
    case 'window',
        if numel(varargin)<2
            win = 'auto';
        elseif numel(varargin{2})==2
            win = varargin{2};
        end
        for i=valid_handles(varargin{1})
            st.vols{i}.window = win;
        end
        redraw(varargin{1});
        
    case 'delete',
        my_delete(varargin{1});
        
    case 'move',
        move(varargin{1},varargin{2});
        % redraw_all;
        
    case 'reset',
        my_reset;
        
    case 'pos',
        if isempty(varargin)
            H = st.centre(:);
        else
            H = pos(varargin{1});
        end
        varargout{1} = H;
        
    case 'interp',
        st.hld = varargin{1};
        redraw_all;
        
    case 'xhairs',
        xhairs(varargin{1});
        
    case 'register',
        register(varargin{1});
        
    case 'addblobs',
        addblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'addcolouredblobs',
        addcolouredblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'addimage',
        addimage(varargin{1}, varargin{2});
        % redraw(varargin{1});
        
    case 'addcolouredimage',
        addcolouredimage(varargin{1}, varargin{2},varargin{3});
        % redraw(varargin{1});
        
    case 'addtruecolourimage',
        if nargin < 2
            varargin(1) = {1};
        end
        if nargin < 3
            varargin(2) = {spm_select(1, 'image', 'Image with activation signal')};
        end
        if nargin < 4
            actc = [];
            while isempty(actc)
                actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
            end
            varargin(3) = {actc};
        end
        if nargin < 5
            varargin(4) = {0.4};
        end
        if nargin < 6
            actv = spm_vol(varargin{2});
            varargin(5) = {max([eps maxval(actv)])};
        end
        if nargin < 7
            varargin(6) = {min([0 minval(actv)])};
        end
        
        addtruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
            varargin{5}, varargin{6});
        % redraw(varargin{1});
        
    case 'addcolourbar',
        addcolourbar(varargin{1}, varargin{2});
        
    case {'removeblobs','rmblobs'},
        rmblobs(varargin{1});
        redraw(varargin{1});
        
    case 'addcontext',
        if nargin == 1
            handles = 1:24;
        else
            handles = varargin{1};
        end
        addcontexts(handles);
        
    case {'removecontext','rmcontext'},
        if nargin == 1
            handles = 1:24;
        else
            handles = varargin{1};
        end
        rmcontexts(handles);
        
    case 'context_menu',
        c_menu(varargin{:});
        
    case 'valid_handles',
        if nargin == 1
            handles = 1:24;
        else
            handles = varargin{1};
        end
        varargout{1} = valid_handles(handles);

    case 'zoom',
        zoom_op(varargin{:});
        
    case 'zoommenu',
        if isempty(zoomlist)
            zoomlist = [NaN 0 5    10  20 40 80 Inf];
            reslist  = [1   1 .125 .25 .5 .5 1  1  ];
        end
        if nargin >= 3
            if all(cellfun(@isnumeric,varargin(1:2))) && ...
                    numel(varargin{1})==numel(varargin{2})
                zoomlist = varargin{1}(:);
                reslist  = varargin{2}(:);
            else
                warning('nirs_orthviews:zoom',...
                        'Invalid zoom or resolution list.')
            end
        end
        if nargout > 0
            varargout{1} = zoomlist;
        end
        if nargout > 1
            varargout{2} = reslist;
        end
        
    otherwise,
        addonaction = strcmpi(st.plugins,action);
        if any(addonaction)
            feval(['spm_ov_' st.plugins{addonaction}],varargin{:});
        end
end

spm('Pointer','Arrow');
return;

%_______________________________________________________________________
function [Fgraph hErrTexthd] = clearErrorMsg
global opt_info
try
    delete(opt_info.hText.he);
end
Fgraph = spm_figure('GetWin','Graphics');
hErrTexthd = axes('Parent',Fgraph,'Position',[0.02 0.85 0.96 0.02],'Visible','off');

%_______________________________________________________________________
function recoreg1 = LIOMrecoreg(x,y,Q,Sp_rom,Dp_rom,Pvoid,Cid,render_template,m2v,OPidx,OPXYZmm,type)
%codes from nirs_run_coreg_new
% Compute optimal rigid transformation to match fiducial
% positions from one coordinate system to the other
[s R t] = abs_orientation(x,y); % y = s*R(x) + t
estY = zeros(size(y));
for Fi = 1:3 %assume 3 fiducials
    estY(:,Fi) = s*R*x(:,Fi) + t;
end;
% Store coregistration error
err = y - estY;
errVal = sum(err(:).^2);
errValMax = (max(sum(err.^2,1))^0.5)/10;
err %show error on each fiducial
disp(['Error Value for subject: ' num2str(errVal)]);
disp(['Worst coregistration error: ' num2str(errValMax) ' cm']);
% Apply the same transformation to all points in order to
% achieve coregistration
Pp_rom = [Sp_rom Dp_rom];
NP = size(Pp_rom,2);
Ns = size(Sp_rom,2);
Nd = size(Dp_rom,2);
Pp_rmm = zeros(size(Pp_rom));
for Pi = 1:NP
    if ~Pvoid(Pi)
        Pp_rmm(:,Pi) = s*R*Pp_rom(:,Pi) + t; %Optode coordinates, on scalp, unnormalised
    end
end;
% Preserve optode position that have already been re-coregistered
if ~isempty(OPidx)
    if strcmp(type, 'Source')
        for i0 = 1:length(OPidx)
            idx = OPidx(i0);
            Pp_rmm(:,idx) = OPXYZmm(:,i0);
        end
    elseif strcmp(type, 'Detector')
        for i0 = 1:length(OPidx)
            idx = OPidx(i0);
            Pp_rmm(:,Ns+idx) = OPXYZmm(:,i0);
        end
    end
end
% unnormalized -> normalized, for optodes
Pp_wmm = Q * [Pp_rmm;ones(1,NP)];          %% unit : mm %Optode coordinates, on scalp, normalised
% mm -> voxels, for optodes
Pp_rvx = m2v.m\[Pp_rmm;ones(1,NP)];
% PROJECT OPTODE POSITIONS ON CORTEX SURFACE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%render to subject
Nch0 = size(Cid,2)/2;
Pch_rmm = zeros(3,Nch0);
Pch_wmm = zeros(4,Nch0); %different convention of keeping the 4th row
for i=1:Nch0
    %indices of source and detector
    Si = Cid(2,i);
    Di = Cid(3,i)+Ns;
    Pch_rmm(:,i) = (Pp_rmm(:,Si)+Pp_rmm(:,Di))/2;
    Pch_wmm(:,i) = (Pp_wmm(:,Si)+Pp_wmm(:,Di))/2; %Channel coordinates, on scalp
end
%Keep the channel position that have been already coregistered
if strcmp(type, 'Channel')
    if ~isempty(OPidx)
        OPXYZwmm = Q * OPXYZmm;
        for i0 = 1:length(OPidx)
            OPtmp = OPXYZmm(:,i0);
            idx = OPidx(i0);
            Pch_rmm(:,idx) = OPtmp(1:3,:);
            Pch_wmm(:,idx) = OPXYZwmm(:,i0);
            Cwarn = 'Channel positions have been changed without re-placing corresponding sources or detectors. Use source or detector coordinates with caution.';
        end
    end
end
Pch_rvx = m2v.m\[Pch_rmm; ones(1,Nch0)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 2: Cortex projection method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pp_c1_wmm = projection_CS(Pp_wmm);%Using Korean template Optode coordinates, on cortex, normalised
Pp_c1_rmm = zeros(4,NP);
for Pi = 1:NP
    if ~Pvoid(Pi)
        %inversion: normalized -> unnormalized
        Pp_c1_rmm(:,Pi) = Q\Pp_c1_wmm(:,Pi); %Optode coordinates, on cortex, unnormalised
    end
end

coreg_projected_channels_on_cortex_rather_than_midpoint = 1;
if coreg_projected_channels_on_cortex_rather_than_midpoint
    Pch_c1_wmm = projection_CS(Pch_wmm); %Korean's method
    Pch_c1_rmm = Q\Pch_c1_wmm;
else
    for i0 = 1 : Nc
        %indices of source and detector
        Si = Cid(2,i);
        Di = Cid(3,i)+Ns;
        Pch_c1_rmm(:,i) = (Pp_c1_rmm(:,Si)+Pp_c1_rmm(:,Di))/2;
        Pch_c1_wmm(:,i) = (Pp_c1_wmm(:,Si)+Pp_c1_wmm(:,Di))/2; %Channel coordinates, on scalp
    end
end

Pp_c1_rvx = m2v.m\Pp_c1_rmm;
Pch_c1_rvx = m2v.m\Pch_c1_rmm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recoreg1.err.errValofCoreg_mm2 = errVal;
recoreg1.err.errValofCoreg_mm2_all = err;
recoreg1.err.errValofCoreg_cm_worst = errValMax;

recoreg1.corr.F.r.m.mm.p = y;
recoreg1.corr.F.r.o.mm.p = x;
recoreg1.corr.S.N = Ns;
recoreg1.corr.S.r.o.mm.p = Sp_rom;
recoreg1.corr.S.r.m.mm.fp = Pp_rmm(:,1:Ns);
recoreg1.corr.S.r.m.mm.c1.p = Pp_c1_rmm(:,1:Ns);
recoreg1.corr.D.N = Nd;
recoreg1.corr.D.r.o.mm.p = Dp_rom;
recoreg1.corr.D.r.m.mm.fp = Pp_rmm(:,(Ns+1):(Ns+Nd));
recoreg1.corr.D.r.m.mm.c1.p = Pp_c1_rmm(:,(Ns+1):(Ns+Nd));
recoreg1.corr.P.void = Pvoid;
recoreg1.corr.P.r.m.mm.p = Pp_rmm;
recoreg1.corr.P.r.m.mm.fp = Pp_rmm;
recoreg1.corr.P.r.m.mm.c1.p = Pp_c1_rmm(1:3,:);
recoreg1.corr.P.w.m.mm.fp = Pp_wmm;
recoreg1.corr.P.w.m.mm.c1.p = Pp_c1_wmm;
recoreg1.corr.C.N = Nch0*2;
recoreg1.corr.C.id = Cid;
recoreg1.corr.C.r.m.mm.fp = [Pch_rmm; ones(1,size(Pch_rmm, 2))];
recoreg1.corr.C.r.m.mm.c1.p = Pch_c1_rmm;
if exist('Cwarn','var')
    recoreg1.corr.C.warning = Cwarn;
end

recoreg1.corr.S.r.m.vx.fp = Pp_rvx(1:3,1:Ns);
recoreg1.corr.S.r.m.vx.c1.p = Pp_c1_rvx(1:3,1:Ns);
recoreg1.corr.D.r.m.vx.fp = Pp_rvx(1:3,(Ns+1):(Ns+Nd));
recoreg1.corr.D.r.m.vx.c1.p = Pp_c1_rvx(1:3,(Ns+1):(Ns+Nd));
recoreg1.corr.C.r.m.vx.fp = Pch_rvx;
recoreg1.corr.C.r.m.vx.c1.p = Pch_c1_rvx;

recoreg1.Q = Q;
if render_template
    recoreg1.corr.C.w.m.mm.fp = Pch_wmm;
    recoreg1.corr.C.w.m.mm.c1.p = Pch_c1_wmm;
    recoreg1.corr.S.w.m.mm.fp = Pp_wmm(:,1:Ns);
    recoreg1.corr.S.w.m.mm.c1.p = Pp_c1_wmm(:,1:Ns);
    recoreg1.corr.D.w.m.mm.fp = Pp_wmm(:,(Ns+1):(Ns+Nd));
    recoreg1.corr.D.w.m.mm.c1.p = Pp_c1_wmm(:,(Ns+1):(Ns+Nd));
    recoreg1.corr.P.w.m.vx.fp = m2v.w\Pp_wmm;
    recoreg1.corr.P.w.m.vx.c1.p = m2v.w\Pp_c1_wmm;
    recoreg1.corr.C.w.m.vx.fp = m2v.w\Pch_wmm;
    recoreg1.corr.C.w.m.vx.c1.p = m2v.w\Pch_c1_wmm;
    recoreg1.corr.S.w.m.vx.fp = recoreg1.corr.P.w.m.vx.fp(:,1:Ns);
    recoreg1.corr.S.w.m.vx.c1.p = recoreg1.corr.P.w.m.vx.c1.p(:,1:Ns);
    recoreg1.corr.D.w.m.vx.fp = recoreg1.corr.P.w.m.vx.fp(:,(Ns+1):(Ns+Nd));
    recoreg1.corr.D.w.m.vx.c1.p = recoreg1.corr.P.w.m.vx.c1.p(:,(Ns+1):(Ns+Nd));    
end
return;
%_______________________________________________________________________
%_______________________________________________________________________
%function recoreg1 = fill_recoreg(recoreg, render_template, type, new_p)

%_______________________________________________________________________
%_______________________________________________________________________
function coreg = ReplaceChannelCoordinates(coreg,P,type,m2v,render_template)
Cid = coreg.corr.C.id;
Nch0 = coreg.corr.C.N/2;
Ns = coreg.corr.S.N;
Q = coreg.Q;
Cid = Cid(:,1:Nch0);
coreg_projected_channels_on_cortex_rather_than_midpoint = 1;
switch type
    case 'Source'
        Cidx = find(Cid(2,:) == P);
    case 'Detector'
        Cidx = find(Cid(3,:) == P);
    otherwise
end
if render_template %If rendered on template
    for i0 = 1:length(Cidx)
        c0 = Cidx(i0);
        s0 = Cid(2,c0);
        d0 = Cid(3,c0);
        coreg.corr.C.w.m.mm.fp(:,c0) = (coreg.corr.P.w.m.mm.fp(:,s0) + coreg.corr.P.w.m.mm.fp(:,d0+Ns))/2;
        coreg.corr.C.w.m.vx.fp(:,c0) = m2v.w\coreg.corr.C.w.m.mm.fp(:,c0);
        coreg.corr.C.r.m.mm.fp(:,c0) = Q\coreg.corr.C.w.m.mm.fp(:,c0);
        coreg.corr.C.r.m.vx.fp(:,c0) = m2v.m\coreg.corr.C.r.m.mm.fp(:,c0);
        if coreg_projected_channels_on_cortex_rather_than_midpoint
            coreg.corr.C.w.m.mm.c1.p(:,c0) = projection_CS(coreg.corr.C.w.m.mm.fp(:,c0));
            coreg.corr.C.w.m.vx.c1.p(:,c0) = m2v.w\coreg.corr.C.w.m.mm.c1.p(:,c0);
            coreg.corr.C.r.m.mm.c1.p(:,c0) = Q\coreg.corr.C.w.m.mm.c1.p(:,c0);
            coreg.corr.C.r.m.vx.c1.p(:,c0) = m2v.m\coreg.corr.C.r.m.mm.c1.p(:,c0);
        else
            coreg.corr.C.w.m.mm.c1.p(:,c0) = ([coreg.corr.P.w.m.mm.c1.p(:,s0);1] + [coreg.corr.P.w.m.mm.c1.p(:,d0+Ns);1])/2;
            coreg.corr.C.w.m.vx.c1.p(:,c0) = m2v.w\coreg.corr.C.w.m.mm.c1.p(:,c0);
            coreg.corr.C.r.m.mm.c1.p(:,c0) = ([coreg.corr.P.r.m.mm.c1.p(:,s0);1] + [coreg.corr.P.r.m.mm.c1.p(:,d0+Ns);1])/2;
            coreg.corr.C.r.m.vx.c1.p(:,c0) = m2v.m\coreg.corr.C.r.m.mm.c1.p(:,c0);
        end
    end
else %if rendered on patient T1 image
    for i0 = 1:length(Cidx)
        c0 = Cidx(i0);
        s0 = Cid(2,c0);
        d0 = Cid(3,c0);
        coreg.corr.C.r.m.mm.fp(:,c0) = ([coreg.corr.P.r.m.mm.fp(:,s0);1] + [coreg.corr.P.r.m.mm.fp(:,d0+Ns);1])/2;
        coreg.corr.C.r.m.vx.fp(:,c0) = m2v.m\coreg.corr.C.r.m.mm.fp(:,c0);
        if coreg_projected_channels_on_cortex_rather_than_midpoint
            Cw = Q * coreg.corr.C.r.m.mm.fp(:,c0);
            Cw_c1 = projection_CS(Cw); %Korean's method
            coreg.corr.C.r.m.mm.c1.p(:,c0) = Q\Cw_c1;
            coreg.corr.C.r.m.vx.c1.p(:,c0) = m2v.m\coreg.corr.C.r.m.mm.c1.p(:,c0);
        else
            coreg.corr.C.r.m.mm.c1.p(:,c0) = ([coreg.corr.P.r.m.mm.c1.p(:,s0);1] + [coreg.corr.P.r.m.mm.c1.p(:,d0+Ns);1])/2;
            coreg.corr.C.r.m.vx.c1.p(:,c0) = m2v.m\coreg.corr.C.r.m.mm.c1.p(:,c0);
        end
    end
end
return;

%_______________________________________________________________________
%_______________________________________________________________________
function coreg = FillOptodePos(coreg,Q,m2v,P,type,CSlp,CCtx,render_template)
if render_template
    P_wmm = CSlp;
    P_c1_wmm = CCtx;
    P_wmm = [P_wmm;1];
    P_c1_wmm = [P_c1_wmm;1];
    P_wvx = m2v.w \ P_wmm;
    P_c1_wvx = m2v.w \ P_c1_wmm;
    switch type
        case 'Source'
            if isfield(coreg.corr.S, 'w')
                coreg.corr.S.w.m.mm.fp(:,P) = P_wmm;
                coreg.corr.S.w.m.mm.c1.p(:,P) = P_c1_wmm;
                coreg.corr.S.w.m.vx.fp(:,P) = P_wvx;
                coreg.corr.S.w.m.vx.c1.p(:,P) = P_c1_wvx;
            end
            if isfield(coreg.corr.S.r,'m')
                coreg.corr.S.r.m.mm.fp(:,P) = Q\P_wmm;
                coreg.corr.S.r.m.mm.c1.p(:,P) = Q\P_c1_wmm;
                coreg.corr.S.r.m.vx.fp(:,P) = m2v.m\coreg.corr.S.r.m.mm.fp(:,P);
                coreg.corr.S.r.m.vx.c1.p(:,P) = m2v.m\coreg.corr.S.r.m.mm.c1.p(:,P);
            end
            coreg.corr.P.w.m.mm.fp(:,P) = P_wmm;
            coreg.corr.P.w.m.mm.p(:,P) = P_wmm;
            coreg.corr.P.w.m.mm.c1.p(:,P) = P_c1_wmm;
            coreg.corr.P.r.m.mm.fp(:,P) = Q\[P_wmm; 1];
            coreg.corr.P.r.m.mm.c1.p(:,P) = Q\[P_c1_wmm;1];
            coreg = ReplaceChannelCoordinates(coreg,P,type,m2v,1);
        case 'Detector'
            Ns = coreg.corr.S.N;
            if isfield(coreg.corr.D, 'w')
                coreg.corr.D.w.m.mm.fp(:,P) = P_wmm;
                coreg.corr.D.w.m.mm.c1.p(:,P) = P_c1_wmm;
                coreg.corr.D.w.m.vx.fp(:,P) = P_wvx;
                coreg.corr.D.w.m.vx.c1.p(:,P) = P_c1_wvx;
            end
            if isfield(coreg.corr.D.r,'m')
                coreg.corr.D.r.m.mm.fp(:,P) = Q\P_wmm;
                coreg.corr.D.r.m.mm.c1.p(:,P) = Q\P_c1_wmm;
                coreg.corr.D.r.m.vx.fp(:,P) = m2v.m\coreg.corr.S.r.m.mm.fp(:,P);
                coreg.corr.D.r.m.vx.c1.p(:,P) = m2v.m\coreg.corr.S.r.m.mm.c1.p(:,P);
            end
            coreg.corr.P.w.m.mm.fp(:,Ns+P) = P_wmm;
            coreg.corr.P.w.m.mm.p(:,Ms+P) = P_wmm;
            coreg.corr.P.w.m.mm.c1.p(:,Ns+P) = P_c1_wmm;
            coreg.corr.P.r.m.mm.fp(:,Ns+P) = Q\P_wmm;
            coreg.corr.P.r.m.mm.c1.p(:,Ns+P) = Q\P_c1_wmm;
            coreg = ReplaceChannelCoordinates(coreg,P,type,m2v,1);
        case 'Channel'
            coreg.corr.C.w.m.mm.fp(:,P) = P_wmm;
            coreg.corr.C.w.m.mm.c1.p(:,P) = P_c1_wmm;
            coreg.corr.C.w.m.vx.fp(:,P) = P_wvx;
            coreg.corr.C.w.m.vx.c1.p(:,P) = P_c1_wvx;
            coreg.corr.C.r.m.mm.fp(:,P) = Q\P_wmm;
            coreg.corr.C.r.m.mm.c1.p(:,P) = Q\P_c1_wmm;
            coreg.corr.C.r.m.vx.fp(:,P) = m2v.m \ coreg.corr.C.r.m.mm.fp(:,P);
            coreg.corr.C.r.m.vx.c1.p(:,P) = m2v.m \ coreg.corr.C.r.m.mm.c1.p(:,P);
            coreg.corr.C.warning = 'Channel positions have been changed without re-placing corresponding sources or detectors. Use source or detector coordinates with caution.';
        otherwise
            %Should not come to otherwise
    end
else
    P_mm = CSlp;
    P_c1_mm = CCtx;
    P_vx = m2v.m \ [P_mm; 1];
    P_c1_vx = m2v.m \ [P_c1_mm; 1];
    switch type
        case 'Source'
            if isfield(coreg.corr.S.r, 'm')
                coreg.corr.S.r.m.mm.fp(1:3,P) = P_mm;
                coreg.corr.S.r.m.mm.c1.p(1:3,P) = P_c1_mm;
                coreg.corr.S.r.m.vx.fp(1:3,P) = P_vx(1:3,:);
                coreg.corr.S.r.m.vx.c1.p(1:3,P) = P_c1_vx(1:3,:);
            end
            coreg.corr.P.r.m.mm.fp(:,P) = P_mm;
            coreg.corr.P.r.m.mm.p(:,P) = P_mm;
            coreg.corr.P.r.m.mm.c1.p(:,P) = P_c1_mm;
            coreg.corr.P.w.m.mm.fp(:,P) = Q*[P_mm; 1];
            coreg.corr.P.w.m.mm.c1.p(:,P) = Q*[P_c1_mm;1];
            coreg = ReplaceChannelCoordinates(coreg,P,type,m2v,0);
        case 'Detector'
            Ns = coreg.corr.S.N;
            if isfield(coreg.corr.D.r, 'm')
                coreg.corr.D.r.m.mm.fp(1:3,P) = P_mm;
                coreg.corr.D.r.m.mm.c1.p(1:3,P) = P_c1_mm;
                coreg.corr.D.r.m.vx.fp(1:3,P) = P_vx(1:3,:);
                coreg.corr.D.r.m.vx.c1.p(1:3,P) = P_c1_vx(1:3,:);
            end
            coreg.corr.P.r.m.mm.fp(:,Ns+P) = P_mm;
            coreg.corr.P.r.m.mm.p(:,Ns+P) = P_mm;
            coreg.corr.P.r.m.mm.c1.p(:,Ns+P) = P_c1_mm;
            coreg.corr.P.w.m.mm.fp(:,Ns+P) = Q*[P_mm; 1];
            coreg.corr.P.w.m.mm.c1.p(:,Ns+P) = Q*[P_c1_mm;1];
            coreg = ReplaceChannelCoordinates(coreg,P,type,m2v,0);
        case 'Channel'
            coreg.corr.C.r.m.mm.fp(:,P) = [P_mm;1];
            coreg.corr.C.r.m.mm.c1.p(:,P) = [P_c1_mm;1];
            coreg.corr.C.r.m.vx.fp(:,P) = P_vx;
            coreg.corr.C.r.m.vx.c1.p(:,P) = P_c1_vx;
            coreg.corr.C.warning = 'Channel positions have been changed without re-placing corresponding sources or detectors. Use source or detector coordinates with caution.';
        otherwise
            %Should not come to otherwise
    end
end
return;


%_______________________________________________________________________
function setallbutton(on)
global opt_info

if on
    set(opt_info.hText.hb.cr, 'Enable', 'On');
    set(opt_info.hText.hb.bsave, 'Enable', 'On');
    set(opt_info.hText.hb.blpa, 'Enable', 'On');
    set(opt_info.hText.hb.bnas, 'Enable', 'On');
    set(opt_info.hText.hb.brpa, 'Enable', 'On');
    set(opt_info.hText.hb.bnplus, 'Enable', 'On');
    set(opt_info.hText.hb.bnminus, 'Enable', 'On');
else
    set(opt_info.hText.hb.cr, 'Enable', 'Off');
    set(opt_info.hText.hb.bsave, 'Enable', 'Off');
    set(opt_info.hText.hb.blpa, 'Enable', 'Off');
    set(opt_info.hText.hb.bnas, 'Enable', 'Off');
    set(opt_info.hText.hb.brpa, 'Enable', 'Off');
    set(opt_info.hText.hb.bnplus, 'Enable', 'Off');
    set(opt_info.hText.hb.bnminus, 'Enable', 'Off');
end
%_______________________________________________________________________
%_______________________________________________________________________
function addblobs(handle, xyz, t, mat, name)
global st
if nargin < 5
    name = '';
end;
for i=valid_handles(handle),
    if ~isempty(xyz),
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        st.vols{i}.blobs=cell(1,1);
        mx = max([eps max(t)]);
        mn = min([0 min(t)]);
        st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx, 'min',mn,'name',name);
        addcolourbar(handle,1);
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addimage(handle, fname)
global st
for i=valid_handles(handle),
    if isstruct(fname),
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end;
    mat = vol.mat;
    st.vols{i}.blobs=cell(1,1);
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx,'min',mn);
    addcolourbar(handle,1);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredblobs(handle, xyz, t, mat, colour, name)
if nargin < 6
    name = '';
end;
global st
for i=valid_handles(handle),
    if ~isempty(xyz),
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        if ~isfield(st.vols{i},'blobs'),
            st.vols{i}.blobs=cell(1,1);
            bset = 1;
        else
            bset = numel(st.vols{i}.blobs)+1;
        end;
        mx = max([eps maxval(vol)]);
        mn = min([0 minval(vol)]);
        st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
            'max',mx, 'min',mn, ...
            'colour',colour, 'name',name);
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle),
    if isstruct(fname),
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end;
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs'),
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end;
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',colour);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image
global st
for i=valid_handles(handle),
    if isstruct(fname),
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end;
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs'),
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end;
    c = struct('cmap', colourmap,'prop',prop);
    st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx, ...
        'min',mn,'colour',c);
    addcolourbar(handle,bset);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolourbar(vh,bh)
global st
if st.mode == 0,
    axpos = get(st.vols{vh}.ax{2}.ax,'Position');
else
    axpos = get(st.vols{vh}.ax{1}.ax,'Position');
end;
st.vols{vh}.blobs{bh}.cbar = axes('Parent',st.fig,...
    'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
    'Box','on', 'YDir','normal', 'XTickLabel',[], 'XTick',[]);
if isfield(st.vols{vh}.blobs{bh},'name')
    ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmblobs(handle)
global st
for i=valid_handles(handle),
    if isfield(st.vols{i},'blobs'),
        for j=1:numel(st.vols{i}.blobs),
            if isfield(st.vols{i}.blobs{j},'cbar') && ishandle(st.vols{i}.blobs{j}.cbar),
                delete(st.vols{i}.blobs{j}.cbar);
            end;
        end;
        st.vols{i} = rmfield(st.vols{i},'blobs');
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function register(hreg)
global st
tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles(1:24);
if ~isempty(h),
    tmp = st.vols{h(1)}.ax{1}.ax;
    st.registry = struct('hReg',hreg,'hMe', tmp);
    spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'nirs_orthviews');
else
    warning('Nothing to register with');
end;
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairs(arg1)
global st
st.xhairs = 0;
opt = 'on';
if ~strcmp(arg1,'on'),
    opt = 'off';
else
    st.xhairs = 1;
end;
for i=valid_handles(1:24),
    for j=1:3,
        set(st.vols{i}.ax{j}.lx,'Visible',opt);
        set(st.vols{i}.ax{j}.ly,'Visible',opt);
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = pos(arg1)
global st
H = [];
for arg1=valid_handles(arg1),
    is = inv(st.vols{arg1}.premul*st.vols{arg1}.mat);
    H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_reset
global st
if ~isempty(st) && isfield(st,'registry') && ishandle(st.registry.hMe),
    delete(st.registry.hMe); st = rmfield(st,'registry');
end;
my_delete(1:24);
reset_st;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_delete(arg1)
global st
% remove blobs (and colourbars, if any)
rmblobs(arg1);
% remove displayed axes
for i=valid_handles(arg1),
    kids = get(st.fig,'Children');
    for j=1:3,
        if any(kids == st.vols{i}.ax{j}.ax),
            set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
            delete(st.vols{i}.ax{j}.ax);
        end;
    end;
    st.vols{i} = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function resolution(arg1)
global st
if nargin == 0
    res = 1; % Default minimum resolution 1mm
else
    res = arg1;
end
for i=valid_handles(1:24)
    % adapt resolution to smallest voxel size of displayed images
    res = min([res,sqrt(sum((st.vols{i}.mat(1:3,1:3)).^2))]);
end
res      = res/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function move(handle,pos)
global st
for handle = valid_handles(handle),
    st.vols{handle}.area = pos;
end;
bbox;
% redraw(valid_handles(handle));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles(1:24)
    premul = st.Space \ st.vols{i}.premul;
    bb = spm_get_bbox(st.vols{i}, 'fv', premul);
    mx = max([bb ; mx]);
    mn = min([bb ; mn]);
end;
bb = [mn ; mx];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function space(arg1,M,dim)
global st
if ~isempty(st.vols{arg1})
    num = arg1;
    if nargin < 2
        M = st.vols{num}.mat;
        dim = st.vols{num}.dim(1:3);
    end;
    Mat = st.vols{num}.premul(1:3,1:3)*M(1:3,1:3);
    vox = sqrt(sum(Mat.^2));
    if det(Mat(1:3,1:3))<0, vox(1) = -vox(1); end;
    Mat = diag([vox 1]);
    Space = (M)/Mat;
    bb = [1 1 1; dim];
    bb = [bb [1;1]];
    bb=bb*Mat';
    bb=bb(:,1:3);
    bb=sort(bb);
    st.Space  = Space;
    st.bb = bb;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function zoom_op(varargin)
global st
if nargin > 0
    fov = varargin{1};
else
    fov = Inf;
end
if nargin > 1
    res = varargin{2};
else
    res = Inf;
end
if isinf(fov)
    st.bb = maxbb;
elseif isnan(fov) || fov == 0
    current_handle = valid_handles(1:24);
    if numel(current_handle) > 1 % called from check reg context menu
        current_handle = get_current_handle;
    end
    if fov == 0
        % zoom to bounding box of current image ~= 0
        thr = 'nz';
    else
        % zoom to bounding box of current image > chosen threshold
        thr = spm_input('Threshold (Y > ...)', '+1', 'r', '0', 1);
    end
    premul = st.Space \ st.vols{current_handle}.premul;
    st.bb = spm_get_bbox(st.vols{current_handle}, thr, premul);
else
    vx = sqrt(sum(st.Space(1:3,1:3).^2));
    vx = vx.^(-1);
    pos = nirs_orthviews('pos');
    pos = st.Space\[pos ; 1];
    pos = pos(1:3)';
    st.bb = [pos-fov*vx; pos+fov*vx];
end;
resolution(res);
bbox;
redraw_all;
if isfield(st.vols{1},'sdip')
    spm_eeg_inv_vbecd_disp('RedrawDip');
end
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = specify_image(arg1)
global st
H=[];
if isstruct(arg1),
    V = arg1(1);
else
    try
        V = spm_vol(arg1);
    catch
        fprintf('Can not use image "%s"\n', arg1);
        return;
    end;
end;
if numel(V)>1, V=V(1); end

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['nirs_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
    ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
        'YDir','normal','ButtonDownFcn',@repos_start);
    d  = image(0,'Tag','Transverse','Parent',ax,...
        'DeleteFcn',DeleteFcn);
    set(ax,'Ydir','normal','ButtonDownFcn',@repos_start);
    
    lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
    ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
    if ~st.xhairs,
        set(lx,'Visible','off');
        set(ly,'Visible','off');
    end;
    V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
V.mapping   = 'linear';
st.vols{ii} = V;

H = ii;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function repos_start(varargin)
% don't use right mouse button to start reposition
if ~strcmpi(get(gcbf,'SelectionType'),'alt')
    set(gcbf,'windowbuttonmotionfcn',@repos_move, 'windowbuttonupfcn',@repos_end);
    nirs_orthviews('reposition');
end
%_______________________________________________________________________
%_______________________________________________________________________
function repos_move(varargin)
nirs_orthviews('reposition');
%_______________________________________________________________________
%_______________________________________________________________________
function repos_end(varargin)
set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
%_______________________________________________________________________
%_______________________________________________________________________
function addcontexts(handles)
for ii = valid_handles(handles),
    addcontext(ii);
end;
nirs_orthviews('reposition',nirs_orthviews('pos'));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmcontexts(handles)
global st
for ii = valid_handles(handles),
    for i=1:3,
        set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
        st.vols{ii}.ax{i} = rmfield(st.vols{ii}.ax{i},'cm');
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bbox
global st
Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else SD = Dims([2 3])'; end;

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');
sz    = get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles(1:24),
    area = st.vols{i}.area(:);
    area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
    if st.mode == 0,
        sx   = area(3)/(Dims(1)+Dims(3))/1.02;
    else
        sx   = area(3)/(Dims(1)+Dims(2))/1.02;
    end;
    sy   = area(4)/(Dims(2)+Dims(3))/1.02;
    s    = min([sx sy]);
    
    offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
    sky = s*(Dims(2)+Dims(3))*0.02;
    if st.mode == 0,
        offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
        skx = s*(Dims(1)+Dims(3))*0.02;
    else
        offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
        skx = s*(Dims(1)+Dims(2))*0.02;
    end;
    
    % Transverse
    set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
        'Position',[offx offy s*Dims(1) s*Dims(2)],...
        'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Coronal
    set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
        'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
        'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Sagittal
    if st.mode == 0,
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    else
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_all
redraw(1:24);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function mx = maxval(vol)
if isstruct(vol),
    mx = -Inf;
    for i=1:vol.dim(3),
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imx = max(tmp(isfinite(tmp)));
        if ~isempty(imx),mx = max(mx,imx);end
    end;
else
    mx = max(vol(isfinite(vol)));
end;
%_______________________________________________________________________
%_______________________________________________________________________
function mn = minval(vol)
if isstruct(vol),
    mn = Inf;
    for i=1:vol.dim(3),
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imn = min(tmp(isfinite(tmp)));
        if ~isempty(imn),mn = min(mn,imn);end
    end;
else
    mn = min(vol(isfinite(vol)));
end;

%_______________________________________________________________________
%_______________________________________________________________________
function redraw(arg1)
global st
bb   = st.bb;
Dims = round(diff(bb)'+1);
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1),
    M = st.Space\st.vols{i}.premul*st.vols{i}.mat;
    TM0 = [ 1 0 0 -bb(1,1)+1
        0 1 0 -bb(1,2)+1
        0 0 1 -cent(3)
        0 0 0 1];
    TM = inv(TM0*M);
    TD = Dims([1 2]);
    
    CM0 = [ 1 0 0 -bb(1,1)+1
        0 0 1 -bb(1,3)+1
        0 1 0 -cent(2)
        0 0 0 1];
    CM = inv(CM0*M);
    CD = Dims([1 3]);
    
    if st.mode ==0,
        SM0 = [ 0 0 1 -bb(1,3)+1
            0 1 0 -bb(1,2)+1
            1 0 0 -cent(1)
            0 0 0 1];
        SM = inv(SM0*M); 
        SD = Dims([3 2]);
    else
        SM0 = [ 0 -1 0 +bb(2,2)+1
            0  0 1 -bb(1,3)+1
            1  0 0 -cent(1)
            0  0 0 1];
        SM = inv(SM0*M);
        SD = Dims([2 3]);
    end;
    
    try
        imgt  = spm_slice_vol(st.vols{i},TM,TD,st.hld)';
        imgc  = spm_slice_vol(st.vols{i},CM,CD,st.hld)';
        imgs  = spm_slice_vol(st.vols{i},SM,SD,st.hld)';
        ok    = true;
    catch
        fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
        ok     = false;
    end
    if ok,
        % get min/max threshold
        if strcmp(st.vols{i}.window,'auto')
            mn = -Inf;
            mx = Inf;
        else
            mn = min(st.vols{i}.window);
            mx = max(st.vols{i}.window);
        end;
        % threshold images
        imgt = max(imgt,mn); imgt = min(imgt,mx);
        imgc = max(imgc,mn); imgc = min(imgc,mx);
        imgs = max(imgs,mn); imgs = min(imgs,mx);
        % compute intensity mapping, if histeq is available
        if license('test','image_toolbox') == 0
            st.vols{i}.mapping = 'linear';
        end;
        switch st.vols{i}.mapping,
            case 'linear',
            case 'histeq',
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'quadhisteq',
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'loghisteq',
                sw = warning('off','MATLAB:log:logOfZero');
                imgt = log(imgt-min(imgt(:)));
                imgc = log(imgc-min(imgc(:)));
                imgs = log(imgs-min(imgs(:)));
                warning(sw);
                imgt(~isfinite(imgt)) = 0;
                imgc(~isfinite(imgc)) = 0;
                imgs(~isfinite(imgs)) = 0;
                % scale log images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
        end;
        % recompute min/max for display
        if strcmp(st.vols{i}.window,'auto')
            mx = -inf; mn = inf;
        end;
        if ~isempty(imgt),
            tmp = imgt(isfinite(imgt));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end;
        if ~isempty(imgc),
            tmp = imgc(isfinite(imgc));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end;
        if ~isempty(imgs),
            tmp = imgs(isfinite(imgs));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end;
        if mx==mn, mx=mn+eps; end;
        
        if isfield(st.vols{i},'blobs'),
            if ~isfield(st.vols{i}.blobs{1},'colour'),
                % Add blobs for display using the split colourmap
                scal = 64/(mx-mn);
                dcoff = -mn*scal;
                imgt = imgt*scal+dcoff;
                imgc = imgc*scal+dcoff;
                imgs = imgs*scal+dcoff;
                
                if isfield(st.vols{i}.blobs{1},'max'),
                    mx = st.vols{i}.blobs{1}.max;
                else
                    mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
                    st.vols{i}.blobs{1}.max = mx;
                end;
                if isfield(st.vols{i}.blobs{1},'min'),
                    mn = st.vols{i}.blobs{1}.min;
                else
                    mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
                    st.vols{i}.blobs{1}.min = mn;
                end;
                
                vol  = st.vols{i}.blobs{1}.vol;
                M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
                tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                
                %tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
                %tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
                %tmps_z = find(tmps==0);tmps(tmps_z) = NaN;
                
                sc   = 64/(mx-mn);
                off  = 65.51-mn*sc;
                msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
                msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
                msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;
                
                cmap = get(st.fig,'Colormap');
                if size(cmap,1)~=128
                    figure(st.fig)
                    spm_figure('Colormap','gray-hot')
                end;
                redraw_colourbar(i,1,[mn mx],(1:64)'+64);
            elseif isstruct(st.vols{i}.blobs{1}.colour),
                % Add blobs for display using a defined
                % colourmap
                
                % colourmaps
                gryc = (0:63)'*ones(1,3)/63;
                actc = ...
                    st.vols{1}.blobs{1}.colour.cmap;
                actp = ...
                    st.vols{1}.blobs{1}.colour.prop;
                
                % scale grayscale image, not isfinite -> black
                imgt = scaletocmap(imgt,mn,mx,gryc,65);
                imgc = scaletocmap(imgc,mn,mx,gryc,65);
                imgs = scaletocmap(imgs,mn,mx,gryc,65);
                gryc = [gryc; 0 0 0];
                
                % get max for blob image
                if isfield(st.vols{i}.blobs{1},'max'),
                    cmx = st.vols{i}.blobs{1}.max;
                else
                    cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
                end;
                if isfield(st.vols{i}.blobs{1},'min'),
                    cmn = st.vols{i}.blobs{1}.min;
                else
                    cmn = -cmx;
                end;
                
                % get blob data
                vol  = st.vols{i}.blobs{1}.vol;
                M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
                tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                
                % actimg scaled round 0, black NaNs
                topc = size(actc,1)+1;
                tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
                tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
                tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
                actc = [actc; 0 0 0];
                
                % combine gray and blob data to
                % truecolour
                imgt = reshape(actc(tmpt(:),:)*actp+ ...
                    gryc(imgt(:),:)*(1-actp), ...
                    [size(imgt) 3]);
                imgc = reshape(actc(tmpc(:),:)*actp+ ...
                    gryc(imgc(:),:)*(1-actp), ...
                    [size(imgc) 3]);
                imgs = reshape(actc(tmps(:),:)*actp+ ...
                    gryc(imgs(:),:)*(1-actp), ...
                    [size(imgs) 3]);
                csz   = size(st.vols{i}.blobs{1}.colour.cmap);
                cdata = reshape(st.vols{i}.blobs{1}.colour.cmap, [csz(1) 1 csz(2)]);
                redraw_colourbar(i,1,[cmn cmx],cdata);
                
            else
                % Add full colour blobs - several sets at once
                scal  = 1/(mx-mn);
                dcoff = -mn*scal;
                
                wt = zeros(size(imgt));
                wc = zeros(size(imgc));
                ws = zeros(size(imgs));
                
                imgt  = repmat(imgt*scal+dcoff,[1,1,3]);
                imgc  = repmat(imgc*scal+dcoff,[1,1,3]);
                imgs  = repmat(imgs*scal+dcoff,[1,1,3]);
                
                cimgt = zeros(size(imgt));
                cimgc = zeros(size(imgc));
                cimgs = zeros(size(imgs));
                
                colour = zeros(numel(st.vols{i}.blobs),3);
                for j=1:numel(st.vols{i}.blobs) % get colours of all images first
                    if isfield(st.vols{i}.blobs{j},'colour'),
                        colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
                    else
                        colour(j,:) = [1 0 0];
                    end;
                end;
                %colour = colour/max(sum(colour));
                
                for j=1:numel(st.vols{i}.blobs),
                    if isfield(st.vols{i}.blobs{j},'max'),
                        mx = st.vols{i}.blobs{j}.max;
                    else
                        mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
                        st.vols{i}.blobs{j}.max = mx;
                    end;
                    if isfield(st.vols{i}.blobs{j},'min'),
                        mn = st.vols{i}.blobs{j}.min;
                    else
                        mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
                        st.vols{i}.blobs{j}.min = mn;
                    end;
                    
                    vol  = st.vols{i}.blobs{j}.vol;
                    M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                    tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                    tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                    tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                    % check min/max of sampled image
                    % against mn/mx as given in st
                    tmpt(tmpt(:)<mn) = mn;
                    tmpc(tmpc(:)<mn) = mn;
                    tmps(tmps(:)<mn) = mn;
                    tmpt(tmpt(:)>mx) = mx;
                    tmpc(tmpc(:)>mx) = mx;
                    tmps(tmps(:)>mx) = mx;
                    tmpt = (tmpt-mn)/(mx-mn);
                    tmpc = (tmpc-mn)/(mx-mn);
                    tmps = (tmps-mn)/(mx-mn);
                    tmpt(~isfinite(tmpt)) = 0;
                    tmpc(~isfinite(tmpc)) = 0;
                    tmps(~isfinite(tmps)) = 0;
                    
                    cimgt = cimgt + cat(3,tmpt*colour(j,1),tmpt*colour(j,2),tmpt*colour(j,3));
                    cimgc = cimgc + cat(3,tmpc*colour(j,1),tmpc*colour(j,2),tmpc*colour(j,3));
                    cimgs = cimgs + cat(3,tmps*colour(j,1),tmps*colour(j,2),tmps*colour(j,3));
                    
                    wt = wt + tmpt;
                    wc = wc + tmpc;
                    ws = ws + tmps;
                    cdata=permute(shiftdim((1/64:1/64:1)'* ...
                        colour(j,:),-1),[2 1 3]);
                    redraw_colourbar(i,j,[mn mx],cdata);
                end;
                
                imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;
                imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;
                imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;
                
                imgt(imgt<0)=0; imgt(imgt>1)=1;
                imgc(imgc<0)=0; imgc(imgc>1)=1;
                imgs(imgs<0)=0; imgs(imgs>1)=1;
            end;
        else
            scal = 64/(mx-mn);
            dcoff = -mn*scal;
            imgt = imgt*scal+dcoff;
            imgc = imgc*scal+dcoff;
            imgs = imgs*scal+dcoff;
        end;
        
        set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
        set(st.vols{i}.ax{1}.lx,'HitTest','off',...
            'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
        set(st.vols{i}.ax{1}.ly,'HitTest','off',...
            'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
        
        set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
        set(st.vols{i}.ax{2}.lx,'HitTest','off',...
            'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
        set(st.vols{i}.ax{2}.ly,'HitTest','off',...
            'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
        
        set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
        if st.mode ==0,
            set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
            set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
        else
            set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
            set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
        end;
        
        if ~isempty(st.plugins) % process any addons
            for k = 1:numel(st.plugins),
                if isfield(st.vols{i},st.plugins{k}),
                    feval(['spm_ov_', st.plugins{k}], ...
                        'redraw', i, TM0, TD, CM0, CD, SM0, SD);
                end;
            end;
        end;
    end;
end;
drawnow;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_colourbar(vh,bh,interval,cdata)
global st
if isfield(st.vols{vh}.blobs{bh},'cbar')
    if st.mode == 0,
        axpos = get(st.vols{vh}.ax{2}.ax,'Position');
    else
        axpos = get(st.vols{vh}.ax{1}.ax,'Position');
    end;
    % only scale cdata if we have out-of-range truecolour values
    if ndims(cdata)==3 && max(cdata(:))>1
        cdata=cdata./max(cdata(:));
    end;
    image([0 1],interval,cdata,'Parent',st.vols{vh}.blobs{bh}.cbar);
    set(st.vols{vh}.blobs{bh}.cbar, ...
        'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1)...
        (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
        'YDir','normal','XTickLabel',[],'XTick',[]);
    if isfield(st.vols{vh}.blobs{bh},'name')
        ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
    end;
end;
%_______________________________________________________________________
%_______________________________________________________________________
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles(1:24),
    for j=1:3,
        if ~isempty(obj),
            if (st.vols{i}.ax{j}.ax == obj),
                cp = get(obj,'CurrentPoint');
            end;
        end;
        if ~isempty(cp),
            cp   = cp(1,1:2);
            is   = inv(st.Space);
            cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
            switch j,
                case 1,
                    cent([1 2])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,2)-1];
                case 2,
                    cent([1 3])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,3)-1];
                case 3,
                    if st.mode ==0,
                        cent([3 2])=[cp(1)+st.bb(1,3)-1 cp(2)+st.bb(1,2)-1];
                    else
                        cent([2 3])=[st.bb(2,2)+1-cp(1) cp(2)+st.bb(1,3)-1];
                    end;
            end;
            break;
        end;
    end;
    if ~isempty(cent), break; end;
end;
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function handles = valid_handles(handles)
global st;
if isempty(st) || ~isfield(st,'vols')
    handles = [];
else
    handles = handles(:)';
    handles = handles(handles<=24 & handles>=1 & ~rem(handles,1));
    for h=handles,
        if isempty(st.vols{h}), handles(handles==h)=[]; end;
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'hld',1,'fig',fig,'mode',1,'plugins',{{}},'snap',[]);
st.vols = cell(24,1);

xTB = spm('TBs');
if ~isempty(xTB)
    pluginbase = {spm('Dir') xTB.dir};
else
    pluginbase = {spm('Dir')};
end
for k = 1:numel(pluginbase)
    pluginpath = fullfile(pluginbase{k},'nirs_orthviews');
    if isdir(pluginpath)
        pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
        if ~isempty(pluginfiles)
            if ~isdeployed, addpath(pluginpath); end
            % fprintf('nirs_orthviews: Using Plugins in %s\n', pluginpath);
            for l = 1:numel(pluginfiles)
                [p, pluginname, e, v] = spm_fileparts(pluginfiles(l).name);
                st.plugins{end+1} = strrep(pluginname, 'spm_ov_','');
                % fprintf('%s\n',st.plugins{k});
            end;
        end;
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1;end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1;
img(img>cml) = cml;
img(~isfinite(img))  = miscol;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname),
    cmap = evalin('base',acmapname,'[]');
    if isempty(cmap), % not a matrix, is .mat file?
        [p, f, e] = fileparts(acmapname);
        acmat     = fullfile(p, [f '.mat']);
        if exist(acmat, 'file'),
            s    = struct2cell(load(acmat));
            cmap = s{1};
        end;
    end;
end;
if size(cmap, 2)~=3,
    warning('Colormap was not an N by 3 matrix')
    cmap = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function item_parent = addcontext(volhandle)
global st;
%create context menu
fg = spm_figure('Findwin','Graphics');set(0,'CurrentFigure',fg);
%contextmenu
item_parent = uicontextmenu;

%contextsubmenu 0
item00  = uimenu(item_parent, 'Label','unknown image');
nirs_orthviews('context_menu','image_info',item00,volhandle);
item0a    = uimenu(item_parent, 'UserData','pos_mm',     'Callback','nirs_orthviews(''context_menu'',''repos_mm'');','Separator','on');
item0b    = uimenu(item_parent, 'UserData','pos_vx',     'Callback','nirs_orthviews(''context_menu'',''repos_vx'');');
item0c    = uimenu(item_parent, 'UserData','v_value');

%contextsubmenu 1
item1     = uimenu(item_parent,'Label','Zoom','Separator','on');
[zl, rl]  = nirs_orthviews('ZoomMenu');
for cz = numel(zl):-1:1
    if isinf(zl(cz))
        czlabel = 'Full Volume';
    elseif isnan(zl(cz))
        czlabel = 'BBox, this image > ...';
    elseif zl(cz) == 0
        czlabel = 'BBox, this image nonzero';
    else
        czlabel = sprintf('%dx%d mm', 2*zl(cz), 2*zl(cz));
    end
    item1_x = uimenu(item1, 'Label',czlabel,...
        'Callback', sprintf(...
        'nirs_orthviews(''context_menu'',''zoom'',%d,%d)',...
        zl(cz),rl(cz)));
    if isinf(zl(cz)) % default display is Full Volume
        set(item1_x, 'Checked','on');
    end
end

%contextsubmenu 2
checked={'off','off'};
checked{st.xhairs+1} = 'on';
item2     = uimenu(item_parent,'Label','Crosshairs');
item2_1   = uimenu(item2,      'Label','on',  'Callback','nirs_orthviews(''context_menu'',''Xhair'',''on'');','Checked',checked{2});
item2_2   = uimenu(item2,      'Label','off', 'Callback','nirs_orthviews(''context_menu'',''Xhair'',''off'');','Checked',checked{1});

%contextsubmenu 3
if st.Space == eye(4)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end;
item3     = uimenu(item_parent,'Label','Orientation');
item3_1   = uimenu(item3,      'Label','World space', 'Callback','nirs_orthviews(''context_menu'',''orientation'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Voxel space (1st image)', 'Callback','nirs_orthviews(''context_menu'',''orientation'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Voxel space (this image)', 'Callback','nirs_orthviews(''context_menu'',''orientation'',1);','Checked','off');

%contextsubmenu 3
if isempty(st.snap)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end;
item3     = uimenu(item_parent,'Label','Snap to Grid');
item3_1   = uimenu(item3,      'Label','Don''t snap', 'Callback','nirs_orthviews(''context_menu'',''snap'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Snap to 1st image', 'Callback','nirs_orthviews(''context_menu'',''snap'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Snap to this image', 'Callback','nirs_orthviews(''context_menu'',''snap'',1);','Checked','off');

%contextsubmenu 4
if st.hld == 0,
    checked = {'off', 'off', 'on'};
elseif st.hld > 0,
    checked = {'off', 'on', 'off'};
else
    checked = {'on', 'off', 'off'};
end;
item4     = uimenu(item_parent,'Label','Interpolation');
item4_1   = uimenu(item4,      'Label','NN',    'Callback','nirs_orthviews(''context_menu'',''interpolation'',3);', 'Checked',checked{3});
item4_2   = uimenu(item4,      'Label','Bilin', 'Callback','nirs_orthviews(''context_menu'',''interpolation'',2);','Checked',checked{2});
item4_3   = uimenu(item4,      'Label','Sinc',  'Callback','nirs_orthviews(''context_menu'',''interpolation'',1);','Checked',checked{1});

%contextsubmenu 5
% item5     = uimenu(item_parent,'Label','Position', 'Callback','nirs_orthviews(''context_menu'',''position'');');

%contextsubmenu 6
item6       = uimenu(item_parent,'Label','Image','Separator','on');
item6_1     = uimenu(item6,      'Label','Window');
item6_1_1   = uimenu(item6_1,    'Label','local');
item6_1_1_1 = uimenu(item6_1_1,  'Label','auto',       'Callback','nirs_orthviews(''context_menu'',''window'',2);');
item6_1_1_2 = uimenu(item6_1_1,  'Label','manual',     'Callback','nirs_orthviews(''context_menu'',''window'',1);');
item6_1_2   = uimenu(item6_1,    'Label','global');
item6_1_2_1 = uimenu(item6_1_2,  'Label','auto',       'Callback','nirs_orthviews(''context_menu'',''window_gl'',2);');
item6_1_2_2 = uimenu(item6_1_2,  'Label','manual',     'Callback','nirs_orthviews(''context_menu'',''window_gl'',1);');
if license('test','image_toolbox') == 1
    offon = {'off', 'on'};
    checked = offon(strcmp(st.vols{volhandle}.mapping, ...
        {'linear', 'histeq', 'loghisteq', 'quadhisteq'})+1);
    item6_2     = uimenu(item6,      'Label','Intensity mapping');
    item6_2_1   = uimenu(item6_2,    'Label','local');
    item6_2_1_1 = uimenu(item6_2_1,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping'',''linear'');');
    item6_2_1_2 = uimenu(item6_2_1,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping'',''histeq'');');
    item6_2_1_3 = uimenu(item6_2_1,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping'',''loghisteq'');');
    item6_2_1_4 = uimenu(item6_2_1,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping'',''quadhisteq'');');
    item6_2_2   = uimenu(item6_2,    'Label','global');
    item6_2_2_1 = uimenu(item6_2_2,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping_gl'',''linear'');');
    item6_2_2_2 = uimenu(item6_2_2,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping_gl'',''histeq'');');
    item6_2_2_3 = uimenu(item6_2_2,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping_gl'',''loghisteq'');');
    item6_2_2_4 = uimenu(item6_2_2,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','nirs_orthviews(''context_menu'',''mapping_gl'',''quadhisteq'');');
end;
%contextsubmenu 7
item7     = uimenu(item_parent,'Label','Blobs');
item7_1   = uimenu(item7,      'Label','Add blobs');
item7_1_1 = uimenu(item7_1,    'Label','local',  'Callback','nirs_orthviews(''context_menu'',''add_blobs'',2);');
item7_1_2 = uimenu(item7_1,    'Label','global', 'Callback','nirs_orthviews(''context_menu'',''add_blobs'',1);');
item7_2   = uimenu(item7,      'Label','Add image');
item7_2_1 = uimenu(item7_2,    'Label','local',  'Callback','nirs_orthviews(''context_menu'',''add_image'',2);');
item7_2_2 = uimenu(item7_2,    'Label','global', 'Callback','nirs_orthviews(''context_menu'',''add_image'',1);');
item7_3   = uimenu(item7,      'Label','Add colored blobs','Separator','on');
item7_3_1 = uimenu(item7_3,    'Label','local',  'Callback','nirs_orthviews(''context_menu'',''add_c_blobs'',2);');
item7_3_2 = uimenu(item7_3,    'Label','global', 'Callback','nirs_orthviews(''context_menu'',''add_c_blobs'',1);');
item7_4   = uimenu(item7,      'Label','Add colored image');
item7_4_1 = uimenu(item7_4,    'Label','local',  'Callback','nirs_orthviews(''context_menu'',''add_c_image'',2);');
item7_4_2 = uimenu(item7_4,    'Label','global', 'Callback','nirs_orthviews(''context_menu'',''add_c_image'',1);');
item7_5   = uimenu(item7,      'Label','Remove blobs',        'Visible','off','Separator','on');
item7_6   = uimenu(item7,      'Label','Remove colored blobs','Visible','off');
item7_6_1 = uimenu(item7_6,    'Label','local', 'Visible','on');
item7_6_2 = uimenu(item7_6,    'Label','global','Visible','on');

for i=1:3,
    set(st.vols{volhandle}.ax{i}.ax,'UIcontextmenu',item_parent);
    st.vols{volhandle}.ax{i}.cm = item_parent;
end;

% process any plugins
for k = 1:numel(st.plugins),
    feval(['spm_ov_', st.plugins{k}], ...
        'context_menu', volhandle, item_parent);
    if k==1
        h = get(item_parent,'Children');
        set(h(1),'Separator','on'); 
    end
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function c_menu(varargin)
global st

switch lower(varargin{1}),
    case 'image_info',
        if nargin <3,
            current_handle = get_current_handle;
        else
            current_handle = varargin{3};
        end;
        if isfield(st.vols{current_handle},'fname'),
            [p,n,e,v] = spm_fileparts(st.vols{current_handle}.fname);
            if isfield(st.vols{current_handle},'n')
                v = sprintf(',%d',st.vols{current_handle}.n);
            end;
            set(varargin{2}, 'Label',[n e v]);
        end;
        delete(get(varargin{2},'children'));
        if exist('p','var')
            item1 = uimenu(varargin{2}, 'Label', p);
        end;
        if isfield(st.vols{current_handle},'descrip'),
            item2 = uimenu(varargin{2}, 'Label',...
                st.vols{current_handle}.descrip);
        end;
        dt = st.vols{current_handle}.dt(1);
        item3 = uimenu(varargin{2}, 'Label', sprintf('Data type: %s', spm_type(dt)));
        str   = 'Intensity: varied';
        if size(st.vols{current_handle}.pinfo,2) == 1,
            if st.vols{current_handle}.pinfo(2),
                str = sprintf('Intensity: Y = %g X + %g',...
                    st.vols{current_handle}.pinfo(1:2)');
            else
                str = sprintf('Intensity: Y = %g X', st.vols{current_handle}.pinfo(1)');
            end;
        end;
        item4  = uimenu(varargin{2}, 'Label',str);
        item5  = uimenu(varargin{2}, 'Label', 'Image dims', 'Separator','on');
        item51 = uimenu(varargin{2}, 'Label',...
            sprintf('%dx%dx%d', st.vols{current_handle}.dim(1:3)));
        prms   = spm_imatrix(st.vols{current_handle}.mat);
        item6  = uimenu(varargin{2}, 'Label','Voxel size', 'Separator','on');
        item61 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', prms(7:9)));
        item7  = uimenu(varargin{2}, 'Label','Origin', 'Separator','on');
        item71 = uimenu(varargin{2}, 'Label',...
            sprintf('%.2f %.2f %.2f', prms(1:3)));
        R      = spm_matrix([0 0 0 prms(4:6)]);
        item8  = uimenu(varargin{2}, 'Label','Rotations', 'Separator','on');
        item81 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(1,1:3)));
        item82 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(2,1:3)));
        item83 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(3,1:3)));
        item9  = uimenu(varargin{2},...
            'Label','Specify other image...',...
            'Callback','nirs_orthviews(''context_menu'',''swap_img'');',...
            'Separator','on');
        
    case 'repos_mm',
        oldpos_mm = nirs_orthviews('pos');
        newpos_mm = spm_input('New Position (mm)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_mm),3);
        nirs_orthviews('reposition',newpos_mm);
        
    case 'repos_vx'
        current_handle = get_current_handle;
        oldpos_vx = nirs_orthviews('pos', current_handle);
        newpos_vx = spm_input('New Position (voxels)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_vx),3);
        newpos_mm = st.vols{current_handle}.mat*[newpos_vx;1];
        nirs_orthviews('reposition',newpos_mm(1:3));
        
    case 'zoom'
        zoom_all(varargin{2:end});
        bbox;
        redraw_all;
        
    case 'xhair',
        nirs_orthviews('Xhairs',varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles),
            z_handle = get(findobj(cm_handles(i),'label','Crosshairs'),'Children');
            set(z_handle,'Checked','off'); %reset check
            if strcmp(varargin{2},'off'), op = 1; else op = 2; end
            set(z_handle(op),'Checked','on');
        end;
        
    case 'orientation',
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles),
            z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
            set(z_handle,'Checked','off');
        end;
        if varargin{2} == 3,
            nirs_orthviews('Space');
            for i = 1:numel(cm_handles),
                z_handle = findobj(cm_handles(i),'label','World space');
                set(z_handle,'Checked','on');
            end;
        elseif varargin{2} == 2,
            nirs_orthviews('Space',1);
            for i = 1:numel(cm_handles),
                z_handle = findobj(cm_handles(i),'label',...
                    'Voxel space (1st image)');
                set(z_handle,'Checked','on');
            end;
        else
            nirs_orthviews('Space',get_current_handle);
            z_handle = findobj(st.vols{get_current_handle}.ax{1}.cm, ...
                'label','Voxel space (this image)');
            set(z_handle,'Checked','on');
            return;
        end;
        
    case 'snap',
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles),
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle,'Checked','off');
        end;
        if varargin{2} == 3,
            st.snap = [];
        elseif varargin{2} == 2,
            st.snap = 1;
        else
            st.snap = get_current_handle;
            z_handle = get(findobj(st.vols{get_current_handle}.ax{1}.cm,'label','Snap to Grid'),'Children');
            set(z_handle(1),'Checked','on');
            return;
        end;
        for i = 1:numel(cm_handles),
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle(varargin{2}),'Checked','on');
        end;
        
    case 'interpolation',
        tmp        = [-4 1 0];
        st.hld     = tmp(varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles),
            z_handle = get(findobj(cm_handles(i),'label','Interpolation'),'Children');
            set(z_handle,'Checked','off');
            set(z_handle(varargin{2}),'Checked','on');
        end;
        redraw_all;
        
    case 'window',
        current_handle = get_current_handle;
        if varargin{2} == 2,
            nirs_orthviews('window',current_handle);
        else
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%.2f %.2f', st.vols{current_handle}.window);
            else
                defstr = '';
            end;
            [w yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            nirs_orthviews('window',current_handle,w);
        end;
        
    case 'window_gl',
        if varargin{2} == 2,
            for i = 1:numel(get_cm_handles),
                st.vols{i}.window = 'auto';
            end;
        else
            current_handle = get_current_handle;
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%d %d', st.vols{current_handle}.window);
            else
                defstr = '';
            end;
            [w yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            for i = 1:numel(get_cm_handles),
                st.vols{i}.window = w;
            end;
        end;
        redraw_all;
        
    case 'mapping',
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', ...
            'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        current_handle = get_current_handle;
        cm_handles = get_cm_handles;
        st.vols{current_handle}.mapping = varargin{2};
        z_handle = get(findobj(cm_handles(current_handle), ...
            'label','Intensity mapping'),'Children');
        for k = 1:numel(z_handle)
            c_handle = get(z_handle(k), 'Children');
            set(c_handle, 'checked', 'off');
            set(c_handle(checked), 'checked', 'on');
        end;
        redraw_all;
        
    case 'mapping_gl',
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        cm_handles = get_cm_handles;
        for k = valid_handles(1:24),
            st.vols{k}.mapping = varargin{2};
            z_handle = get(findobj(cm_handles(k), ...
                'label','Intensity mapping'),'Children');
            for l = 1:numel(z_handle)
                c_handle = get(z_handle(l), 'Children');
                set(c_handle, 'checked', 'off');
                set(c_handle(checked), 'checked', 'on');
            end;
        end;
        redraw_all;
        
    case 'swap_img',
        current_handle = get_current_handle;
        newimg = spm_select(1,'image','select new image');
        if ~isempty(newimg)
            new_info = spm_vol(newimg);
            fn = fieldnames(new_info);
            for k=1:numel(fn)
                st.vols{current_handle}.(fn{k}) = new_info.(fn{k});
            end;
            nirs_orthviews('context_menu','image_info',get(gcbo, 'parent'));
            redraw_all;
        end
        
    case 'add_blobs',
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle; end;
        spm_figure('Clear','Interactive');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            for i = 1:numel(cm_handles),
                addblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M);
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','nirs_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1,
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','nirs_orthviews(''context_menu'',''remove_blobs'',1);');
                end;
            end;
            redraw_all;
        end
        
    case 'remove_blobs',
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle; end;
        for i = 1:numel(cm_handles),
            rmblobs(cm_handles(i));
            c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
            delete(get(c_handle,'Children'));
            set(c_handle,'Visible','off');
        end;
        redraw_all;
        
    case 'add_image',
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle; end;
        spm_figure('Clear','Interactive');
        fname = spm_select(1,'image','select image');
        if ~isempty(fname)
            for i = 1:numel(cm_handles),
                addimage(cm_handles(i),fname);
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','nirs_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1,
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','nirs_orthviews(''context_menu'',''remove_blobs'',1);');
                end;
            end;
            redraw_all;
        end
        
    case 'add_c_blobs',
        % Add blobs to the image - in full colour
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle; end;
        spm_figure('Clear','Interactive');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            c         = spm_input('Colour','+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours   = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names   = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',xSPM.title,c_names{c});
            for i = 1:numel(cm_handles),
                addcolouredblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M,colours(c,:),xSPM.title);
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1   = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');nirs_orthviews(''context_menu'',''remove_c_blobs'',2,c);',...
                    'UserData',c);
                if varargin{2} == 1,
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');nirs_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end;
            end;
            redraw_all;
        end
    case 'remove_c_blobs',
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle; end;
        colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
        for i = 1:numel(cm_handles),
            if isfield(st.vols{cm_handles(i)},'blobs'),
                for j = 1:numel(st.vols{cm_handles(i)}.blobs),
                    if all(st.vols{cm_handles(i)}.blobs{j}.colour == colours(varargin{3},:));
                        if isfield(st.vols{cm_handles(i)}.blobs{j},'cbar')
                            delete(st.vols{cm_handles(i)}.blobs{j}.cbar);
                        end
                        st.vols{cm_handles(i)}.blobs(j) = [];
                        break;
                    end;
                end;
                rm_c_menu = findobj(st.vols{cm_handles(i)}.ax{1}.cm,'Label','Remove colored blobs');
                delete(gcbo);
                if isempty(st.vols{cm_handles(i)}.blobs),
                    st.vols{cm_handles(i)} = rmfield(st.vols{cm_handles(i)},'blobs');
                    set(rm_c_menu, 'Visible', 'off');
                end;
            end;
        end;
        redraw_all;
        
    case 'add_c_image',
        % Add truecolored image
        cm_handles = valid_handles(1:24);
        if varargin{2} == 2, cm_handles = get_current_handle;end;
        spm_figure('Clear','Interactive');
        fname   = spm_select([1 Inf],'image','select image(s)');
        for k = 1:size(fname,1)
            c       = spm_input(sprintf('Image %d: Colour',k),'+1','m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',fname(k,:),c_names{c});
            for i = 1:numel(cm_handles),
                addcolouredimage(cm_handles(i),fname(k,:),colours(c,:));
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1 = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');nirs_orthviews(''context_menu'',''remove_c_blobs'',2,c);','UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');nirs_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
end;
%_______________________________________________________________________
%_______________________________________________________________________
function current_handle = get_current_handle
cm_handle      = get(gca,'UIContextMenu');
cm_handles     = get_cm_handles;
current_handle = find(cm_handles==cm_handle);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cm_pos
global st
for i = 1:numel(valid_handles(1:24)),
    if isfield(st.vols{i}.ax{1},'cm')
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_mm'),...
            'Label',sprintf('mm:  %.1f %.1f %.1f',nirs_orthviews('pos')));
        pos = nirs_orthviews('pos',i);
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_vx'),...
            'Label',sprintf('vx:  %.1f %.1f %.1f',pos));
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','v_value'),...
            'Label',sprintf('Y = %g',spm_sample_vol(st.vols{i},pos(1),pos(2),pos(3),st.hld)));
    end
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cm_handles = get_cm_handles
global st
cm_handles = [];
for i=valid_handles(1:24),
    cm_handles = [cm_handles st.vols{i}.ax{1}.cm];
end
return;
%_______________________________________________________________________
%_______________________________________________________________________
function zoom_all(zoom,res)
global st
cm_handles = get_cm_handles;
zoom_op(zoom,res);
for i = 1:numel(cm_handles)
    z_handle = get(findobj(cm_handles(i),'label','Zoom'),'Children');
    set(z_handle,'Checked','off');
    if isinf(zoom)
        set(findobj(z_handle,'Label','Full Volume'),'Checked','on');
    elseif zoom > 0
        set(findobj(z_handle,'Label',sprintf('%dx%d mm', 2*zoom, 2*zoom)),'Checked','on');
    end % leave all unchecked if either bounding box option was chosen
end
return;
