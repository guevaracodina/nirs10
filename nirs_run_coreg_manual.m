function out = nirs_run_coreg_manual(job)
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

%Usage: using previously created NIRS.mat containing information of 
%optode and fiducial positions in one system of coordinates, this function
%will

%Load NIRS.mat information
try
    %standalone coregistration
    WanatT1 = job.Coreg_choice.Coreg_standalone.WanatT1{1,1};
    NormParamsFile = job.Coreg_choice.Coreg_standalone.NormParams{1,1};
    NormParams = load(NormParamsFile);
    Vw = spm_vol(WanatT1);
    volume = spm_read_vols(Vw);
    V = spm_vol(job.Vsegmented{1,1});
    standalone = 1;
catch
    %coregistration based on NIRS.mat
    load(job.Coreg_choice.CoregFromNIRS.NIRSmatSingle{1,1});
    V = spm_vol(job.Vsegmented{1,1});
    volume = spm_read_vols(V);
    standalone = 0;
end

try 
    subjectID = NIRS.subjectID;
catch
    subjectID = '';
end

try 
    initSrcPos = NIRS.initSrcPos;
    initDetPos = NIRS.initDetPos;
catch
    
end

% CLICKOPTODESPOS  Allow the user to manually select the positions of
%                 optodes relative to the volume previous to Monte-Carlo
%                 simulations. This GUI is called by the main GUI
%                 "prepMCsim". 
%
%       Comments displayed at the command line in response 
%       to the help command. 
% (Leave a blank line following the help.)


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Initialization tasks (before components creation)  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Data initializations
    
    % Input arguments :
    % volume : 3D matrix of head MRI volume
    % subjectID : name used for window title
    % initSrc/DetPos : initial opotdes positions (those already entered
    % before, the newly entered ones will be added to the list)
    if ~exist('initSrcPos','var') || isempty(initSrcPos)
        initSrcPos = zeros(0,3);
    end
    if ~exist('initDetPos','var') || isempty(initDetPos)
        initDetPos = zeros(0,3);
    end

    % Output arguments : will contain user-defined optodes positions, or an
    % empty matrix if no position is entered or figure is closed
    outputPos = {};    
    
    % Optodes positions (to be updated to match with user-defined positions
    % (from clicks or text editing)
    optodes.srcPos = initSrcPos;
    optodes.detPos = initDetPos;
    
    % Volume min and max values (for color mapping for display)
    clims = [min(volume(:)) max(volume(:))];


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Construct the components  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  Create and then hide the GUI as it is being constructed.
   hfig = figure('Visible','off','Position',[100 200 1100 700],...
          'Color',[.925 .914 .847],...
          'DeleteFcn',{@hfig_DeleteFcn});

% --- Project title and instructions for the user
   htitle = uicontrol('Style','text','String',char(subjectID),...
       'Position',[0 670 250 30],'FontSize',14,...
       'Min',0,'Max',2,...
       'BackgroundColor',[0 0 0],'ForegroundColor',get(hfig,'Color'));
   hinstructions = uicontrol('Style','text',...
       'Position',[50 270 1000 20],'FontSize',12,...
       'Min',0,'Max',2,...
       'BackgroundColor',[.8 .8 .75],'ForegroundColor',[0 0 0]);
   set(hinstructions,'String',['Repérez les optodes',...
       ' en cliquant sur un des trois axes,',...
       ' puis à la position de la source (clic gauche)',...
       ' ou du détecteur (clic droit).']);

% --- Axes and sliders for volume display and graphical inputs
   haxes_1axial = axes('Parent',hfig,'NextPlot','replacechildren',...
                  'Units','pixels','Position',[50 370 300 250],...
                  'Tag','axes_1axial',...
                  'ButtonDownFcn',{@axes_ButtonDownFcn});
                  title({'Bas -> haut'; '  -  Axiale  -  '},...
                      'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   haxes_2coronal = axes('Parent',hfig,'NextPlot','replacechildren',...
                    'Units','pixels','Position',[400 370 300 250],...
                    'Tag','axes_2coronal',...
                    'ButtonDownFcn',{@axes_ButtonDownFcn});
                    title({'Derrière -> Devant'; ' - Coronale - '},...
                     'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   haxes_3sagittal = axes('Parent',hfig,'NextPlot','replacechildren',...
                     'Units','pixels','Position',[750 370 300 250],...
                     'Tag','axes_3sagittal',...
                     'ButtonDownFcn',{@axes_ButtonDownFcn});
                     title({'Gauche -> droite'; ' - Sagittale - '},...
                      'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   set([haxes_1axial haxes_2coronal haxes_3sagittal],...
       'DataAspectRatio',[1 1 1],...
       'PlotBoxAspectRatio',[1 1 1]);
   
   hslider_1axial = uicontrol(hfig,'Style','slider',...
                    'Position',[50 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_1axial_Callback});
   hslider_2coronal = uicontrol(hfig,'Style','slider',...
                    'Position',[400 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_2coronal_Callback});
   hslider_3sagittal = uicontrol(hfig,'Style','slider',...
                    'Position',[750 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_3sagittal_Callback});  
                
   % Volume dimension for initializing sliders
   dim = size(volume);
          set(hslider_1axial,'Enable','on',...
              'Min',1,'Max',dim(3),...
              'SliderStep',[1/dim(3) 15/dim(3)]);
          set(hslider_2coronal,'Enable','on',...
              'Min',1,'Max',dim(2),...
              'SliderStep',[1/dim(2) 15/dim(2)]);
          set(hslider_3sagittal,'Enable','on',...
              'Min',1,'Max',dim(1),...
              'SliderStep',[1/dim(1) 15/dim(1)]);

% --- Edit texts for text inputs
   hpanel_edits = uipanel('Parent',hfig,'Title','Positions entrées',...
                  'BackgroundColor',[.8 .8 .75],...
                  'Units','pixels','Position',[50 25 500 225]);   
     hedit_srcPos = uicontrol(hpanel_edits,'Style','edit',...
                       'Units','normalized','Position',[.1 .05 .35 .9],...
                       'Min',0,'Max',10,'String',num2str(optodes.srcPos),...
                       'Tag','s',...
                       'Callback',{@edit_Pos_Callback});
     hedit_detPos = uicontrol(hpanel_edits,'Style','edit',...
                       'Units','normalized','Position',[.55 .05 .35 .9],...
                       'Min',0,'Max',10,'String',num2str(optodes.detPos),...
                       'Tag','d',...
                       'Callback',{@edit_Pos_Callback});


% --- Pushbuttons for user interactions
   hpushbutton_OK = uicontrol(hfig,'Style','pushbutton',...
                  'String','Terminé',...
                  'Units','pixels','Position',[800 110 150 130],...
                  'Callback',{@pushbutton_OK_Callback});
   hpushbutton_cancel = uicontrol(hfig,'Style','pushbutton',...
                  'String','Annuler',...
                  'Units','pixels','Position',[800 25 150 50],...
                  'Callback',{@pushbutton_cancel_Callback});

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Initialization tasks (after components creation)  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% --- GUI initializations
   set([hfig htitle haxes_1axial...
       haxes_2coronal haxes_3sagittal hslider_1axial...
       hslider_2coronal hslider_3sagittal],...
   'Units','normalized'); % So components resize automatically

   % General figure appearance 
   set(hfig,'Name','Repérage des optodes',... % Window title
       'NumberTitle','off'); % Do not display figure number
   % Move the GUI to the center of the screen, then up a little
   movegui(hfig,'center');
   set(hfig,'Position',get(hfig,'Position')+[0 .1 0 0]);
   
   % Initialize each plot and invert vertical axes
   % (see explanation below in "draw" function
    updateDisplay(1:3);
    initializeAxesTicks;
    
   
   % Make the GUI visible.
   set(hfig,'Visible','on');
   
   % Block execution until user has interacted with or closed the figure
    uiwait(hfig);

    % Return the coordinates selected by the user if it is requested
    outputPos{1} = optodes.srcPos;
    outputPos{2} = optodes.detPos;
    if nargout>0
        [varargout{1:nargout}] = outputPos{:};
    end
    
    if standalone
        %outputPos{1} = NIRS.SrcPos_MNIv;
        %outputPos{2} = NIRS.DetPos_MNIv;
        NIRS.SrcPos_MNIv = outputPos{1};
        NIRS.DetPos_MNIv = outputPos{2};
        temp_SrcPos_MNIv = [outputPos{1} ones(size(outputPos{1},1),1)];
        tempSrc = (Vw.mat * temp_SrcPos_MNIv')';       
        NIRS.SrcPos_MNIw = tempSrc(:,1:3); 
        temp_DetPos_MNIv = [outputPos{2} ones(size(outputPos{2},1),1)];
        tempDet = (Vw.mat * temp_DetPos_MNIv')';
        NIRS.DetPos_MNIw = tempDet(:,1:3); 
        Affine = NormParams.Affine;
        %Naive!!! RIGID transformation is not good enough
%         tempSrc = (((NormParams.VG.mat/Affine)/NormParams.VF.mat) * tempSrc')';
%         NIRS.SrcPos = tempSrc(:,1:3);
%         tempDet = (((NormParams.VG.mat/Affine)/NormParams.VF.mat) * tempDet')';
%         NIRS.DetPos = tempDet(:,1:3);
        NIRS.SrcPos_MNI = spm_get_orig_coord(NIRS.SrcPos_MNIw,NormParamsFile);
        NIRS.DetPos_MNI = spm_get_orig_coord(NIRS.DetPos_MNIw,NormParamsFile);
    
        %
        
        try
            out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
            save(job.NIRSmat{1,1},'NIRS');
        catch
            out.NIRSmat{1} = 'NIRS.mat';
            save('NIRS','NIRS');
        end
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Callbacks for CLICKOPTODESPOS  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 

 
% --- GRAPHICAL INPUTS

   % --- Extracting optode coordinates upon user click on a location
    function axes_ButtonDownFcn(hObject,eventdata)
        % When the user double-clicks, the first click will call this
        % ButtonDown function and the second one will respond to the
        % "ginput" command:
        [horizClicked vertClicked SorD] = ginput(1);
        %get(hObject,'ButtonDownFcn')
        % The newly selected optode will be added to the list
        updateClick(hObject,round(horizClicked),round(vertClicked),SorD);
%         % Re-define the ButtonDown function which has been reset (why?)
%         set(hObject,'ButtonDownFcn',{@axes_ButtonDownFcn})
    end

    % --- Slider callbacks
     % --- Navigate thourgh 3D current volume with sliders
   function slider_1axial_Callback(hObject,eventdata)
       % Update display on 1st plot
       updateDisplay(1);
   end

   function slider_2coronal_Callback(hObject,eventdata)
       % Update display on second plot
       updateDisplay(2);
   end

   function slider_3sagittal_Callback(hObject,eventdata)
       % Update display on 3rd plot
       updateDisplay(3);
   end


% --- EDIT TEXT INPUTS
      % Edit texts are updated whenever a graphical input is given
      % Graphics are updated to match new entered position

    function edit_Pos_Callback(hObject,eventdata)
        SorD = get(hObject,'Tag'); % source or detector edit text
        updateEdit(hObject,SorD);
    end


% --- FINISHING

   % --- Terminate selections
    function pushbutton_OK_Callback(hObject,eventdata)
        % Values in optodes.srcPos and optodes.detPos will be output when
        % execution resumes
       uiresume;
        % End
        delete(hfig);
    end

     function pushbutton_cancel_Callback(hObject,eventdata)
        % Cancel : return to initial values (possibly empty)
        optodes.srcPos = initSrcPos;
        optodes.detPos = initDetPos;
        uiresume;
        % End
        delete(hfig);
      end

    % --- If figure is closed prematurely
    function hfig_DeleteFcn(hObject,eventdata)
        % Return the positions anyway
        uiresume;
    end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Utility functions for PREPMCSIM2  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % --- Initialize vertical axes ticks (inverted because of imrotate)
   % (see exmplanation in "draw" function below)
   % Must have initialized display with updateDisplay prior to this!
    function initializeAxesTicks
        axesHandles = [haxes_1axial; haxes_2coronal; haxes_3sagittal];
        for view = 1:3 % for each axes
            hax = axesHandles(view);
            if ~isempty(get(hax,'Children')) % make sure display was initialized
                ticksYdef = get(hax,'YTick'); % get current ticks (ex: 0-256; i=5,255)
                %ticksYnew = size(imrotate(image2D,90),1) - ticksYdef;
                ticksYnew = max(get(hax,'ylim')) - ticksYdef; % calculate new ticks : y -> -y + max(y) (invert) (i=251,1)
                set(hax,'YTick',sort(ticksYnew,'ascend')) % set new ticks (at i=251,1)
                set(hax,'YTickLabel',num2str(sort(ticksYdef','descend'))) % but label them (label i=5,255)
            end
        end
    end

   
   % --- Update display to match current matrix choice and display options
    function updateDisplay(viewNumbers)
        % viewNumbers indicate which plot(s) to update : numbers 1, 2 or 3
        
        % Display each view (axial, coronal, sagittal) in the corresponding axes
        sliderHandles = [hslider_1axial; hslider_2coronal; hslider_3sagittal];
        axesHandles = [haxes_1axial; haxes_2coronal; haxes_3sagittal];
        if ~isempty(volume)
            for view = viewNumbers
                % Slice number is indicated by slider
                sliceNo = ceil(get(sliderHandles(view),'Value'));
                switch view
                    case 1
                        theSlice = squeeze(volume(:,:,sliceNo));
                        coordHoriz = 1;
                        coordVert = 2;
                    case 2
                        theSlice = squeeze(volume(:,sliceNo,:));
                        coordHoriz = 1;
                        coordVert = 3;
                    case 3
                        theSlice = squeeze(volume(sliceNo,:,:));
                        coordHoriz = 2;
                        coordVert = 3;
                    otherwise
                        % Plots are not updated
                end
                
                % Clear axes of previous children (image, rectangle...)
                cla(axesHandles(view));
                
                % Draw the image of the slice on the plot, indicate slice
                % number on top and label colorbar on the left
                draw(theSlice,axesHandles(view),clims,...
                    sliceNo);

                if any(optodes.srcPos) % inialized to [0 0 0]
                % Display sources and detectors which are in this plane
                    for nsrc = 1:size(optodes.srcPos,1); % For each source of known position
                        if optodes.srcPos(nsrc,4-view) == sliceNo;
                        % 4 - view : axial/z , coronal/y, sagittal/x
                             axes(axesHandles(view))
                             hold on, plot(optodes.srcPos(nsrc,coordHoriz),...
                                 size(imrotate(theSlice,90),1)-optodes.srcPos(nsrc,coordVert),... % vertical axis is inverted
                                 '*y','MarkerSize',10);
                        end
                    end
                end
                if any(optodes.detPos)
                    for ndet = 1:size(optodes.detPos,1); % For each source of known position
                        if optodes.detPos(ndet,4-view) == sliceNo;
                        % 4-view : axial/z , coronal/y, sagittal/x
                             axes(axesHandles(view))
                             hold on, plot(optodes.detPos(ndet,coordHoriz),...
                                 size(imrotate(theSlice,90),1)-optodes.detPos(ndet,coordVert),... % vertical axis is inverted
                                 'og','MarkerSize',10);
                        end
                    end
                end
                
            end
            
        else % clear (but do not reset) axes
            axes(haxes_1axial); cla;
            axes(haxes_2coronal); cla;
            axes(haxes_3sagittal); cla;
        end
        
    end
   

    % --- Draw an image on specified axes
    function draw(image2D,h_axes,clims,sliceNumber,ticks)

    % inputs :
    % plan : 2D mage to display (could be a slice of a 3D matrix)
    % h_axes : handles to (existing) axes where the image will be drawn
    % optional inputs :
    % clims : [min max] limits for color mapping (see imagesc)
    % sliceNumber : will be displayed on top of axes (ID the slice within
    % the 3D volume)
    % ticks : labels for the colorbar in a cell array of strings
    % ex: {'Air'; 'Mgr'; 'Mbl'; 'CSF'; 'Skull'; 'Scalp'}

    % Remove singleton dimensions
    if size(size(image2D))>2
        image2D = squeeze(image2D);
    end
    
    % Default values for arguments
        switch nargin
            case 2
                clims = [min(image2D(:)) max(image2D(:))];
                sliceNumber = [];
                ticks = [];
            case 3
                sliceNumber = 0;
                ticks = [];
            case 4
                ticks = [];
            case 5
            otherwise
                return
        end

       
        axes(h_axes)
        imagesc(imrotate(image2D,90),clims);
        axis ij; axis tight; colormap gray;
        % Rotate because imagesc plots the first coordinate vertical
        % and the second horizontal. However this rotation causes the
        % vertical axis to be inverted relative to the convention, so reset
        % the graduations so that they reflect the real axis direction.     
       
     
        % Colorbar labelled or not
        if ~isempty(ticks)
            colorbar('Location','EastOutside','Ytick',0:length(ticks)-1,'YTickLabel',ticks)
            % maximum 5 tissues + air
        else
            colorbar('Location','EastOutside')
        end

        % Display slice number within the 3D matrix
        if sliceNumber~=0
            txt = text(0.01,1.04,'','Units','normalized');
            set(txt,'String',strcat('Tranche #',int2str(sliceNumber)),...
                'BackgroundColor',[.8 .8 .75],'Color','k','FontSize',10,'FontWeight','Normal');
        end

        
    end
        

    % --- Update the list and display of optodes positions to include the
    % newly clicked optode (source or detector)
    function updateClick(haxes,horizPos,vertPos,SorD)
        % Optode type : source (1) or detector (3) (if 2, do nothing)
        
        % Coordinates of clicked point
        switch get(haxes,'Tag')
            case 'axes_1axial'
                % Position sélectionnée :
                x = round(horizPos);
                y = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                z = ceil(get(hslider_1axial,'Value'));

            case 'axes_2coronal'
                % Position sélectionnée :
                x = round(horizPos);
                z = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                y = ceil(get(hslider_2coronal,'Value'));
               
            case 'axes_3sagittal'
                % Position sélectionnée :
                y = round(horizPos);
                z = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                x = ceil(get(hslider_3sagittal,'Value'));
        end
        
        % Protection : abort if user clicked outside the axes
        dim = size(volume);
        if ( (1<=x && x<=dim(1)) && (1<=y && y<=dim(2)) && (1<=z && z<=dim(3)) )
            % List of positions
            if SorD == 1 % left click : source
                optodes.srcPos = [optodes.srcPos; x y z];
            elseif SorD == 3 % right click : detector
                optodes.detPos = [optodes.detPos; x y z];
            end % middle click : do nothing

            % Display optode on graphics
            updateDisplay(1:3);

            % Update list of positions in edit text
            set(hedit_srcPos,'String',num2str(optodes.srcPos));
            set(hedit_detPos,'String',num2str(optodes.detPos));
            
        end
        
    end




  % --- Update the list and display of optodes positions to include the
    % newly enetered optode (source or detector)
    function updateEdit(hObject,SorD)
        % hObject : handles to the edit text object that called the
        % function
        % SorD : 's' (source) or 'd' (detector)
        
        % Get the coordinates entered
        coords = str2num(get(hObject,'String'));
        if ~isempty(coords)
            allx = round(coords(:,1));
            ally = round(coords(:,2));
            allz = round(coords(:,3));
            dim = size(volume); % Make sure the new position is inside the volume
            % (and while we are at it, this will also check that the input is
            % really a number (char, Inf or NaN will not respect this condition)
            if isreal(sum(allx(:))+sum(ally(:))+sum(allz)) && ...
               isfinite(sum(allx(:))+sum(ally(:))+sum(allz)) && ...
                (1 <= min(allx)) && (max(allx) <= dim(1)) && ...
                (1 <= min(ally)) && (max(ally) <= dim(2)) && ...
                (1 <= min(allz)) && (max(allz) <= dim(3))
                    % Update list of positions
                    if strcmp(SorD,'s')
                        optodes.srcPos = [allx ally allz];
                    elseif strcmp(SorD,'d')
                        optodes.detPos = [allx ally allz];
                    end
                    % Display optodes on graphics
                    updateDisplay(1:3);
            else  % Do not consider the last change
                if strcmp(SorD,'s')
                    set(hObject,'String',num2str(optodes.srcPos));
                elseif strcmp(SorD,'d')
                    set(hObject,'String',num2str(optodes.detPos));
                end
                    
            end

        elseif ~isempty(get(hObject,'String')) % Unless user did delete all positions,
            % do not consider the last change
            if strcmp(SorD,'s')
                set(hObject,'String',num2str(optodes.srcPos));
            elseif strcmp(SorD,'d')
                set(hObject,'String',num2str(optodes.detPos));
            end
        end            
    end
end