function out = nirs_criugm_getHelmet2(varargin)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Cl�ment Bonn�ry
% 2010-09

persistent PF FS WS iugmwin   % GUI related constants
persistent handles

if (~strcmp(varargin{1},'init') &&...
        ~strcmp(varargin{1},'sav') &&...
        ~strcmp(varargin{1},'fid_na') &&...
        ~strcmp(varargin{1},'fid_le') &&...
        ~strcmp(varargin{1},'fid_re') &&...
        ~strcmp(varargin{1},'src') &&...
        ~strcmp(varargin{1},'det') &&...
        ~strcmp(varargin{1},'poi') &&...
        ~strcmp(varargin{1},'rem_fid_na') &&...
        ~strcmp(varargin{1},'rem_fid_le') &&...
        ~strcmp(varargin{1},'rem_fid_re') &&...
        ~strcmp(varargin{1},'rem_src') &&...
        ~strcmp(varargin{1},'rem_det') &&...
        ~strcmp(varargin{1},'rem_poi') &&...
        ~strcmp(varargin{1},'all_src') &&...
        ~strcmp(varargin{1},'all_det') &&...
        ~strcmp(varargin{1},'all_poi'))
    
    staxp = varargin{1,1}.subj.helmet.staxp;
    sDtp = varargin{1,1}.subj.sDtp;
    ui_action = 'welcome';
    
else
    ui_action = varargin{1};
end

switch lower(ui_action)
    case 'welcome'
        visibility = 'On';
        
        % Get all default values (these may effect GUI)
        ui_action = 'init';
        nirs_criugm_getHelmet2(ui_action,staxp,sDtp);
        
        % Since we are using persistent variables we better make sure
        % there is no-one else out there.
        if ~exist('iugmwin')%isempty(iugmwin)
            figure(iugmwin);
            set(iugmwin,'Visible',visibility);
            return
        end
        
        S0   = spm('WinSize','0',1);
        WS   = spm('WinScale');
        FS   = spm('FontSizes');
        PF   = spm_platform('fonts');
        rect = [100 100 400 400].*WS;
        
        iugmwin = figure('IntegerHandle','off',...
            'Name',sprintf('%s%s','Get positions from Brainsight',...
            spm('GetUser',' (%s)')),...
            'NumberTitle','off',...
            'Tag','FieldMap',...
            'Position',[S0(1),S0(2),0,0] + rect,...
            'Resize','off',...
            'Pointer','Arrow',...
            'Color',[1 1 1]*.8,...
            'MenuBar','none',...
            'DefaultTextFontName',PF.helvetica,...
            'DefaultTextFontSize',FS(10),...
            'DefaultUicontrolFontName',PF.helvetica,...
            'DefaultUicontrolFontSize',FS(10),...
            'HandleVisibility','on',...
            'Visible',visibility,...
            'DeleteFcn','FieldMap(''Quit'');');
        
        % Frames and text
        uicontrol(iugmwin,'Style','Frame','BackgroundColor',spm('Colour'),...
            'Position',[10 10 380 344].*WS);
        uicontrol(iugmwin,'Style','Text','String','Pick positions you want to coregister',...
            'Position',[25 329 355 20].*WS,...
            'ForegroundColor','k','BackgroundColor',spm('Colour'),...
            'FontName',PF.times,'FontAngle','Italic')
        %         uicontrol(iugmwin,'Style','Text','String','Choose anatomical points, sources and detectors',...
        %             'Position',[25 200 120 20].*WS,...
        %             'ForegroundColor','k','BackgroundColor',spm('Colour'),...
        %             'FontName',PF.times)
        
        %display infos from text document
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Anatomical Points',...
            'String',handles.On,...
            'Value',1,...
            'Position',[25 253 120 66].*WS,...
            'Tag','txt_ana');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Sources',...
            'String',handles.Sn,...
            'Value',1,...
            'Position',[25 177 120 66].*WS,...
            'Tag','txt_src');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Detectors',...
            'String',handles.Dn,...
            'Value',1,...
            'Position',[25 101 120 66].*WS,...
            'Tag','txt_det');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Other Points of Interest',...
            'String',handles.On,...
            'Value',1,...
            'Position',[25 25 120 66].*WS,...
            'Tag','txt_poi');
        
        % Select buttons
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 297 22 22].*WS,...
            'ToolTipString','Select nasion fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''fid_na'');',...
            'Tag','button_na');
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 275 22 22].*WS,...
            'ToolTipString','Select left ear fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''fid_le'');',...
            'Tag','button_le');
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 253 22 22].*WS,...
            'ToolTipString','Select right ear fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''fid_re'');',...
            'Tag','button_re');
        
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 221 22 22].*WS,...
            'ToolTipString','Select Sources',...
            'CallBack','nirs_criugm_getHelmet2(''src'');',...
            'Tag','button_src');
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 145 22 22].*WS,...
            'ToolTipString','Select Detectors',...
            'CallBack','nirs_criugm_getHelmet2(''det'');',...
            'Tag','button_det');
        uicontrol(iugmwin,'String','>>',...
            'Position',[170 69 22 22].*WS,...
            'ToolTipString','Select Points of Interest',...
            'CallBack','nirs_criugm_getHelmet2(''poi'');',...
            'Tag','button_poi');
        
        uicontrol(iugmwin,'String','<<',...
            'Position',[192 297 22 22].*WS,...
            'ToolTipString','Remove nasion fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''rem_fid_na'');',...
            'Tag','remove_fid_na');
        uicontrol(iugmwin,'String','<<',...
            'Position',[192 275 22 22].*WS,...
            'ToolTipString','Remove left ear fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''rem_fid_le'');',...
            'Tag','remove_fid_le');
        uicontrol(iugmwin,'String','<<',...
            'Position',[192 253 22 22].*WS,...
            'ToolTipString','Remove right ear fiducial',...
            'CallBack','nirs_criugm_getHelmet2(''rem_fid_re'');',...
            'Tag','remove_fid_re');
        
        uicontrol(iugmwin,'String','<<',...
            'Position',[170 199 22 22].*WS,...
            'ToolTipString','Remove Sources',...
            'CallBack','nirs_criugm_getHelmet2(''rem_src'');',...
            'Tag','remove_src');
        uicontrol(iugmwin,'String','<<',...
            'Position',[170 123 22 22].*WS,...
            'ToolTipString','Remove Detectors',...
            'CallBack','nirs_criugm_getHelmet2(''rem_det'');',...
            'Tag','remove_det');
        uicontrol(iugmwin,'String','<<',...
            'Position',[170 47 22 22].*WS,...
            'ToolTipString','Remove Points of Interest',...
            'CallBack','nirs_criugm_getHelmet2(''rem_poi'');',...
            'Tag','remove_poi');
        
        uicontrol(iugmwin,'String','All',...
            'Position',[170 177 22 22].*WS,...
            'ToolTipString','Select all Sources',...
            'CallBack','nirs_criugm_getHelmet2(''all_src'');',...
            'Tag','selectall_src');
        uicontrol(iugmwin,'String','All',...
            'Position',[170 101 22 22].*WS,...
            'ToolTipString','Select all Detectors',...
            'CallBack','nirs_criugm_getHelmet2(''all_det'');',...
            'Tag','selectall_det');
        uicontrol(iugmwin,'String','All',...
            'Position',[170 25 22 22].*WS,...
            'ToolTipString','Select all Points of Interest',...
            'CallBack','nirs_criugm_getHelmet2(''all_poi'');',...
            'Tag','selectall_poi');
        
        uicontrol(iugmwin,'String','Save',...
            'Position',[35 360 60 22].*WS,...
            'ToolTipString','You are ok with your choices',...
            'CallBack','nirs_criugm_getHelmet2(''sav'');',...
            'Tag','button_gen');
        
        %selected
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Selected nasion fiducial',...
            'String',handles.na_s,...
            'Value',1,...
            'Position',[217 297 120 18].*WS,...
            'Tag','s_fid_na');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Selected left ear fiducial',...
            'String',handles.LE_s,...
            'Value',1,...
            'Position',[217 275 120 18].*WS,...
            'Tag','s_fid_le');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Selected right ear fiducial',...
            'String',handles.RE_s,...
            'Value',1,...
            'Position',[217 253 120 18].*WS,...
            'Tag','s_fid_re');
        
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Selected Sources',...
            'String',handles.Sn_s,...
            'Value',1,...
            'Position',[217 177 120 66].*WS,...
            'Tag','s_src');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Selected Detectors',...
            'String',handles.Dn_s,...
            'Value',1,...
            'Position',[217 101 120 66].*WS,...
            'Tag','s_det');
        uicontrol(iugmwin,'Style','listbox',...
            'ToolTipString','Other Selected Points of Interest',...
            'String',handles.Qn_s,...
            'Value',1,...
            'Position',[217 25 120 66].*WS,...
            'Tag','s_poi');
        out =0;
        
    case 'init'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read Brainsight text file and build NIRS matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jobRB.staxp = varargin{1,2};
        jobRB.sDtp = varargin{1,3};
        handles = nirs_criugm_readbrainsight(jobRB);
        handles.sDtp = varargin{1,3};
        
        % raw lists of points from stereotaxic logiciel (here Brainsight)
        handles.On = handles.Cf.H.O.n';
        handles.Sn = handles.Cf.H.S.n';
        handles.Dn = handles.Cf.H.D.n';
        % coordinates of these points
        handles.Op_rom = handles.Cf.H.O.r.o.mm.p';
        handles.Sp_rom = handles.Cf.H.S.r.o.mm.p';
        handles.Dp_rom = handles.Cf.H.D.r.o.mm.p';
        
        handles.Fn_s = {}; % as it is a catch code, the automatic selection as already failed
        handles.na_s = {};
        handles.LE_s = {};
        handles.RE_s = {};
        handles.Sn_s = handles.Cf.H.S.n';
        handles.Dn_s = handles.Cf.H.D.n';
        handles.Qn_s = {};
        
        handles.Fp_s_rom = {};
        handles.Sp_s_rom = handles.Cf.H.S.r.o.mm.p';
        handles.Dp_s_rom = handles.Cf.H.D.r.o.mm.p';
        handles.Qp_s_rom = {};
        out =0;
        
    case 'fid_na'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_ana'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_na');
        % add selected item in selected list if not already selected
        if ~any(strcmp(handles.na_s,handles.On{htext_value,1}))
            handles.na_s = handles.On{htext_value,1};
            handles.Fp_s_rom{1,1} = handles.Op_rom{htext_value,1};
            handles.Fp_s_rom{1,2} = handles.Op_rom{htext_value,2};
            handles.Fp_s_rom{1,3} = handles.Op_rom{htext_value,3};
            set(hselected,'String',handles.na_s);
        end
        out =0;
    case 'fid_le'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_ana'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_le');
        % add selected item in selected list if not already selected
        if ~any(strcmp(handles.LE_s,handles.On{htext_value,1}))
            handles.LE_s = handles.On{htext_value,1};
            handles.Fp_s_rom{2,1} = handles.Op_rom{htext_value,1};
            handles.Fp_s_rom{2,2} = handles.Op_rom{htext_value,2};
            handles.Fp_s_rom{2,3} = handles.Op_rom{htext_value,3};
            set(hselected,'String',handles.LE_s);
        end
        out =0;
    case 'fid_re'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_ana'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_re');
        % add selected item in selected list if not already selected
        if ~any(strcmp(handles.RE_s,handles.On{htext_value,1}))
            handles.RE_s = handles.On{htext_value,1};
            handles.Fp_s_rom{3,1} = handles.Op_rom{htext_value,1};
            handles.Fp_s_rom{3,2} = handles.Op_rom{htext_value,2};
            handles.Fp_s_rom{3,3} = handles.Op_rom{htext_value,3};
            set(hselected,'String',handles.RE_s);
        end
        out =0;
    case 'src'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_src'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_src');
        % add selected item in selected list
        if ~any(strcmp(handles.Sn_s,handles.Sn{htext_value,1}))
            count = length(handles.Sn_s)+1;
            handles.Sn_s{count,1} = handles.Sn{htext_value,1};
            handles.Sp_s_rom{count,1} = handles.Sp_rom{htext_value,1};
            handles.Sp_s_rom{count,2} = handles.Sp_rom{htext_value,2};
            handles.Sp_s_rom{count,3} = handles.Sp_rom{htext_value,3};
            set(hselected,'String',handles.Sn_s);
        end
        out =0;
    case 'det'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_det'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_det');
        % add selected item in selected list
        if ~any(strcmp(handles.Dn_s,handles.Dn{htext_value,1}))
            count = length(handles.Dn_s)+1;
            handles.Dn_s{count,1} = handles.Dn{htext_value,1};
            handles.Dp_s_rom{count,1} = handles.Dp_rom{htext_value,1};
            handles.Dp_s_rom{count,2} = handles.Dp_rom{htext_value,2};
            handles.Dp_s_rom{count,3} = handles.Dp_rom{htext_value,3};
            set(hselected,'String',handles.Dn_s);
        end
        out =0;
    case 'poi'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','txt_poi'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_poi');
        % add selected item in selected list
        if ~any(strcmp(handles.Qn_s,handles.On{htext_value,1}))
            count = length(handles.Qn_s)+1;
            handles.Qn_s{count,1} = handles.On{htext_value,1};
            handles.Qp_s_rom{count,1} = handles.Op_rom{htext_value,1};
            handles.Qp_s_rom{count,2} = handles.Op_rom{htext_value,2};
            handles.Qp_s_rom{count,3} = handles.Op_rom{htext_value,3};
            set(hselected,'String',handles.Qn_s);
        end
        out =0;
        
    case 'rem_fid_na'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_na');
        handles.na_s = {};
        handles.Fp_s_rom{1,1} = 0;
        handles.Fp_s_rom{1,2} = 0;
        handles.Fp_s_rom{1,3} = 0;
        set(hselected,'String',handles.na_s);
        out =0;
    case 'rem_fid_le'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_le');
        handles.LE_s = {};
        handles.Fp_s_rom{2,1} = 0;
        handles.Fp_s_rom{2,2} = 0;
        handles.Fp_s_rom{2,3} = 0;
        set(hselected,'String',handles.LE_s);
        out =0;
    case 'rem_fid_re'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_fid_re');
        handles.RE_s = {};
        handles.Fp_s_rom{3,1} = 0;
        handles.Fp_s_rom{3,2} = 0;
        handles.Fp_s_rom{3,3} = 0;
        set(hselected,'String',handles.RE_s);
        out =0;
    case 'rem_src'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','s_src'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_src');
        bool =~strcmp(handles.Sn_s,handles.Sn_s{htext_value,1});
        Sn_cor_s = handles.Sn_s(bool);
        count=1;
        for i=1:length(bool)
            if bool(i)==1
                Sp_cor_s_rom(count,1:3) = handles.Sp_s_rom(i,1:3);
                count = count+1;
            end
        end
        clear handles.Sn_s handles.Sp_s_rom
        handles.Sn_s = Sn_cor_s;
        handles.Sp_s_rom = Sp_cor_s_rom;
        set(hselected,'String',handles.Sn_s);
        out =0;
        
    case 'rem_det'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','s_det'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_det');
        bool = ~strcmp(handles.Dn_s,handles.Dn_s{htext_value,1});
        Dn_cor_s = handles.Dn_s(bool);
        count=1;
        for i=1:length(bool)
            if bool(i)==1
                Dp_cor_s_rom(count,1:3) = handles.Dp_s_rom(i,1:3);
                count = count+1;
            end
        end
        clear handles.Dn_s handles.Dp_s_rom
        handles.Dn_s = Dn_cor_s;
        handles.Dp_s_rom = Dp_cor_s_rom;
        set(hselected,'String',handles.Dn_s);
        out =0;
        
    case 'rem_poi'
        htext_value = get(findobj(get(iugmwin,'Children'),'Tag','s_poi'),'value');
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_poi');
        bool = ~strcmp(handles.Qn_s,handles.Qn_s{htext_value,1});
        Qn_cor_s = handles.Qn_s(bool);
        count=1;
        for i=1:length(bool)
            if bool(i)==1
                Qp_cor_s_rom(count,1:3) = handles.Qp_s_rom(i,1:3);
                count = count+1;
            end
        end
        clear handles.Qn_s handles.Qp_s_rom
        handles.Qn_s = Qn_cor_s;
        handles.Qp_s_rom = Qp_cor_s_rom;
        set(hselected,'String',handles.Qn_s);
        out =0;
        
    case 'all_src'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_src');
        handles.Sn_s = handles.Sn;
        handles.Sp_s_rom = handles.Sp_rom;
        set(hselected,'String',handles.Sn_s);
        out =0;
    case 'all_det'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_det');
        handles.Dn_s = handles.Dn;
        handles.Dp_s_rom = handles.Dp_rom;
        set(hselected,'String',handles.Dn_s);
        out =0;
    case 'all_poi'
        hselected = findobj(get(iugmwin,'Children'),'Tag','s_poi');
        handles.Qn_s = handles.On;
        handles.Qp_s_rom = handles.Op_rom;
        set(hselected,'String',handles.Qn_s);
        out =0;
        
    case 'sav'
        handles.Fn_s{1,1} = handles.na_s;
        handles.Fn_s{2,1} = handles.LE_s;
        handles.Fn_s{3,1} = handles.RE_s;
        
        %NIRS
        NIRS =[];
        NIRS.Cf.H.F.n = handles.Fn_s';   % Fiducials
        NIRS.Cf.H.S.n = handles.Sn_s';   % Sources
        NIRS.Cf.H.D.n = handles.Dn_s';   % Detectors
        NIRS.Cf.H.Q.n = handles.Qn_s';   % Points of interest
        
        NIRS.Cf.H.F.N = size(handles.Fn_s,1);
        NIRS.Cf.H.S.N = size(handles.Sn_s,1);
        NIRS.Cf.H.D.N = size(handles.Dn_s,1);
        NIRS.Cf.H.Q.N = size(handles.Qn_s,1);
        
        NIRS.Cf.H.F.r.o.mm.p = cell2mat(handles.Fp_s_rom');
        NIRS.Cf.H.S.r.o.mm.p = cell2mat(handles.Sp_s_rom');
        NIRS.Cf.H.D.r.o.mm.p = cell2mat(handles.Dp_s_rom');
        NIRS.Cf.H.Q.r.o.mm.p = cell2mat(handles.Qp_s_rom');
        
        % saving names of fiducials for next subjects
        load(fullfile(fileparts(which('nirs10')),'nirs10_templates','fid_names.mat'));
        
        %Nasion
        bool =strcmpi(handles.Fn_s{1,1},fid_names.names_na);
        if sum(bool)==0
            fid_names.names_na{1,size(fid_names.names_na,2)+1} = handles.Fn_s{1,1};
        end
        %LeftEar
        bool =strcmpi(handles.Fn_s{2,1},fid_names.names_le);
        if sum(bool)==0
            fid_names.names_le{1,size(fid_names.names_le,2)+1} = handles.Fn_s{2,1};
        end
        %RightEar
        bool =strcmpi(handles.Fn_s{3,1},fid_names.names_re);
        if sum(bool)==0
            fid_names.names_re{1,size(fid_names.names_re,2)+1} = handles.Fn_s{3,1};
        end
        
        
        %         NIRS.Cf.H.O.n = On';
        %         NIRS.Cf.H.O.r.o.mm.p = Op_rom';
        %         NIRS.Cf.H.O.N = size(NIRS.Cf.H.O.n,2);
        
        save(fullfile(handles.sDtp,'NIRS_Cf.mat'),'NIRS');
        out = 1;
        delete(gcf);% window is closed once the button is pushed
    otherwise
        a=1;
end
return