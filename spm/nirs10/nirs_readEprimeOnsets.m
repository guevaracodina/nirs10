function out = nirs_readEprimeOnsets(varargin)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-09

persistent PF FS WS eprimewin   % GUI related constants
persistent handles

if (~strcmp(varargin{1},'init') &&...
        ~strcmp(varargin{1},'sav') &&...
        ~strcmp(varargin{1},'trig') &&...
        ~strcmp(varargin{1},'stim') &&...
        ~strcmp(varargin{1},'cond') &&...
        ~strcmp(varargin{1},'rem_trig') &&...
        ~strcmp(varargin{1},'rem_stim') &&...
        ~strcmp(varargin{1},'rem_cond'))
    
    excelp = varargin{1,1}.subj.excelp;
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
        nirs_readEprimeOnsets(ui_action,excelp,sDtp);
        
        % Since we are using persistent variables we better make sure
        % there is no-one else out there.
        if ~exist('eprimewin')%isempty(eprimewin)
            figure(eprimewin);
            set(eprimewin,'Visible',visibility);
            return
        end
        
        S0   = spm('WinSize','0',1);
        WS   = spm('WinScale');
        FS   = spm('FontSizes');
        PF   = spm_platform('fonts');
        rect = [100 100 400 400].*WS;
        
        eprimewin = figure('IntegerHandle','off',...
            'Name',sprintf('%s%s','Get needful columns from Eprime',...
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
        uicontrol(eprimewin,'Style','Frame','BackgroundColor',spm('Colour'),...
            'Position',[10 10 380 344].*WS);
        uicontrol(eprimewin,'Style','Text','String','Pick columns for:',...
            'Position',[25 329 355 20].*WS,...
            'ForegroundColor','k','BackgroundColor',spm('Colour'),...
            'FontName',PF.times,'FontAngle','Italic')
        %         uicontrol(eprimewin,'Style','Text','String','Choose anatomical points, sources and detectors',...
        %             'Position',[25 200 120 20].*WS,...
        %             'ForegroundColor','k','BackgroundColor',spm('Colour'),...
        %             'FontName',PF.times)
        
        %display infos from text document
        uicontrol(eprimewin,'Style','listbox',...
            'ToolTipString','Header list',...
            'String',handles.Headn,...
            'Value',1,...
            'Position',[25 50 120 218].*WS,...
            'Tag','txt_header');
        
        % Select buttons
        uicontrol(eprimewin,'String','>>',...
            'Position',[170 221 22 22].*WS,...
            'ToolTipString','Select Trigger',...
            'CallBack','nirs_readEprimeOnsets(''trig'');',...
            'Tag','button_trig');
        uicontrol(eprimewin,'String','>>',...
            'Position',[170 145 22 22].*WS,...
            'ToolTipString','Select StimOnsets',...
            'CallBack','nirs_readEprimeOnsets(''stim'');',...
            'Tag','button_stim');
        uicontrol(eprimewin,'String','>>',...
            'Position',[170 69 22 22].*WS,...
            'ToolTipString','Select ConditionTags',...
            'CallBack','nirs_readEprimeOnsets(''cond'');',...
            'Tag','button_cond');
        
        uicontrol(eprimewin,'String','<<',...
            'Position',[170 199 22 22].*WS,...
            'ToolTipString','Remove Trigger',...
            'CallBack','nirs_readEprimeOnsets(''rem_trig'');',...
            'Tag','remove_trig');
        uicontrol(eprimewin,'String','<<',...
            'Position',[170 123 22 22].*WS,...
            'ToolTipString','Remove StimOnsets',...
            'CallBack','nirs_readEprimeOnsets(''rem_stim'');',...
            'Tag','remove_stim');
        uicontrol(eprimewin,'String','<<',...
            'Position',[170 47 22 22].*WS,...
            'ToolTipString','Remove ConditionTags',...
            'CallBack','nirs_readEprimeOnsets(''rem_cond'');',...
            'Tag','remove_cond');
        
        uicontrol(eprimewin,'String','Save',...
            'Position',[35 360 60 22].*WS,...
            'ToolTipString','You are ok with your choices',...
            'CallBack','nirs_readEprimeOnsets(''sav'');',...
            'Tag','button_gen');
        
        %selected
        uicontrol(eprimewin,'Style','Text','String','Triggers',...
            'Position',[217 250 120 20].*WS,...
            'ForegroundColor','k','BackgroundColor',spm('Colour'),...
            'FontName',PF.times,'FontAngle','Italic')
        uicontrol(eprimewin,'Style','listbox',...
            'ToolTipString','Selected Trigger',...
            'String',handles.Trign_s,...
            'Value',1,...
            'Position',[217 210 150 40].*WS,...
            'Tag','s_trig');
        uicontrol(eprimewin,'Style','Text','String','Stim Onsets',...
            'Position',[217 170 120 20].*WS,...
            'ForegroundColor','k','BackgroundColor',spm('Colour'),...
            'FontName',PF.times,'FontAngle','Italic')
        uicontrol(eprimewin,'Style','listbox',...
            'ToolTipString','Selected StimOnsets',...
            'String',handles.Stimn_s,...
            'Value',1,...
            'Position',[217 130 150 40].*WS,...
            'Tag','s_stim');
        uicontrol(eprimewin,'Style','Text','String','Condition Tags',...
            'Position',[217 90 120 20].*WS,...
            'ForegroundColor','k','BackgroundColor',spm('Colour'),...
            'FontName',PF.times,'FontAngle','Italic')
        uicontrol(eprimewin,'Style','listbox',...
            'ToolTipString','Selected ConditionTags',...
            'String',handles.Condn_s,...
            'Value',1,...
            'Position',[217 50 150 40].*WS,...
            'Tag','s_cond');
        out =0;
        
    case 'init'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read Brainsight text file and build NIRS matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        excelp = varargin{1,2};
        sDtp = varargin{1,3};
        %         handles = nirs_criugm_readbrainsight(jobRB);
        [data,header]=xlsread(excelp);
        handles.sDtp = varargin{1,3};
        
        % raw lists of header from excelEprime
        handles.Headn = header(1,:)';
        % address of columns for each header in excel sheet
        %         handles.Head_clp = xlsColNum2Str(ind(handles.Headn));
        
        
        handles.Trign_s = {'Start5Trig.OnsetTime'};
        handles.Stimn_s = {'Stimulus.OnsetTime'};
        handles.Condn_s = {'Tag'};
        
        handles.Trig_s_clp = {};
        handles.Stim_s_clp = {};
        handles.Cond_s_clp = {};
        out =0;
        
    case 'trig'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','txt_header'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_trig');
        % add selected item in selected list
        if ~any(strcmp(handles.Trign_s,handles.Headn{htext_value,1}))
            count = length(handles.Trign_s)+1;
            handles.Trign_s{count,1} = handles.Headn{htext_value,1};
            handles.Trig_s_clp{count,1} = xlsColNum2Str(find(strcmp(handles.Headn,handles.Headn{htext_value,1})));
            set(hselected,'String',handles.Trign_s);
        end
        out =0;
    case 'stim'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','txt_header'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_stim');
        % add selected item in selected list
        if ~any(strcmp(handles.Stimn_s,handles.Headn{htext_value,1}))
            count = length(handles.Stimn_s)+1;
            handles.Stimn_s{count,1} = handles.Headn{htext_value,1};
            handles.Stim_s_clp{count,1} = xlsColNum2Str(find(strcmp(handles.Headn,handles.Headn{htext_value,1})));
            set(hselected,'String',handles.Stimn_s);
        end
        out =0;
    case 'cond'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','txt_header'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_cond');
        % add selected item in selected list
        if ~any(strcmp(handles.Condn_s,handles.Headn{htext_value,1}))
            count = length(handles.Condn_s)+1;
            handles.Condn_s{count,1} = handles.Headn{htext_value,1};
            handles.Cond_s_clp{count,1} = xlsColNum2Str(find(strcmp(handles.Headn,handles.Headn{htext_value,1})));
            set(hselected,'String',handles.Condn_s);
        end
        out =0;
        
    case 'rem_trig'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','s_trig'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_trig');
        if size(handles.Trign_s,1)==1
            bool =~strcmp(handles.Trign_s,handles.Trign_s{htext_value,1});
            Trign_cor_s = handles.Trign_s(bool);
            clear handles.Trign_s
            handles.Trign_s = Trign_cor_s;
        else
            handles.Trign_s = {};
        end
        set(hselected,'String',handles.Trign_s);
        out =0;
        
    case 'rem_stim'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','s_stim'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_stim');
        if size(handles.Stimn_s,1)==1
            bool = ~strcmp(handles.Stimn_s,handles.Stimn_s{htext_value,1});
            Stimn_cor_s = handles.Stimn_s(bool);
            clear handles.Stimn_s
            handles.Stimn_s = Stimn_cor_s;
        else
            handles.Stimn_s = {};
        end
        set(hselected,'String',handles.Stimn_s);
        out =0;
        
    case 'rem_cond'
        htext_value = get(findobj(get(eprimewin,'Children'),'Tag','s_cond'),'value');
        hselected = findobj(get(eprimewin,'Children'),'Tag','s_cond');
        if size(handles.Condn_s,1)==1
            bool = ~strcmp(handles.Condn_s,handles.Condn_s{htext_value,1});
            Condn_cor_s = handles.Condn_s(bool);
            clear handles.Condn_s
            handles.Condn_s = Condn_cor_s;
        else
            handles.Condn_s = {};
        end
        set(hselected,'String',handles.Condn_s);
        out =0;
        
    case 'sav'
        %NIRS
        eprime =[];
        eprime.col.trig.n = handles.Trign_s';   % Trigger
        eprime.col.stim.n = handles.Stimn_s';   % StimOnsets
        eprime.col.cond.n = handles.Condn_s';   % ConditionTags
        
        %         NIRS.Cf.H.S.N = size(handles.Sn_s,1);
        %         NIRS.Cf.H.D.N = size(handles.Dn_s,1);
        %         NIRS.Cf.H.Q.N = size(handles.Qn_s,1);
        %
        eprime.col.trig.pos = handles.Trig_s_clp;
        eprime.col.stim.pos = handles.Stim_s_clp;
        eprime.col.cond.pos = handles.Cond_s_clp;
        
        %         nirs.cf.h.d.r.o.mm.p = cell2mat(handles.dp_s_rom');
        %         nirs.cf.h.q.r.o.mm.p = cell2mat(handles.qp_s_rom');
        
        %         % saving names of fiducials for next subjects
        %         load(fullfile(fileparts(which('nirs10')),'nirs10_templates','fid_names.mat'));
        %
        %         %Nasion
        %         bool =strcmpi(handles.Fn_s{1,1},fid_names.names_na);
        %         if sum(bool)==0
        %             fid_names.names_na{1,size(fid_names.names_na,2)+1} = handles.Fn_s{1,1};
        %         end
        %         %LeftEar
        %         bool =strcmpi(handles.Fn_s{2,1},fid_names.names_le);
        %         if sum(bool)==0
        %             fid_names.names_le{1,size(fid_names.names_le,2)+1} = handles.Fn_s{2,1};
        %         end
        %         %RightEar
        %         bool =strcmpi(handles.Fn_s{3,1},fid_names.names_re);
        %         if sum(bool)==0
        %             fid_names.names_re{1,size(fid_names.names_re,2)+1} = handles.Fn_s{3,1};
        %         end
        
        
        %         NIRS.Cf.H.O.n = On';
        %         NIRS.Cf.H.O.r.o.mm.p = Op_rom';
        %         NIRS.Cf.H.O.N = size(NIRS.Cf.H.O.n,2);
        
        save(fullfile(handles.sDtp,'NIRS_eprime.mat'),'eprime');
        out = 1;
        delete(gcf);% window is closed once the button is pushed
    otherwise
        a=1;
        
end
return