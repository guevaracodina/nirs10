function out = nirs_run_liom_orth_coreg(job)
%==========================================================================
%This function intends to co-register the optode position as well as 
%the activations onto orthogonal T1 image
%Ke Peng,
%2013-07-25, version 0.1, Function created
%==========================================================================
%select xSPM template

clear global opt_info %Clear past records

global opt_info

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        ERad = job.radius_channel;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'OrthCoreg_OK') || job.force_redo)
            
            file_render = job.render_image{1};
            if isfield(NIRS.Cf.H.C, 'w')
                opt_info.render_template = 1;
            else
                opt_info.render_template = 0;
            end
            coreg_layer = job.coreg_layer;
            if isempty(file_render)
                file_render = NIRS.Dt.ana.T1;
                disp('Render image not specified. T1 image is chosen by default');
            end
            opsel = job.coreg_type.NIRS_channels_optodes.channel_optode_select;
            coreg.corr = NIRS.Cf.H;
            if ~isfield(NIRS.Dt.ana, 'wT1')
                anatT1 = NIRS.Dt.ana.T1;
                [pth,nam] = spm_fileparts(deblank(anatT1));
                sn_filename  = fullfile(pth,[nam '_sn.mat']);
                if ~exist(sn_filename,'file')
                    sn_filename  = fullfile(pth,['m' nam '_sn.mat']);
                end
                NIRS.Dt.ana.wT1 = load(sn_filename);
            end
            wT1 = NIRS.Dt.ana.wT1;
            Q = (wT1.VG.mat/wT1.Affine)/wT1.VF.mat;
            coreg.Q = Q;
            opt_info.ERad = ERad;
            
            T1 = NIRS.Dt.ana.T1;
            Vt1 = spm_vol(T1);
            m2v.m = Vt1.mat;
            opt_info.m2v.m = m2v.m;
            [dirT1, fil, ext] = fileparts(T1);
            fwT1 = fullfile(dirT1,['w' fil ext(1:4)]);
            VwT1 = spm_vol(fwT1);
            m2v.w = VwT1.mat;
            opt_info.m2v.w = m2v.w;
            wT1_info.mat = VwT1.mat;
            wT1_info.dim = VwT1.dim;
            opt_info.wT1_info = wT1_info;
            
            [coreg xSPM opt] = nirs_orthPrep(coreg,opsel,opt_info.render_template,coreg_layer,m2v);
            opt_info.type = opt.type;
            opt_info.xSPM = xSPM;%xSPM for display
            opt_info.recoreg = coreg;%recoreg for storing re-coregistration results
            if opt.flag %if channel coordinates have been filled, save them
                NIRS.Cf.H = coreg.corr;
            end
            clear opt
            
            %Codes from spm_sections.m
            Fgraph = spm_figure('GetWin','Graphics');
            spm_results_ui('Clear',Fgraph);
            nirs_orthdraw(xSPM, file_render);
            
            switch job.coreg_type.NIRS_channels_optodes.channel_optode_select
                case 0
                    color_ity = 'Red: Channels onto Scalp; Yellow: Channels onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);                
                case 1
                    color_ity = 'Red: Sources onto Scalp; Yellow: Sources onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                case 2
                    color_ity = 'Red: Detectors onto Scalp; Yellow: Detectors onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                case 3
                    color_ity = 'Red-Channels; Orange-Sources; Yellow-Detectors';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                otherwise
            end
           
            %Add button
            if job.coreg_type.NIRS_channels_optodes.channel_optode_select ~= 3
                opt_info.hText.hb.cr = uicontrol(Fgraph,'Style','PushButton','String','Correspondence','FontSize',12,...
                'ToolTipString',...
                'find corresponding place on cortex/scalp of an optode',...
                'Callback','nirs_orthviews(''correspondence'');',...
                'Interruptible','on','Enable','on',...
                'Position',[380 320 125 030]);
            
                opt_info.hText.hb.op = uicontrol(Fgraph,'Style','PushButton','String','Re-coreg(optode)','FontSize',12,...
                'ToolTipString',...
                'Assign another position to an optode to see its projection on the cortex',...
                'Callback','nirs_orthviews(''coregop'');',...
                'Interruptible','on','Enable','on',...
                'Position',[380 290 125 030]); 
            
                opt_info.hText.hb.crgop = uicontrol(Fgraph,'Style','PushButton','String','Done','FontSize',12,...
                'ToolTipString',...
                'Finish choosing an optode on the scalp',...
                'Callback','nirs_orthviews(''coregop_exe'');',...
                'Interruptible','on','Enable','off',...
                'Position',[505 290 050 030]);
            
                opt_info.hText.hb.crgopu = uicontrol(Fgraph,'Style','PushButton','String','Undo','FontSize',12,...
                'ToolTipString',...
                'Undo the optode coregistration',...
                'Callback','nirs_orthviews(''coregop_undo'');',...
                'Interruptible','on','Enable','off',...
                'Position',[555 290 050 030]);
            
                opt_info.hText.hb.crgopr = uicontrol(Fgraph,'Style','PushButton','String','Redo','FontSize',12,...
                'ToolTipString',...
                'Redo the optode coregistration',...
                'Callback','nirs_orthviews(''coregop_redo'');',...
                'Interruptible','on','Enable','off',...
                'Position',[605 290 050 030]);
            
                opt_info.hText.hb.crgopcf = uicontrol(Fgraph,'Style','PushButton','String','Confirm op recoreg','FontSize',12,...
                'ToolTipString',...
                'Confirm the optode coregistration (non invertable)',...
                'Callback','nirs_orthviews(''coregop_confirm'');',...
                'Interruptible','on','Enable','off',...
                'Position',[505 260 150 030]);            
            
                opt_info.hText.hb.fd = uicontrol(Fgraph,'Style','PushButton','String','Re-coreg(fiducial)','FontSize',12,...
                'ToolTipString',...
                'Assign another position to a fiducial to see the projection of all optodes/channels on the cortex',...
                'Callback','nirs_orthviews(''coregfd'');',...
                'Interruptible','on','Enable','on',...
                'Position',[380 230 125 030]);   
            
                opt_info.hText.hb.crgfd = uicontrol(Fgraph,'Style','PushButton','String','Done','FontSize',12,...
                'ToolTipString',...
                'find corresponding place on cortex/scalp of an optode',...
                'Callback','nirs_orthviews(''coregfd_exe'');',...
                'Interruptible','on','Enable','off',...
                'Position',[505 230 050 030]);
            
                opt_info.hText.hb.crgfdu = uicontrol(Fgraph,'Style','PushButton','String','Undo','FontSize',12,...
                'ToolTipString',...
                'Undo the fiducial coregistration',...
                'Callback','nirs_orthviews(''coregfd_undo'');',...
                'Interruptible','on','Enable','off',...
                'Position',[555 230 050 030]);
            
                opt_info.hText.hb.crgfdr = uicontrol(Fgraph,'Style','PushButton','String','Redo','FontSize',12,...
                'ToolTipString',...
                'Redo the fiducial coregistration',...
                'Callback','nirs_orthviews(''coregfd_redo'');',...
                'Interruptible','on','Enable','off',...
                'Position',[605 230 050 030]);
            
                opt_info.hText.hb.crgfdcf = uicontrol(Fgraph,'Style','PushButton','String','Confirm fd recoreg','FontSize',12,...
                'ToolTipString',...
                'Confirm the fiducial coregistration (non invertable)',...
                'Callback','nirs_orthviews(''coregfd_confirm'');',...
                'Interruptible','on','Enable','off',...
                'Position',[505 200 150 030]);                
            
                opt_info.hText.hi = uicontrol(Fgraph,'Style','Edit','Tag','ByNum',...
                'CallBack','nirs_orthviews(''dc_num'')',...
                'Interruptible','off','BusyAction','Queue',...
                'Position',[430 050 050 020]);
            
                opt_info.hText.hb.bnplus = uicontrol(Fgraph,'Style','PushButton','String','>>','FontSize',12,...
                'ToolTipString',...
                'Go to the next optode/channel',...
                'Callback','nirs_orthviews(''dc_num_plus'');',...
                'Interruptible','on','Enable','on',...
                'Position',[480 050 050 020]);
            
                opt_info.hText.hb.bnminus = uicontrol(Fgraph,'Style','PushButton','String','<<','FontSize',12,...
                'ToolTipString',...
                'Go to the previous optode/channel',...
                'Callback','nirs_orthviews(''dc_num_minus'');',...
                'Interruptible','on','Enable','on',...
                'Position',[380 050 050 020]);
            
                opt_info.hText.hb.blpa = uicontrol(Fgraph,'Style','PushButton','String','LeftPA','FontSize',12,...
                'ToolTipString',...
                'Go to left peri-auricular fiducial',...
                'Callback','nirs_orthviews(''goleftpa'');',...
                'Interruptible','on','Enable','on',...
                'Position',[530 050 050 020]);            
            
                opt_info.hText.hb.bnas = uicontrol(Fgraph,'Style','PushButton','String','Nasion','FontSize',12,...
                'ToolTipString',...
                'Go to nasion fiducial',...
                'Callback','nirs_orthviews(''gonasion'');',...
                'Interruptible','on','Enable','on',...
                'Position',[580 050 050 020]);
            
                opt_info.hText.hb.brpa = uicontrol(Fgraph,'Style','PushButton','String','RightPA','FontSize',12,...
                'ToolTipString',...
                'Go to right peri-auricular fiducial',...
                'Callback','nirs_orthviews(''gorightpa'');',...
                'Interruptible','on','Enable','on',...
                'Position',[630 050 050 020]);
            
                opt_info.hText.hb.bsave = uicontrol(Fgraph,'Style','PushButton','String','Save current coreg','FontSize',12,...
                'ToolTipString',...
                'Save current modified NIRS',...
                'Callback','nirs_orthviews(''orthsave'');',...
                'Interruptible','on','Enable','on',...
                'Position',[380 070 150 030]);
            
%                 opt_info.hText.hb.bload = uicontrol(Fgraph,'Style','PushButton','String','Load previous coreg','FontSize',12,...
%                 'ToolTipString',...
%                 'Load previously saved coreg file',...
%                 'Callback','nirs_orthviews(''orthload'');',...
%                 'Interruptible','on','Enable','on',...
%                 'Position',[530 070 150 030]);
            end
        end
        NIRS.flags.OrthCoreg_OK = 1;
        opt_info.NIRS = NIRS;
        opt_info.frender = file_render;
        opt_info.newNIRSlocation = newNIRSlocation;
        opt_info.OPrecoreg.idx = [];
        opt_info.OPrecoreg.SlpXYZmm = [];
        opt_info.OPrecoreg.CtxXYZmm = [];
        opt_info.OPrecoreg.undo = [];
        save(newNIRSlocation,'NIRS');
        
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not Coregister for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
    
    
    
end
out.NIRSmat = job.NIRSmat;
end
