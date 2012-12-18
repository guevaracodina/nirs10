function out = nirs_run_liom_projection(job)
%This function intends to project the contrasts onto 2-D planes, try to
%locate the clinic epileptic focus.
%Ke Peng, 
%2012-08-08, version 0.1, Function created
%**************************************************************************
%   Detailed explanation goes here

 
SPMmat = {};
%Load SPM mat file
try
    SPMmat = job.SPMmat;
    if ~isempty(SPMmat{1})
        nsubj = length(SPMmat);
        SPMok = 1;
    else
        nsubj = 1;
        SPMok = 0;
    end
catch
    SPMok = 0;
    nsubj = 1;
end

%Read configurations from cfg file

projection_cfg.contrasts = job.contrasts_selected;

if isfield(job.mask_option, 'mask_none')
    projection_cfg.mask.mask_choice = 0;
else
    if isfield(job.mask_option, 'mask_contrast')
        projection_cfg.mask.mask_choice = 1;
        projection_cfg.mask.type = 'contrast';
        projection_cfg.mask.contrast = job.mask_option.mask_contrast.mask_contrast_select;
        projection_cfg.mask.unc_p_value = job.mask_option.mask_contrast.mask_contrast_unc_p_value;
        projection_cfg.mask.natureofmask = job.mask_option.mask_contrast.mask_contrast_nature_mask;
    elseif isfield(job.mask_option, 'mask_image')
        projection_cfg.mask.mask_choice = 2;
        projection_cfg.mask.type = 'image';
        projection_cfg.mask.image = job.mask_option.mask_image.mask_image_select;
        projection_cfg.mask.natureofmask = job.mask_option.mask_image.mask_image_nature_mask;
    else
    end
end

if isfield(job.title_comparison, 'title_default')
    projection_cfg.title.option = 1;
elseif isfield(job.title_comparison, 'title_input')
    projection_cfg.title.option = 2;
    projection_cfg.title.TitleName = job.title_comparison.title_input.title_comparison_input;
else
end

if isfield(job.p_value_adjustment, 'p_value_FWE')
    projection_cfg.p_value_adjustment.option = 1;
    projection_cfg.p_value_adjustment.name = 'FWE'; 
    projection_cfg.p_value_adjustment.p_value = job.p_value_adjustment.p_value_FWE.p_value_FWE_p_input;
elseif isfield(job.p_value_adjustment, 'p_value_none')
    projection_cfg.p_value_adjustment.option = 2;
    projection_cfg.p_value_adjustment.name = 'none'; 
    projection_cfg.p_value_adjustment.threshold = job.p_value_adjustment.p_value_none.p_value_none_threshold_input;
else
end

projection_cfg.extent_threshold = job.extent_threshold;

%**************************************************************************
%Excecute

for Idx = 1 : nsubj

    try
        [swd,swd_name,swd_ext] = fileparts(SPMmat{Idx});
        clear swd_name swd_ext;
        try
            load(fullfile(swd,'SPM.mat'));
            SPM.swd = swd;
        catch
            SPMok = 0;
            error(['Cannot read ' fullfile(swd,'SPM.mat')]);
        end
        
        %generate xSPM structure
        try
            [SPM,xSPM] = nirs_liom_get_xSPM(SPM,projection_cfg);
        catch
            SPMok = 0;
            error(['Cannot generate xSPM structure']);
        end
        
        if SPMok
           
            if isempty(xSPM) 
                out = {[],[],[]};
                return;
            end

            %-Ensure pwd = swd so that relative filenames are valid
            %----------------------------------------------------------------------
            cd(SPM.swd)

            %-Get space information
            %======================================================================
            M         = SPM.xVol.M;
            DIM       = SPM.xVol.DIM;

            %-Space units
            %----------------------------------------------------------------------
            try
                try
                    units = SPM.xVol.units;
                catch
                    units = xSPM.units;
                end
            catch
                try
                    if strcmp(spm('CheckModality'),'EEG')
                        datatype = {...
                            'Volumetric (2D/3D)',...
                            'Scalp-Time',...
                            'Scalp-Frequency',...
                            'Time-Frequency',...
                            'Frequency-Frequency'};
                        selected = spm_input('Data Type: ','+1','m',datatype);
                        datatype = datatype{selected};
                    else
                        datatype = 'Volumetric (2D/3D)';
                    end
                catch
                    datatype     = 'Volumetric (2D/3D)';
                end

                switch datatype
                    case 'Volumetric (2D/3D)'
                        units    = {'mm' 'mm' 'mm'};
                    case 'Scalp-Time'
                        units    = {'mm' 'mm' 'ms'};
                    case 'Scalp-Frequency'
                        units    = {'mm' 'mm' 'Hz'};
                    case 'Time-Frequency'
                        units    = {'Hz' 'ms' ''};
                    case 'Frequency-Frequency'
                        units    = {'Hz' 'Hz' ''};
                    otherwise
                        error('Unknown data type.');
                end
            end
            if DIM(3) == 1, units{3} = ''; end
            xSPM.units      = units;
            SPM.xVol.units  = units;


            %-Setup Results User Interface; Display MIP, design matrix & parameters
            %======================================================================
            
            SVNid = '$Rev: 4209 $'; 
            SPMid = spm('FnBanner',mfilename,SVNid);
            [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
            spm_clf('Satellite')
            FS    = spm('FontSizes');

            %-Setup Maximum intensity projection (MIP) & register
            %----------------------------------------------------------------------
            hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
            hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);

            %spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
            if xSPM.STAT == 'P'
                str = xSPM.STATstr;
            else
                str = ['SPM\{',xSPM.STATstr,'\}'];
            end
            text(240,260,str,...
                'Interpreter','TeX',...
                'FontSize',FS(14),'Fontweight','Bold',...
                'Parent',hMIPax)


            %-Print comparison title
            %----------------------------------------------------------------------
            hTitAx = axes('Parent',Fgraph,...
                'Position',[0.02 0.95 0.96 0.02],...
                'Visible','off');

            text(0.5,0,xSPM.title,'Parent',hTitAx,...
                'HorizontalAlignment','center',...
                'VerticalAlignment','baseline',...
                'FontWeight','Bold','FontSize',FS(14))


            %-Print SPMresults: Results directory & thresholding info
            %----------------------------------------------------------------------
            hResAx = axes('Parent',Fgraph,...
                'Position',[0.05 0.55 0.45 0.05],...
                'DefaultTextVerticalAlignment','baseline',...
                'DefaultTextFontSize',FS(9),...
                'DefaultTextColor',[1,1,1]*.7,...
                'Units','points',...
                'Visible','off');
            AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
            h     = text(0,24,'SPMresults:','Parent',hResAx,...
                'FontWeight','Bold','FontSize',FS(14));
            text(get(h,'Extent')*[0;0;1;0],24,spm_str_manip(SPM.swd,'a30'),'Parent',hResAx)
            try
                thresDesc = xSPM.thresDesc;
                text(0,12,sprintf('Height threshold %c = %0.6f  {%s}',xSPM.STAT,xSPM.u,thresDesc),'Parent',hResAx)
            catch
                text(0,12,sprintf('Height threshold %c = %0.6f',xSPM.STAT,xSPM.u),'Parent',hResAx)
            end
            text(0,00,sprintf('Extent threshold k = %0.0f voxels',xSPM.k), 'Parent',hResAx)


            %-Plot design matrix
            %----------------------------------------------------------------------
                        
            hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
            hDesMtxIm = image((SPM.xX.nKX + 1)*32);
            xlabel('Design matrix')
            set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')',...
                'UserData',struct(...
                'X',        SPM.xX.xKXs.X,...
                'fnames',   {reshape({SPM.xY.VY.fname},size(SPM.xY.VY))},...
                'Xnames',   {SPM.xX.name}))

            %-Plot contrasts
            %----------------------------------------------------------------------
            nPar   = size(SPM.xX.X,2);
            xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
            nCon   = length(xSPM.Ic);
            xCon   = SPM.xCon;
            if nCon
                dy     = 0.15/max(nCon,2);
                hConAx = axes('Position',[0.65 (0.80 + dy*.1) 0.25 dy*(nCon-.1)],...
                    'Tag','ConGrphAx','Visible','off');
                title('contrast(s)')
                htxt   = get(hConAx,'title');
                set(htxt,'Visible','on','HandleVisibility','on')
            end

            for ii = nCon:-1:1
                axes('Position',[0.65 (0.80 + dy*(nCon - ii +.1)) 0.25 dy*.9])
                if xCon(xSPM.Ic(ii)).STAT == 'T' && size(xCon(xSPM.Ic(ii)).c,2) == 1

                    %-Single vector contrast for SPM{t} - bar
                    %--------------------------------------------------------------
                    yy = [zeros(1,nPar);repmat(xCon(xSPM.Ic(ii)).c',2,1);zeros(1,nPar)];
                    h  = patch(xx,yy,[1,1,1]*.5);
                    set(gca,'Tag','ConGrphAx',...
                        'Box','off','TickDir','out',...
                        'XTick',spm_DesRep('ScanTick',nPar,10) - 0.5,'XTickLabel','',...
                        'XLim', [0,nPar],...
                        'YTick',[-1,0,+1],'YTickLabel','',...
                        'YLim',[min(xCon(xSPM.Ic(ii)).c),max(xCon(xSPM.Ic(ii)).c)] +...
                        [-1 +1] * max(abs(xCon(xSPM.Ic(ii)).c))/10  )

                else

                    %-F-contrast - image
                    %--------------------------------------------------------------
                    h = image((xCon(xSPM.Ic(ii)).c'/max(abs(xCon(xSPM.Ic(ii)).c(:)))+1)*32);
                    set(gca,'Tag','ConGrphAx',...
                        'Box','on','TickDir','out',...
                        'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
                        'XLim', [0,nPar]+0.5,...
                        'YTick',[0:size(SPM.xCon(xSPM.Ic(ii)).c,2)]+0.5,...
                        'YTickLabel','',...
                        'YLim', [0,size(xCon(xSPM.Ic(ii)).c,2)]+0.5 )

                end
                ylabel(num2str(xSPM.Ic(ii)))
                set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
                    'UserData', struct( 'i',        xSPM.Ic(ii),...
                    'h',        htxt,...
                    'xCon',     xCon(xSPM.Ic(ii))))
            end


            %-Store handles of results section Graphics window objects
            %----------------------------------------------------------------------
            H  = get(Fgraph,'Children');
            H  = findobj(H,'flat','HandleVisibility','on');
            H  = findobj(H);
            Hv = get(H,'Visible');
            set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))

            %-Finished results setup
            %----------------------------------------------------------------------
            out.SPMmat = SPM;
            out.xSPMmat = xSPM;
            spm('Pointer','Arrow')
            
        end
        
    catch
    end
end

end

