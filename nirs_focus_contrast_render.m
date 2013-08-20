function nirs_focus_contrast_render(NIRS,dat,brt,TOPO,rend,disp_op)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render(dat,brt,rendfile)
%
% dat      - a struct array of length 1 to 3
%            each element is a structure containing:
%            - XYZ - the x, y & z coordinates of the transformed SPM{.}
%                    values in units of voxels.
%            - t   - the SPM{.} values.
%            - mat - affine matrix mapping from XYZ voxels to MNI.
%            - dim - dimensions of volume from which XYZ is drawn.
% brt      - brightness control:
%            If NaN, then displays using the old style with hot
%            metal for the blobs, and grey for the brain.
%            Otherwise, it is used as a ``gamma correction'' to
%            optionally brighten the blobs up a little.
% rendfile - the file containing the images to render on to (see also
%            spm_surf.m) or a surface mesh file.
%
% TOPO     - the TOPO file generated by NIRS_SPM
%
% 
%__________________________________________________________________________
% 
% spm_render prompts for details of up to three SPM{.}s that are then
% displayed superimposed on the surface of a 'standard' brain.
%
% The first is shown in red, then green then blue.
%
% The blobs which are displayed are the integral of all transformed t
% values, exponentially decayed according to their depth. Voxels that
% are 10mm behind the surface have half the intensity of ones at the
% surface.
%__________________________________________________________________________
% Molecular and Optical Imaging Laboratory
% Ke Peng
% 2012-08-20

%============================Configurations================================
try

thres = disp_op.thres_sel;
sessions = disp_op.sessions;
views = disp_op.views;
chromophore = disp_op.chromophore;
save_dir = disp_op.save_dir;
activation = 1;
    
num     = length(dat);
%define the colours for focus and contrasts display
col_contrast_neg = [0,0,1;0,1,0;0,0,1];%Blue for decrease
col_contrast = eye(3); %red for increase

spm('Pointer','Watch');

%Display Calibration
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

nrow = ceil(length(views)/2);
hght = 0.95;
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

%Create structure to restore focus render map and contrast render map
try
    render_proj.focus = cell(1,length(dat));
    render_proj.all_focus = cell(1,length(rend));
    for z0 = 1 : length(dat)
        render_proj.focus{z0} = struct();
        render_proj.focus{z0}.num = z0;
        render_proj.focus{z0}.f_map = cell(1,length(rend));

    end    
    
    render_proj.contrast = cell(1,length(TOPO.v{4}.s)); %Take Dorsal view to obtain the session number
    for i = 1 : length(render_proj.contrast)
        render_proj.contrast{i} = struct();
        render_proj.contrast{i}.session = i;
        render_proj.contrast{i}.c_map = cell(1,length(rend));
    end
    
catch
    disp('Error when creating render_proj');
end

clear z0 z00

%===================Render the epileptic focus (sixviews)==================

for brain_view=1:length(rend),
    rend{brain_view}.max=0;
    rend{brain_view}.data = cell(size(dat,1),1);
    if issparse(rend{brain_view}.ren),
        % Assume that images have been DCT compressed
        % - the SPM99 distribution was originally too big.
        d = size(rend{brain_view}.ren);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        rend{brain_view}.ren = B1*rend{brain_view}.ren*B2';
        % the depths did not compress so well with
        % a straight DCT - therefore it was modified slightly
        rend{brain_view}.dep = exp(B1*rend{brain_view}.dep*B2')-1;
    end
    rend{brain_view}.ren(rend{brain_view}.ren>=1) = 1;
    rend{brain_view}.ren(rend{brain_view}.ren<=0) = 0;
   
end

mx = zeros(length(rend),1)+eps;
mn = zeros(length(rend),1);

for j=1:length(dat),
    XYZ = dat(j).XYZ;
    t   = dat(j).Z'; %modified for NIRS
    dim = dat(j).DIM; %modified for NIRS
    mat = dat(j).M; %modified for NIRS

    for brain_view=1:length(rend),

        % transform from Talairach space to space of the rendered image
        %------------------------------------------------------------------
        M1  = rend{brain_view}.M*mat;
        zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
        M2  = diag([zm' 1 1]);
        M  = M2*M1;
        cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
               1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
        tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
        off = min(tcor(1:2,:)');
        M2  = spm_matrix(-off+1)*M2;
        M  = M2*M1;
        xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
        d2  = ceil(max(xyz(1:2,:)'));

        % Calculate 'depth' of values
        %------------------------------------------------------------------
        if ~isempty(d2)
            dep = spm_slice_vol(rend{brain_view}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
            z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));

            if ~isfinite(brt), msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
            else msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end
        else
            msk = [];
        end
        
        if ~isempty(msk),

            % Generate an image of the integral of the blob values.
            %--------------------------------------------------------------
            xyz = xyz(:,msk);
            if ~isfinite(brt), t0  = t(msk);
            else
                dst = xyz(3,:) - z1(msk);
                dst = max(dst,0);
                t0  = t(msk).*exp((log(0.5)/10)*dst)';
            end
            X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
            hld = 1; if ~isfinite(brt), hld = 0; end
            X_pos   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{brain_view}.dep),hld);
            msk = find(X_pos<0);
            X_pos(msk) = 0;
        else
            X_pos = zeros(size(rend{brain_view}.dep));
        end

        % Brighten the blobs
        %------------------------------------------------------------------
        if isfinite(brt), X_pos = X_pos.^brt; end

        mx(j) = max([mx(j) max(max(X_pos))]);
        mn(j) = min([mn(j) min(min(X_pos))]);

        rend{brain_view}.data_focus{j} = X_pos;
        
        %Write the view name
        %------------------------------------------------------------------
        render_proj.focus{j}.f_map{brain_view}.data = X_pos;
        [side_hemi spec_hemi] = nirs_get_brain_view(brain_view);
        render_proj.focus{j}.f_map{brain_view}.view = spec_hemi;

    end
end

if ~isempty(dat) %If more than one focus is indicated, add all the focus maps
    for brain_view = 1 : length(rend)
        render_proj.all_focus{brain_view} = struct();
        render_proj.all_focus{brain_view}.view = render_proj.focus{1}.f_map{brain_view}.view;
        render_proj.all_focus{brain_view}.map = zeros(size(render_proj.focus{1}.f_map{brain_view}.data));
        for j=1:length(dat)
            render_proj.all_focus{brain_view}.map = render_proj.all_focus{brain_view}.map + render_proj.focus{j}.f_map{brain_view}.data;
        end
        render_proj.all_focus{brain_view}.map(find(render_proj.all_focus{brain_view}.map > 1)) = 1; %#ok<FNDSB>
    end
end

clear i j
%One focus map (in six views) ready to use

%=====================Render the contrasts (Manually) =====================

for i = sessions
    for j = views
        
        if i == 0 %Group view
            con_XYZ = TOPO.v{j}.g.hb{chromophore}.c{activation}.Tmap;
            con_XYZ = flipud(con_XYZ);
            con_XYZ(find(abs(con_XYZ) < 2.1)) = 0;
            con_name = TOPO.v{j}.g.hb{2}.c{1}.c.name;
            disp(['Projection HbR activation ' con_name ' for view ' j ' in group stage']);
            disp('Only uncorrected results are presented.');
            i0 = 20;
            
        else
            con_XYZ = squeeze(TOPO.v{j}.s{i}.hb{chromophore}.stat_map(activation,:,:)); %To project HbR concentration as a first stage. Only have the first activation.
            activation_type = 1;
            con_XYZ = flipud(con_XYZ); %Upside down

            %Remove massive interpolation
            if isfield(rend{j},'view_mask_2d') % for back-compatibility
                con_XYZ = con_XYZ .* rend{j}.view_mask_2d;
            end


            %Apply the threshold value
            if thres
                if isfield(TOPO.v{j}.s{i}.hb{2}, 'th_z') %Corresponding to HbR
                    if isfield(TOPO.v{j}.s{i}.hb{2}.th_z{1}, 'positive_thz') %Corresponding to the first activation
                        thz_p = TOPO.v{j}.s{i}.hb{2}.th_z{1}.positive_thz;
                        con_p = con_XYZ;
                        con_p((con_XYZ - thz_p) <= 0) = 0;  
                    end
                    if isfield(TOPO.v{j}.s{i}.hb{2}.th_z{1}, 'negative_thz')
                        thz_n = TOPO.v{j}.s{i}.hb{2}.th_z{1}.negative_thz;
                        con_n = con_XYZ;
                        con_n((con_XYZ + thz_n) >= 0) = 0;  
                    end

                    if exist('con_p', 'var') && exist('con_n', 'var')
                        con_XYZ = con_p + con_n;
                        clear con_p con_n
                    elseif exist('con_p', 'var')
                        con_XYZ = con_p;
                        clear con_p
                    elseif exist('con_n', 'var')
                        con_XYZ = con_n;
                        clear con_n;
                    else
                        con_XYZ = zeros(size(con_XYZ));
                    end

                else
                %No corrected image
                    con_XYZ = zeros(size(con_XYZ));
                end
            else
                con_XYZ(find(abs(con_XYZ) < 1.7)) = 0;
            end
            i0 = i;
        end
        
  
        mxmx = max(max(abs(con_XYZ)));
        mnmn = min(min(abs(con_XYZ)));
        
        render_proj.contrast{i0}.c_map{j} = struct();
        render_proj.contrast{i0}.c_map{j}.view_no = j;
        
        [side_hemi spec_hemi] = nirs_get_brain_view(j);
        render_proj.contrast{i0}.c_map{j}.view = spec_hemi;
 
        render_proj.contrast{i0}.c_map{j}.data{1} = con_XYZ;%Modify when adding multiple activations
        render_proj.contrast{i0}.c_map{j}.mxmx{1} = mxmx;
        render_proj.contrast{i0}.c_map{j}.mnmn{1} = mnmn;
        
    end

end

%Contrasts data ready to use

%====================Display both focus and contrasts =====================

if ~isfinite(brt),%not used
    % Old style split colourmap display.
    %----------------------------------------------------------------------
    load Split;
    colormap(split);
    for i=1:length(rend),
        ren = rend{i}.ren;
        X_pos   = (rend{i}.data{1}-mnmn)/(mxmx-mnmn);
        msk = find(X_pos);
        ren(msk) = X_pos(msk)+(1+1.51/64);
        ax=axes('Parent',Fgraph,'units','normalized',...
            'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
            'Visible','off');
        image(ren*64,'Parent',ax);
        set(ax,'DataAspectRatio',[1 1 1], ...
            'PlotBoxAspectRatioMode','auto',...
            'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
    end
else
    % Combine the brain surface renderings with the blobs, and display using
    % 24 bit colour.
    %----------------------------------------------------------------------
    
    switch chromophore
        case 1
            c_s = 'HbO';
        case 2
            c_s =  'HbR';
        case 3
            c_s = 'HbT';
    end
    for h = sessions
        
        
        v_t = 0;
        if h == 0 
            h0 = 20;
            image_tab = ['Group_unc_' c_s];
        else
            h0 = h;
            activation_name = NIRS.Dt.fir.Sess(h).U(activation).name{1};
            if thres
                image_tab = ['Sess' int2str(h) '_EC_' c_s '_' activation_name];
            else
                image_tab = ['Sess' int2str(h) '_unc_' c_s '_' activation_name];
            end
            
        end

        
        for n = views,
            
            if h == 0
                activation_name = TOPO.v{n}.g.hb{chromophore}.c{activation}.c.name;
                [side_hemi spec_hemi] = nirs_get_brain_view(n);
                image_tab = [image_tab '_' spec_hemi activation_name];
            end

            %For old version of render file
            %Get the correct render image. Note that views for 1 and 2 are
            %inversed. Views for 3 and 4 are also inversed.
%             if n == 1
%                 i = 2;
%             elseif n == 2
%                 i = 1;
%             elseif n == 3
%                 i = 4;
%             elseif n == 4
%                 i = 3;
%             else
%                 i = n;
%             end
            %For new version of render file
            i = n;
            
            v_t = v_t + 1;
            n_t = v_t + 1;
            if n_t > length(views);
                n_t = n_t - length(views);
            end
            
            %ren = flipud(rend{i}.ren);
            ren = rend{i}.ren;
%             X = cell(3,1);
%             
%             for j=1:length(render_proj.contrast{h}.c_map{n}.data),
%                 mxmx = render_proj.contrast{h}.c_map{n}.mxmx{j};
%                 mnmn = render_proj.contrast{h}.c_map{n}.mnmn{j};
%                 X{j} = render_proj.contrast{h}.c_map{n}.data{j}/(mxmx-mnmn)-mnmn;%data(j), May need to modify if have multiple activations
%             end
%             for j=(length(render_proj.contrast{h}.c_map{n}.data)+1):3
%                 X{j}=zeros(size(X{1}));
%             end
% 
%             rgb = zeros([size(ren) 3]);
%             tmp = ren.*max(1-X{1}-X{2}-X{3},0);
%             for k = 1:3
%                 rgb(:,:,k) = tmp + X{1}*col_contrast(1,k) + X{2}*col_contrast(2,k) +X{3}*col_contrast(3,k);
%             end
%             rgb(rgb>1) = 1;  
            
            X_pos = cell(3,1);
            X_neg = cell(3,1);
            
            for j=1:length(render_proj.contrast{h0}.c_map{n}.data),
                mxmx = render_proj.contrast{h0}.c_map{n}.mxmx{j};
                mnmn = render_proj.contrast{h0}.c_map{n}.mnmn{j};
                X = render_proj.contrast{h0}.c_map{n}.data{j};
                X_pos_temp = X; X_neg_temp = X;
                X_pos_temp(X_pos_temp < 0) = 0;
                X_neg_temp(X_neg_temp > 0) = 0;
                X_neg_temp = abs(X_neg_temp);
                X_pos{j} = X_pos_temp/(mxmx-mnmn)-mnmn;%data(j), May need to modify if have multiple activations
                X_neg{j} = X_neg_temp/(mxmx-mnmn)-mnmn;%data(j), May need to modify if have multiple activations
            end
            for j=(length(render_proj.contrast{h0}.c_map{n}.data)+1):3
                X_pos{j}=zeros(size(X_pos{1}));
                X_neg{j}=zeros(size(X_neg{1}));
            end

            rgb_pos = zeros([size(ren) 3]);
            rgb_neg = zeros([size(ren) 3]);
            tmp_pos = ren.*max(1-X_pos{1}-X_pos{2}-X_pos{3},0);
            tmp_neg = ren.*max(1-X_neg{1}-X_neg{2}-X_neg{3},0);
            for k = 1:3
                rgb_pos(:,:,k) = tmp_pos + X_pos{1}*col_contrast(1,k) + X_pos{2}*col_contrast(2,k) +X_pos{3}*col_contrast(3,k);
            end
            for k = 1:3
                rgb_neg(:,:,k) = tmp_neg + X_neg{1}*col_contrast_neg(1,k) + X_neg{2}*col_contrast_neg(2,k) +X_neg{3}*col_contrast_neg(3,k);
            end
            rgb = (rgb_pos + rgb_neg)./2;
            %rgb = rgb_neg;
            rgb(rgb>1) = 1;  
            
            %add focus color
            data_focus = render_proj.all_focus{i}.map;
            for i0 = 1 : size(data_focus,1)
                for j0 = 1 : size(data_focus,2)
                    if (data_focus(i0,j0) > 0) && ((data_focus(i0-1, j0) == 0) || (data_focus(i0+1, j0) == 0) || (data_focus(i0, j0-1) == 0) || (data_focus(i0, j0+1) == 0))
                        rgb(i0,j0,2) = 1;%Green color
                        rgb(i0-1,j0,2) = 1;rgb(i0-2,j0,2) = 1;rgb(i0-3,j0,2) = 1;rgb(i0-4,j0,2) = 1;rgb(i0-5,j0,2) = 1;
                        rgb(i0+1,j0,2) = 1;rgb(i0+2,j0,2) = 1;rgb(i0+3,j0,2) = 1;rgb(i0+4,j0,2) = 1;rgb(i0-5,j0,2) = 1;
                        rgb(i0,j0-1,2) = 1;rgb(i0,j0-2,2) = 1;rgb(i0,j0-3,2) = 1;rgb(i0,j0-4,2) = 1;rgb(i0,j0-5,2) = 1;
                        rgb(i0,j0+1,2) = 1;rgb(i0,j0+2,2) = 1;rgb(i0,j0+3,2) = 1;rgb(i0,j0+4,2) = 1;rgb(i0,j0+5,2) = 1;
                    end
                end
            end
            
            %Display

            ax=axes('Parent',Fgraph,'units','normalized',...
                'Position',[rem(n_t-1,2)*0.5, floor((n_t-1)/2)*hght/nrow, 0.5, hght/nrow],...
                'nextplot','add', ...
                'Visible','off');
            image(rgb,'Parent',ax);
            set(ax,'DataAspectRatio',[1 1 1], ...
                'PlotBoxAspectRatioMode','auto',...
                'YTick',[],'XTick',[],...
                'XDir','normal','YDir','normal');
        end
        %Save projected image
        proj_image_location = [save_dir '\' image_tab];
        saveas(gcf, proj_image_location, 'tif');
    end
end




spm('Pointer','Arrow');
catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
end

%==========================================================================
% function surf_rend(dat,rend,col)
%==========================================================================
function surf_rend(dat,rend,col)

%-Setup figure and axis
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

ax0 = axes(...
    'Tag',      'SPMMeshRenderBackground',...
    'Parent',   Fgraph,...
    'Units',    'normalized',...
    'Color',    [1 1 1],...
    'XTick',    [],...
    'YTick',    [],...
    'Position', [-0.05, -0.05, 1.05, 0.555]);

ax = axes(...
    'Parent',   Fgraph,...
    'Units',    'normalized',...
    'Position', [0.05, 0.05, 0.9, 0.4],...
    'Visible',  'off');

H = spm_mesh_render('Disp',rend,struct('parent',ax));
spm_mesh_render('Overlay',H,dat,col);

try
    setAllowAxesRotate(H.rotate3d, setxor(findobj(Fgraph,'Type','axes'),ax), false);
end
    
%-Register with MIP
%--------------------------------------------------------------------------
try % meaningless when called outside spm_results_ui
    hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
    xyz  = spm_XYZreg('GetCoords',hReg);
    hs   = mydispcursor('Create',ax,dat.mat,xyz);
    spm_XYZreg('Add2Reg',hReg,hs,@mydispcursor);
end
  

%==========================================================================
function varargout = mydispcursor(varargin)

switch lower(varargin{1})
    %======================================================================
    case 'create'
    %======================================================================
    % hMe = mydispcursor('Create',ax,M,xyz)
    ax  = varargin{2};
    M   = varargin{3};
    xyz = varargin{4};
    
    [X,Y,Z] = sphere;
    vx = sqrt(sum(M(1:3,1:3).^2));
    X = X*vx(1) + xyz(1);
    Y = Y*vx(2) + xyz(2);
    Z = Z*vx(3) + xyz(3);
    hold(ax,'on');
    hs = surf(X,Y,Z,'parent',ax,...
        'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting', 'phong');
    set(hs,'UserData',xyz);
    
    varargout = {hs};
    
    %=======================================================================
    case 'setcoords'    % Set co-ordinates
    %=======================================================================
    % [xyz,d] = mydispcursor('SetCoords',xyz,hMe,hC)
    hMe  = varargin{3};
    pxyz = get(hMe,'UserData');
    xyz  = varargin{2};
    
    set(hMe,'XData',get(hMe,'XData') - pxyz(1) + xyz(1));
    set(hMe,'YData',get(hMe,'YData') - pxyz(2) + xyz(2));
    set(hMe,'ZData',get(hMe,'ZData') - pxyz(3) + xyz(3));
    set(hMe,'UserData',xyz);
    
    varargout = {xyz,[]};
    
    %=======================================================================
    otherwise
    %=======================================================================
    error('Unknown action string')

end
