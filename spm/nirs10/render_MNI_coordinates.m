function [rend, rendered_MNI] = render_MNI_coordinates(vx_wMNI,vx_MNI,wT1_info,VF,render_subj_info)
Nmark = size(vx_wMNI, 2);
if isfield(render_subj_info,'render_template')
    %use SPM template
    load([spm('dir') filesep 'rend' filesep 'render_single_subj.mat']);
    rnc = -1;
else
    if isfield(render_subj_info,'render_subject')
        try
            rf = render_subj_info.render_subject.render_file{1};
            rnc = render_subj_info.render_subject.render_normalize_choice;
            occipital_shift = render_subj_info.render_subject.occipital_shift;
            [dir0,fil0,ext0] = fileparts(rf);
            try
                %could be a rendered .mat file
                load(rf);
            catch
                %should be an image
                owd = pwd;
                cd(dir0);
                V0 = spm_vol(rf);
                %generate rendered image
                spm_surf(rf,1);
                file0 = fullfile(dir0,['render_' fil0 '.mat']);
                load(file0);
                %reorder the views, and save them
                rend0 = rend;
                rend{3} = rend0{4};
                rend{4} = rend0{3};
                %                 rend{2} = rend0{5};
                %                 rend{5} = rend0{2};
                rend{1} = rend0{2};
                rend{2} = rend0{1};
                save(file0,'rend');
                cd(owd);
            end
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            disp('Problem with render_subject');
        end
    end
end
%values: -1: standard treatment (normalization to template);
%1: normalization of subject; 0: no normalization
switch rnc
    case {-1,1}
        dat.mat = wT1_info.mat;
        dat.dim = wT1_info.dim;
        dat.XYZ = zeros(3,Nmark);
        for kk = 1:Nmark
            dat.t(1,(kk-1)+1:kk) = kk*10;
        end
        dat.XYZ(:,:) = vx_wMNI(1:3,:);
        tmp_mm_MNI = wT1_info.mat * vx_wMNI;
        dat.mm_MNI(:,:) = tmp_mm_MNI(1:3,:);
    case 0
        %specify .dat structure
        dat.mat = VF.mat;
        dat.dim = VF.dim;
        dat.XYZ = zeros(3,Nmark);
        for kk = 1:Nmark
            dat.t(1,(kk-1)+1:kk) = kk*10;
        end
        dat.XYZ(:,:) = vx_MNI(1:3,:);
        tmp_mm_MNI = VF.mat * vx_MNI;
        dat.mm_MNI(:,:) = tmp_mm_MNI(1:3,:);
end
clear tmp_MNI;
brt = NaN;

if (exist('rend','var') ~= 1), % Assume old format...
    rend = cell(size(Matrixes,1),1);
    for i=1:size(Matrixes,1),
        rend{i}=struct('M',eval(Matrixes(i,:)),...
            'ren',eval(Rens(i,:)),...
            'dep',eval(Depths(i,:)));
        rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
    end;
end;


for i=1:length(rend),
    rend{i}.max=0;
    rend{i}.data = cell(size(dat,1),1);
    if issparse(rend{i}.ren),
        % Assume that images have been DCT compressed
        % - the SPM99 distribution was originally too big.
        d = size(rend{i}.ren);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        rend{i}.ren = B1*rend{i}.ren*B2';
        % the depths did not compress so well with
        % a straight DCT - therefore it was modified slightly
        rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
    end;
    msk = rend{i}.ren>1;rend{i}.ren(msk)=1;
    msk = find(rend{i}.ren<0);rend{i}.ren(msk)=0;
end;

mx = zeros(length(rend),1)+eps;
mn = zeros(length(rend),1);

%for j=1:length(dat),
j=1;
XYZ = dat(j).XYZ;
mm_MNI = dat(j).mm_MNI;
t  = dat(j).t;
dim = dat(j).dim;
mat = dat(j).mat;

% transform from Talairach space to space of the rendered image
%---------------------------------------------------------------
for i=1:length(rend),
    M1  = rend{i}.M*mat;
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
    % calculate 'depth' of values
    %-------------------------------------------------------
    dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
    %         z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));
    switch rnc
        case {-1,1}
            switch i
                case 1
                    msk = find(mm_MNI(3,:) < -24 | mm_MNI(3,:) == -24); % ventral view
                case 2
                    msk = find(mm_MNI(3,:) > -24 | mm_MNI(3,:) == -24);  %dorsal view
                case 3
                    msk = find(mm_MNI(1,:) > 0 | mm_MNI(1,:) == 0); % lateral view (right hemisphere)
                case 4
                    msk = find(mm_MNI(1,:) < 0 | mm_MNI(1,:) == 0); % lateral view (left hemisphere)
                case 5
                    msk = find(mm_MNI(2,:) > -44 | mm_MNI(2,:) == -44); % frontal view
                case 6
                    msk = find(mm_MNI(2,:) < -44 | mm_MNI(2,:) == -44); % occipital view
            end
        case 0
            %this will depend on ...
            %                 vd = wT1_info.mat\ [0 0 24 1]';
            %                 vd = vd(1);
            %                 fo = wT1_info.mat\ [0 44 0 0];
            %                 fo = fo(1);
            vd = 24;
            fo = 44-occipital_shift; %44;
            switch i
                case 1
                    msk = find(mm_MNI(3,:) < -vd | mm_MNI(3,:) == -vd); % ventral view
                case 2
                    msk = find(mm_MNI(3,:) > -vd | mm_MNI(3,:) == -vd);  %dorsal view
                case 3
                    msk = find(mm_MNI(1,:) > 0 | mm_MNI(1,:) == 0); % lateral view (right hemisphere)
                case 4
                    msk = find(mm_MNI(1,:) < 0 | mm_MNI(1,:) == 0); % lateral view (left hemisphere)
                case 5
                    msk = find(mm_MNI(2,:) > -fo | mm_MNI(2,:) == -fo); % frontal view
                case 6
                    msk = find(mm_MNI(2,:) < -fo | mm_MNI(2,:) == -fo); % occipital view
            end
    end
    if ~isempty(msk),
        xyz =  xyz(:,msk);
        if ~isfinite(brt)
            t0  = t(msk);
        end;
        
        %%% Mahnoush 2012-01
        % Manually remove potentially duplicated* coordinates after rounding
        % *"Any elements of s in S = sparse(i,j,s,m,n,nzmax) that have duplicate values of i and j are added together."
        
        rxyz = round(xyz);
        [srxyz,ind_rxyz] = sortrows(rxyz');
        srxyz = srxyz'; ind_rxyz = ind_rxyz';
        dxyz(1,:) = diff(srxyz(1,:));
        dxyz(2,:) = diff(srxyz(2,:));
        Pdbl = find(dxyz(1,:)==0 & dxyz(2,:)==0);
        if ~isempty(Pdbl)
            Idbl = ind_rxyz(Pdbl);
            for numIdbl=1:size(Idbl,2)
                %         if dxyz(1,:)==0
                %             coord=1;
                %         else
                %             coord=2;
                %         end
                %         change the treshold of rounding to avoid duplicates:
                %         xyz(coord,Idbl(numIdbl)) = xyz(coord,Idbl(numIdbl)) - 0.5;
                % elimitate the duplicate and add at the end in rend.data
                t0(Idbl(numIdbl)) = 0;
                dupt0 = t(msk) - t0;
                % M "sparse: Any elements of s that are zero are ignored, along with the
                % corresponding values of i and j."
            end
            dupX0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), dupt0, d2(1), d2(2)));
            hld = 1; if ~isfinite(brt), hld = 0; end;
            dupX   = spm_slice_vol(dupX0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
            dupmsk = msk;
            dupmsk = find(dupX<0);
            dupX(dupmsk) = 0;
            Rdxyz{i}.data = rxyz(:,Idbl);
            Rdxyz{i}.ind = Idbl;
        else
            dupX = 0;
            % allocate
            Rdxyz{i}.data = 0;
            Rdxyz{i}.ind = 0;
        end
        
        clear dxyz rxyz srxyz ind_rxyz Pdbl Idbl
        
        %%% updated 2009-02-05
        
        %msk = find(xyz(1,:) > 0  & xyz(2,:) > 0 & xyz(3,:) > 0);
        %xyz = xyz(:,msk);
        %t0  = t0(:,msk);
        %%%
        X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
        hld = 1; if ~isfinite(brt), hld = 0; end;
        X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
        msk = find(X<0);
        X(msk) = 0;
        
        %%% Mahnoush
        % add omited duplicated points
        %         X = X + dupX;
        %%% --------
        
        % Add masking for statistical maps according to possible interpolation
        % Totally heuristic but seems to work well,
        % ISSUE: the 20 below works on the
        % MNI template, have to write in mm to be consistent
        % All this does is filter the current positions of the channel by a
        % circular filter, keeps all voxels over a threshold of 0.01 (since
        % we filter an image with max = 1 this hard coded number is not a
        % problem).
        view_mask_2d = double(X>0);
        hh=fspecial('disk',20);
        view_mask_2d=imfilter(view_mask_2d,hh,'same');
        view_mask_2d=view_mask_2d>0.01;
        % End change
        
    else
        X = zeros(size(rend{i}.dep));
        %%% Mahnoush
        dupX = zeros(size(rend{i}.dep));
        Rdxyz{i}.data = 0;
        Rdxyz{i}.ind = 0;
        %%% --------
        
        view_mask_2d = zeros(size(X));
    end;
    % Brighten the blobs
    if isfinite(brt), X = X.^brt; end;
    mx(j) = max([mx(j) max(max(X))]);
    mn(j) = min([mn(j) min(min(X))]);
    rend{i}.data{j} = X;
    rend{i}.ddata{j} = dupX; %Mahnoush
    % Adding to rendered_MNI since this is what is saved later in TopoData,
    % no need to change code elsewhere
    rendered_MNI{i}.view_mask_2d=view_mask_2d;
end;
%end;

rend{1}.mxmx = max(mx);
rend{1}.mnmn = min(mn);

for i=1:length(rend)
    temp_dep = rend{i}.dep;
    temp_ren = rend{i}.ren;
    temp_data = rend{i}.data{1};
    temp_ddata = rend{i}.ddata{1}; %Mahnoush
    temp_viewmask_2d=rendered_MNI{i}.view_mask_2d;
    rend{i}.dep = temp_dep(size(temp_dep,1):-1:1,:);
    rend{i}.ren = temp_ren(size(temp_ren,1):-1:1,:);
    rend{i}.data{1} = temp_data(size(temp_data,1):-1:1,:);
    rend{i}.ddata{1} = temp_ddata(size(temp_ddata,1):-1:1,:); %Mahnoush
    rendered_MNI{i}.view_mask_2d=temp_viewmask_2d(size(temp_data,1):-1:1,:);
end
for kk = 1:length(rend)
    for i = 1:Nmark
        [rows, cols, vals] = find(rend{kk}.data{1} == i*10);
        [drows, dcols, dvals] = find(rend{kk}.ddata{1} == i*10);
        if isempty(rows) == 1 || isempty(cols) == 1
            %%% -Mahnoush-
            if  ~isempty(find(Rdxyz{kk}.ind == i)) && ~(isempty(drows) == 1 || isempty(dcols) == 1)
                rendered_MNI{kk}.rchn(i,1) = round(mean(drows));
                rendered_MNI{kk}.cchn(i,1) = round(mean(dcols));
                %                 dupCn = Rdxyz{kk}.data(:,find(Rdxyz{kk}.ind == i));
                %                 rendered_MNI{kk}.rchn(i,1) = dupCn(1);
                %                 rendered_MNI{kk}.cchn(i,1) = dupCn(2) ;
            else
                %%% -----------
                rendered_MNI{kk}.rchn(i,1) = -1;
                rendered_MNI{kk}.cchn(i,1) = -1;
            end
            
        elseif isempty(rows) == 0 && isempty(cols) == 0
            rendered_MNI{kk}.rchn(i,1) = round(mean(rows));
            rendered_MNI{kk}.cchn(i,1) = round(mean(cols));
        end
    end
end