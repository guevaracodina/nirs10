function [varargout] = nirs_buildroi(varargin)
% Set boundaries and build ROI
% FORMAT out_buildroi = spm_mc_buildroi(V,b,'CmdLine','vertice');
% V       - a vector of structures containing image volume information.
% b       - boundaries
% bb      - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
% CmdLine - choice to get vertice or build roi once all the vertices have
% been set
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% $Id: spm_mc_buildroi.m 2010-06-11 Clément $
% Mise au format par Clément : 2010-11-17 et 02-2011

V = varargin{1};
Y = spm_read_vols(V);
[dir,name] = fileparts(V.fname);

b = varargin{2};
Action = varargin{3};

if strcmp(Action,'getvertice')
    vertice = varargin{4};
elseif strcmp(Action,'build_ROI')
    output_prefix = varargin{4};
    NIRSmat = varargin{5};
    Ckpt = varargin{6};
end

switch Action
    
    case 'getvertice'
        read_vertice = spm_orthviews('Pos',1);
        switch vertice
            case 'rf'
                b(2,3) = read_vertice(3,1);
                b(1,1) = read_vertice(1,1);
            case 'up'
                b(2,2) = read_vertice(2,1);
            case 'lb'
                b(2,1) = read_vertice(1,1);
                b(1,3) = read_vertice(3,1);
            case 'lo'
                b(1,2) = read_vertice(2,1);
        end
        varargout{1} = b;
        clear b
        
    case 'build_ROI'
        disp('Running ''Build ROI''');
        bb = round(b);
        % What is done :
        % - each volume must be able to live on its own (each must have an
        % origin a V.mat matrix and all must be standardised so that it is
        % an image as any other)
        % - both transformations inbetween origin voxel and origin point
        % (in mm) must be saved in the NIRS matrix
        
        % BIEN COMPRENDRE COMMENT CA SE FAIT...
        % the smallest value must always be in the first raw
        for i_col =1:3
            if bb(1,i_col)>bb(2,i_col)
                tempvar = bb(2,i_col);
                bb(2,i_col) = bb(1,i_col);
                bb(1,i_col) = tempvar;
            end
        end
        
        %%%%%%%
        % roi must contain the selected channels
        load(NIRSmat);
        %%% on sauve les fibres qui sont a la surface dans la ROI
        NS = NIRS.Cf.H.S.N;
        Pfp_rmv = NIRS.Cf.H.P.r.m.vx.fp;
        Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
        Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
        Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
        Cid = NIRS.Cf.H.C.id;
        NC = NIRS.Cf.H.C.N;
        wl = NIRS.Cf.dev.wl;
        
        nc = NC/length(wl);
        Ckpt = Ckpt.keepChannels;
        for iwl=2:length(wl)
            Ckpt = [Ckpt Ckpt+(iwl-1)*nc];
        end
        Skpt = unique(Cid(2,Ckpt));
        Dkpt = unique(Cid(3,Ckpt));
        Pkpt=[Skpt Dkpt+NS];

        for i=1:size(Pkpt,2)
            Pfp_roi_rmm(:,i) = Pfp_rmm(:,Pkpt(i));
            Pfp_roi_rmv(:,i) = Pfp_rmv(:,Pkpt(i));
            Pp_roi_rmm(:,i) = Pp_rmm(:,Pkpt(i));
            Pp_roi_c1_rmm(:,i) = Pp_c1_rmm(:,Pkpt(i));
        end
        %         Pn_segR_fp_rmv =[];
        marge =20;
        bb(1,1) = min(min(Pfp_roi_rmv(1,:))-marge,bb(1,1));
        bb(1,2) = min(min(Pfp_roi_rmv(2,:))-marge,bb(1,2));
        bb(2,1) = max(max(Pfp_roi_rmv(1,:))+marge,bb(2,1));
        bb(2,2) = max(max(Pfp_roi_rmv(2,:))+marge,bb(2,2));
        bb(1,3) = min(min(Pfp_roi_rmv(3,:))-marge,bb(1,3));
        bb(2,3) = max(max(Pfp_roi_rmv(3,:))+marge,bb(2,3));
        bb = round(bb);
        
        % the size of the plotted image can be bigger than the size read in
        % the header, in such case the value kept is the one of header
        for ib =1:3
            bb(1,ib) = max(1,bb(1,ib));
            bb(2,ib) = min(V.dim(ib),bb(2,ib));
        end
        
        % With resizing
        %  mat_roi2raw is the transformation matrix from the voxel space
        %  of the ROI to the voxel space of the raw image
        mat_roi2raw = eye(4) + [zeros(4,3) [bb(1,:)';0]];
        
        % Sources and detectors positions must be updated also
        for i=1:size(Pkpt,2)
            PfpR_roi_rmv(:,i) = [Pfp_roi_rmv(:,i);1] - mat_roi2raw(:,4);
        end
        %then saved
        if ~isfield(NIRS,'Cs'), NIRS.Cs={}; end
        if isfield(NIRS.Cs,'temp'), clear NIRS.Cs.temp; end
        NIRS.Cs.temp.Pkpt = Pkpt;
        NIRS.Cs.temp.NSkpt = size(Skpt,2);
        NIRS.Cs.temp.NDkpt = size(Dkpt,2);
        NIRS.Cs.temp.Pfp_roi_rmv = PfpR_roi_rmv(1:3,:);
        NIRS.Cs.temp.Pfp_roi_rmm = Pfp_roi_rmm;
        NIRS.Cs.temp.Pp_roi_rmm = Pp_roi_rmm;
        NIRS.Cs.temp.Pp_roi_c1_rmm = Pp_roi_c1_rmm;
        NIRS.Cs.temp.segR = fullfile(dir,[output_prefix,name,'.nii']);
        save(NIRSmat,'NIRS');
        
        V_roi.dim = [bb(2,1)-bb(1,1)+1 bb(2,2)-bb(1,2)+1 bb(2,3)-bb(1,3)+1];
        % V_roi.mat is the mat from voxel ROI to the raw image space in mm
        V_roi.mat = V.mat*mat_roi2raw;
        
        V_roi = struct('fname',fullfile(dir,[output_prefix,name,'.nii']),...
            'dim',  V_roi.dim,...
            'dt',   V.dt,...
            'pinfo',V.pinfo,...
            'mat',V_roi.mat);
        
        Y_roi = Y(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3));
        
        V_roi = spm_create_vol(V_roi);
        V_roi = spm_write_vol(V_roi, Y_roi);
        
        %         % rajout pour que tout roule...
        %         NIRS = [];
        %         load(NIRSmat);
        %         %%% on sauve les fibres qui sont a la surface dans la ROI
        %         NP = NIRS.Cf.H.P.N;
        %         Pfp_rmv = NIRS.Cf.H.P.r.m.vx.fp;
        %         %         Pfp_segR_rmv = [];
        %         Pn_segR_fp_rmv =[];
        %         for i=1:NP
        %             if (Pfp_rmv(1,i)>b(1,1) && Pfp_rmv(1,i)<b(2,1) &&...
        %                     Pfp_rmv(2,i)>b(1,2) && Pfp_rmv(2,i)<b(2,2) &&...
        %                     Pfp_rmv(3,i)>b(1,3) && Pfp_rmv(3,i)<b(2,3))
        %                 %                 Pfp_segR_rmv = [Pfp_segR_rmv Pfp_rmv(:,i)] ;
        %                 Pn_segR_fp_rmv = [Pn_segR_fp_rmv i];
        %             end
        %         end
        
        %         if isfield(NIRS.Cs,'temp'), clear NIRS.Cs.temp; end
        %         NIRS.Cs.temp.P_segR.n = Pn_segR_fp_rmv;
        %         %%% il va falloir avoir le choix pour la simulation qu'on fait
        %         %%% rouler ... donc NIRS.Cs.mcs(imcs).segR
        %         NIRS.Cs.temp.segR = fullfile(dir,[output_prefix,name,'.nii']);
        %         save(NIRSmat,'NIRS');
        
        varargout{1} = V_roi;
        disp('Done    ''Build ROI''');
        disp('Done');
end
%         if c==0 % no crop_image
%             % No resizing
%             Y_roi = zeros(V.dim(1:3));
%             Y_roi(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)) = Y(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3));
%
%             V_roi = struct('fname',fullfile(dir,[output_prefix,name,'.nii']),...
%                 'dim',  V.dim,...
%                 'dt',   V.dt,...
%                 'pinfo',V.pinfo,...
%                 'mat',V.mat);
%         else