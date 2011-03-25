function out = nirs_run_view3d(job)
% View optodes en 3D
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% based on Michèle Desjardins version of prepMCsim2
%
% Clément Bonnéry
% 2010-06

image_in = job.image_in;
V = spm_vol(image_in{:});
Y = spm_read_vols(V);

NIRSmat = job.NIRSmatSingle;
load(NIRSmat{:});

if isfield(NIRS.Cf.H.P.r.m.mm,'p')
    Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;% positions on head surface
    Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;% positions on cortex
    Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;% positions fitted on skin
end

NP = NIRS.Cf.H.P.N;
NS = NIRS.Cf.H.S.N;
ND = NIRS.Cf.H.D.N;

% % Fitted positions (positions are fitted with respect to the scalp)
% jobe.NIRS = NIRS;
% jobe.image_in = image_in;
% out = nirs_fit_probe(jobe);
% Pfp_rmm = out{1};

% from MNI real space (mm) to MNI voxel space
for i=1:NP
    Pp_rmv(:,i) = V.mat\[Pp_rmm(:,i);1];
    Pp_c1_rmv(:,i) = V.mat\[Pp_c1_rmm(:,i);1];
    Pfp_rmv(:,i) = V.mat\[Pfp_rmm(:,i);1];
end
% return to 3D space (V.mat is (:,4))
Pp_rmv = Pp_rmv(1:3,:);
Pp_c1_rmv = Pp_c1_rmv(1:3,:);

%% Display 3D image of the ROI with optodes pointing towards their
% respective directions, allowing a final check of the setup
% Surface de l'IRM
set(0,'defaultfigurevisible','on');
hfig_all3D = figure;

patch(isosurface(smooth3(Y),4),...
    'FaceColor',[1,.75,.65],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);

% Views adjustments
view(90,-90)
daspect([1,1,1])
lightangle(70,-70);
set(hfig_all3D,'Renderer','zbuffer');
lighting phong

% Display optodes and their directions
hold on

% adding tags
for Pi = 1:NP
    if Pi<=NS
        list{1,1} = 'r';
        list{2,1} = 10;
        list{3,1} = 'S#';
    elseif (Pi>NS && Pi<=NS+ND)
        list{1,1} = 'b';
        list{2,1} = 2;
        list{3,1} = 'D#';
    elseif Pi>NS+ND
        list{1,1} = 'g';
        list{2,1} = 2;
        list{3,1} = 'Q#';
    end
    
    xp = Pp_rmv(2,Pi);
    yp = Pp_rmv(1,Pi);
    zp = Pp_rmv(3,Pi);
    text(xp,yp,zp,[list{3,1} int2str(Pi)],'FontWeight','bold','Color',list{1,1});
    
    xfp = Pfp_rmv(2,Pi);
    yfp = Pfp_rmv(1,Pi);
    zfp = Pfp_rmv(3,Pi);
    text(xfp,yfp,zfp,'X','FontWeight','bold','Color','black');%list{1,1}
    
    plot3([Pp_c1_rmv(2,Pi),Pp_rmv(2,Pi)],...
        [Pp_c1_rmv(1,Pi),Pp_rmv(1,Pi)],...
        [Pp_c1_rmv(3,Pi),Pp_rmv(3,Pi)], 'Linewidth',list{2,1},'Color',list{1,1});
end
hold off

clear NIRS
load(NIRSmat{:});
NIRS.Cf.H.P.r.m.mm.fp =Pfp_rmm;
NIRS.Cf.H.P.r.m.vx.fp = Pfp_rmv(1:3,:);
save(NIRSmat{:},'NIRS');

out =1;
return