function out = nirs_buildroi_SDcorrected(job)
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
% 02-2011


% ce module opere a partir d'une premiere ROI et des channels demandes, il
% faut connaitre le cote des channels......................................

load(job.NIRSmat{:});

Ckpt = job.Ckpt;
Cid = NIRS.Cf.H.C.id;
NC = NIRS.Cf.H.C.N;

wl = NIRS.Cf.dev.wl;
nc = NC/length(wl);
for iwl=2:length(wl)
    Ckpt = [Ckpt Ckpt+(iwl-1)*nc];
end

Skpt = unique(Cid(2,Ckpt));
Dkpt = unique(Cid(3,Ckpt));
NS = NIRS.Cf.H.S.N;
Pkpt=[Skpt Dkpt+NS];
Pp_rmm =[];
Pfp_rmv =[];
Pp_c1_rmm =[];
for i=1:length(Pkpt)
    Pp_rmm = [Pp_rmm NIRS.Cf.H.P.r.m.mm.p(:,Pkpt(i))];
    Pfp_rmv = [Pfp_rmv NIRS.Cf.H.P.r.m.vx.fp(:,Pkpt(i))];
    Pp_c1_rmm = [Pp_c1_rmm NIRS.Cf.H.P.r.m.mm.c1.p(:,Pkpt(i))];
end

bb = job.bb;
max(Pfp_rmv(1,:))
max(Pfp_rmv(2,:))
max(Pfp_rmv(3,:))


end