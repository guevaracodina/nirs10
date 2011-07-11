function out = nirs_criugm_getHelmet_autosave(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2011-06

staxp = job.subj.helmet.staxp;
sDtp = job.subj.sDtp;

load(fullfile(fileparts(which('nirs10')),'nirs10_templates','fid_names.mat'));
names_na = fid_names.names_na;
names_le = fid_names.names_le;
names_re = fid_names.names_re;

%% case initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Brainsight text file and build NIRS matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobRB.staxp = staxp;
jobRB.sDtp = sDtp;
handles = nirs_criugm_readbrainsight(jobRB);

% raw lists of points from stereotaxic logiciel (here Brainsight)
handles.On = handles.Cf.H.O.n';% *** distinguer les deux noms si modification de la liste de gauche (par exemple si on retire les points déjà selectionnés)
handles.Sn = handles.Cf.H.S.n';
handles.Dn = handles.Cf.H.D.n';
% coordinates of these points
handles.Op_rom = handles.Cf.H.O.r.o.mm.p';
handles.Sp_rom = handles.Cf.H.S.r.o.mm.p';
handles.Dp_rom = handles.Cf.H.D.r.o.mm.p';

% extract Nasion, LeftEar, RightEar
icn=1;ind_na = 0;
while ind_na==0
ind_na = sum(strcmpi(handles.Cf.H.O.n',names_na{1,icn}));
icn =icn+1;
end
ind_na = strcmpi(handles.Cf.H.O.n',names_na{1,icn-1});
ind_na_n = sum(ind_na.*(1:length(handles.Cf.H.O.n'))');

icL=1;ind_le = 0;
while ind_le==0
ind_le = sum(strcmpi(handles.Cf.H.O.n',names_le{1,icL}));
icL =icL+1;
end
ind_le = strcmpi(handles.Cf.H.O.n',names_le{1,icL-1});
ind_le_n = sum(ind_le.*(1:length(handles.Cf.H.O.n'))');

icR=1;ind_re = 0;
while ind_re==0
ind_re = sum(strcmpi(handles.Cf.H.O.n',names_re{1,icR}));
icR =icR+1;
end
ind_re = strcmpi(handles.Cf.H.O.n',names_re{1,icR-1});
ind_re_n = sum(ind_re.*(1:length(handles.Cf.H.O.n'))');

%pre-selected points
handles.Fn_s{1,1} = handles.Cf.H.O.n{1,ind_na_n};
handles.Fn_s{2,1} = handles.Cf.H.O.n{1,ind_le_n};
handles.Fn_s{3,1} = handles.Cf.H.O.n{1,ind_re_n};

handles.Sn_s = handles.Cf.H.S.n';
handles.Dn_s = handles.Cf.H.D.n';
handles.Qn_s = {};

% coordinates of these pre-selected points
handles.Fp_s_rom{1,1} = handles.Cf.H.O.r.o.mm.p{1,ind_na};
handles.Fp_s_rom{1,2} = handles.Cf.H.O.r.o.mm.p{2,ind_na};
handles.Fp_s_rom{1,3} = handles.Cf.H.O.r.o.mm.p{3,ind_na};

handles.Fp_s_rom{2,1} = handles.Cf.H.O.r.o.mm.p{1,ind_le};
handles.Fp_s_rom{2,2} = handles.Cf.H.O.r.o.mm.p{2,ind_le};
handles.Fp_s_rom{2,3} = handles.Cf.H.O.r.o.mm.p{3,ind_le};

handles.Fp_s_rom{3,1} = handles.Cf.H.O.r.o.mm.p{1,ind_re};
handles.Fp_s_rom{3,2} = handles.Cf.H.O.r.o.mm.p{2,ind_re};
handles.Fp_s_rom{3,3} = handles.Cf.H.O.r.o.mm.p{3,ind_re};

handles.Sp_s_rom = handles.Cf.H.S.r.o.mm.p';
handles.Dp_s_rom = handles.Cf.H.D.r.o.mm.p';
handles.Qp_s_rom = {};

%% case save
%NIRS
NIRS =[];
NIRS.Cf.H.F.n = handles.Fn_s';    % Fiducials
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

save(fullfile(sDtp,'NIRS_Cf.mat'),'NIRS');
out = 1;
end