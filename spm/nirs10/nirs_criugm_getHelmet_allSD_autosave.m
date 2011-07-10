function out = nirs_criugm_getHelmet_allSD_autosave(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2011-06

staxp = job.subj.helmet.staxp;
sDtp = job.subj.sDtp;

%% case initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Brainsight text file and build NIRS matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobRB.staxp = staxp;
jobRB.sDtp = sDtp;
handles = nirs_criugm_readbrainsight(jobRB);

% raw lists of points from stereotaxic logiciel (here Brainsight)
handles.On = handles.Cf.H.O.n';% ***  NIRS
handles.Sn = handles.Cf.H.S.n';%NIRS
handles.Dn = handles.Cf.H.D.n';%NIRS
% *** distinguer les deux noms si modification de la liste de gauche
% (par exemple si on retire les points déjà selectionnés)
% coordinates of these points
handles.Op_rom = handles.Cf.H.O.r.o.mm.p';%NIRS
handles.Sp_rom = handles.Cf.H.S.r.o.mm.p';%NIRS
handles.Dp_rom = handles.Cf.H.D.r.o.mm.p';%NIRS

%%% extract Nasion, LeftEar, RightEar
ind_nasion = strcmp(handles.Cf.H.O.n','Nasion');
ind_nasion_n = sum(ind_nasion.*(1:length(handles.Cf.H.O.n'))');
ind_AL = strcmp(handles.Cf.H.O.n','LeftEar');
ind_AL_n = sum(ind_AL.*(1:length(handles.Cf.H.O.n'))');
ind_AR = strcmp(handles.Cf.H.O.n','RightEar');
ind_AR_n = sum(ind_AR.*(1:length(handles.Cf.H.O.n'))');

%pre-selected points
handles.Fn_s{1,1} = handles.Cf.H.O.n{1,ind_nasion_n};
try
handles.Fn_s{2,1} = handles.Cf.H.O.n{1,ind_AL_n};
catch
   handles.Fn_s{2,1} = handles.Cf.H.O.n{1,3}; 
end
try
handles.Fn_s{3,1} = handles.Cf.H.O.n{1,ind_AR_n};
catch
handles.Fn_s{3,1} = handles.Cf.H.O.n{1,4};    
end
handles.Sn_s = handles.Cf.H.S.n';
handles.Dn_s = handles.Cf.H.D.n';
handles.Qn_s = {};

% coordinates of these pre-selected points
try
handles.Fp_s_rom{1,1} = handles.Cf.H.O.r.o.mm.p{1,ind_nasion};
handles.Fp_s_rom{1,2} = handles.Cf.H.O.r.o.mm.p{2,ind_nasion};
handles.Fp_s_rom{1,3} = handles.Cf.H.O.r.o.mm.p{3,ind_nasion};

handles.Fp_s_rom{2,1} = handles.Cf.H.O.r.o.mm.p{1,ind_AL};
handles.Fp_s_rom{2,2} = handles.Cf.H.O.r.o.mm.p{2,ind_AL};
handles.Fp_s_rom{2,3} = handles.Cf.H.O.r.o.mm.p{3,ind_AL};

handles.Fp_s_rom{3,1} = handles.Cf.H.O.r.o.mm.p{1,ind_AR};
handles.Fp_s_rom{3,2} = handles.Cf.H.O.r.o.mm.p{2,ind_AR};
handles.Fp_s_rom{3,3} = handles.Cf.H.O.r.o.mm.p{3,ind_AR};

handles.Sp_s_rom = handles.Cf.H.S.r.o.mm.p';
handles.Dp_s_rom = handles.Cf.H.D.r.o.mm.p';
handles.Qp_s_rom = {};

%% case save
% sort Fiducials so that they respect order used in run_coreg
%nasion
ind_nasion = strcmp(handles.Fn_s,'Nasion');
ind_nasion_n = sum(ind_nasion.*(1:length(handles.Fn_s))');
Fn_sorted_s{1,1} = handles.Fn_s(ind_nasion_n,1);
Fp_sorted_s_rom(1,:) = handles.Fp_s_rom(ind_nasion,:);
%Left Ear
ind_AL = strcmp(handles.Fn_s,'LeftEar');
ind_AL_n = sum(ind_AL.*(1:length(handles.Fn_s))');
Fn_sorted_s{2,1} = handles.Fn_s(ind_AL_n,1);
Fp_sorted_s_rom(2,:) = handles.Fp_s_rom(ind_AL,:);
%Right Ear
ind_AR = strcmp(handles.Fn_s,'RightEar');
ind_AR_n = sum(ind_AR.*(1:length(handles.Fn_s))');
Fn_sorted_s{3,1} = handles.Fn_s(ind_AR_n,1);
Fp_sorted_s_rom(3,:) = handles.Fp_s_rom(ind_AR,:);
end
%NIRS
NIRS =[];
try
NIRS.Cf.H.F.n = Fn_sorted_s';    % Fiducials
NIRS.Cf.H.S.n = handles.Sn_s';   % Sources
NIRS.Cf.H.D.n = handles.Dn_s';   % Detectors
NIRS.Cf.H.Q.n = handles.Qn_s';   % Points of interest

NIRS.Cf.H.F.N = size(Fn_sorted_s,1);
NIRS.Cf.H.S.N = size(handles.Sn_s,1);
NIRS.Cf.H.D.N = size(handles.Dn_s,1);
NIRS.Cf.H.Q.N = size(handles.Qn_s,1);

NIRS.Cf.H.F.r.o.mm.p = cell2mat(Fp_sorted_s_rom');
NIRS.Cf.H.S.r.o.mm.p = cell2mat(handles.Sp_s_rom');
NIRS.Cf.H.D.r.o.mm.p = cell2mat(handles.Dp_s_rom');
NIRS.Cf.H.Q.r.o.mm.p = cell2mat(handles.Qp_s_rom');
catch
    disp('Problem with NIRS_Cf');
end
save(fullfile(sDtp,'NIRS_Cf.mat'),'NIRS');
out = 1;
end
