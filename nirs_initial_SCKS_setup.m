function SCKS = nirs_initial_SCKS_setup(O,DO,IC,S,dt,dir1)
SCKS = [];
%Various options
SCKS.O = O;
SCKS.DO = DO;
SCKS.IC = IC;
SCKS.O.include_HbR = IC.include_HbR; %for nirs_gx
SCKS.O.include_HbO = IC.include_HbO;
SCKS.O.include_HbT = IC.include_HbT;
SCKS.dt = dt;
SCKS.dir1 = dir1;
%Simulation option
SCKS.S = S;