function HDM = nirs_initial_hdm_setup(O,DO,EM,IC,S,dt,dir1)
HDM = [];
%Various options
HDM.O = O;
HDM.DO = DO;
HDM.EM = EM;
HDM.IC = IC;
HDM.dt = dt;
HDM.dir1 = dir1;
%Simulation option
HDM.S = S;