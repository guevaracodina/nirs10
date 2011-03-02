function out = nirs_fit_probe(job)

%prec = 0.001;
prec=1;

P = job.image_in{:};
V = spm_vol(P);
Y = spm_read_vols(V);

% NIRS = job.NIRS;
%
% Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
% Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
% NP = NIRS.Cf.H.P.N;
%%%%% a changer eventuellement
Pp_rmm = job.Pp_rmm;
Pp_c1_rmm = job.Pp_c1_rmm;
NP = job.NP;
%%%%%

Pd_rmm = -(Pp_rmm - Pp_c1_rmm)/norm(Pp_rmm - Pp_c1_rmm);% pointing outside of the head
Pfp_rmm = Pp_c1_rmm;

% method inspired from article: Qianqian Fang and David A. Boas, "Monte
% Carlo simulation of photon migration in 3D turbid media accelerated by
% graphics processing units," Opt. Express 17, 20178-20190 (2009)

%1) find the first intersection interface (direction of propagation must be known)
%2) neighbouring vowel is an air voxel or not
%3) C3

% handmade method
count=1;
max_size = V.dim';
min_size = [1;1;1];

for Pi = 1:NP
    test_MNI = Pfp_rmm(:,Pi);
    test_vx = V.mat\[test_MNI;1];
    test_vx = round(test_vx(1:3));
    
    test_vx = min(max(test_vx,min_size),max_size);
    
    %previous version
    while sum(Y(test_vx(1,1),test_vx(2,1),test_vx(3,1)))>0%==0
        Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) - prec*Pd_rmm(:,Pi);
        %     while sum(Y(test_vx(1,1),test_vx(2,1),test_vx(3,1)))>0%==0
        %         Pfp_rmm(:,Pi) = Pp_c1_rmm(:,Pi) + count*prec*Pd_rmm(:,Pi);
        
        test_MNI = Pfp_rmm(:,Pi);
        test_vx = V.mat\[test_MNI;1];
        test_vx = round(test_vx(1:3));
        test_vx = min(max(test_vx,min_size),max_size);
    end
    Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) + 20*prec*Pd_rmm(:,Pi);
end

out{1} = Pfp_rmm;
end