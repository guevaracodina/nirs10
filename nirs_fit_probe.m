function out = nirs_fit_probe(job)

P = job.image_in{:};
V = spm_vol(P);
Y = spm_read_vols(V);

Pp_rmm = job.Pp_rmm;
Pp_c1_rmm = job.Pp_c1_rmm;
NP = job.NP;

Pd_rmm = -(Pp_rmm - Pp_c1_rmm)/norm(Pp_rmm - Pp_c1_rmm);% pointing outside of the head
Pfp_rmm = Pp_c1_rmm;

%%%%%
% method inspired from article: Qianqian Fang and David A. Boas, "Monte
% Carlo simulation of photon migration in 3D turbid media accelerated by
% graphics processing units," Opt. Express 17, 20178-20190 (2009)

%1) find the first intersection interface (direction of propagation must be known)
%2) neighbouring vowel is an air voxel or not
%3) C3
%%%%%

max_size = V.dim';
min_size = [1;1;1];

prec=1;

switch job.lby
    case 'coreg'
        for Pi = 1:NP
            Pfp_rmv_i = V.mat\[Pfp_rmm(:,Pi);1];
            Pfp_rmv_i = round(Pfp_rmv_i(1:3));
            Pfp_rmv_i = min(max(Pfp_rmv_i,min_size),max_size);
            
            while Y(Pfp_rmv_i(1),Pfp_rmv_i(2),Pfp_rmv_i(3))>0
                Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) - prec*Pd_rmm(:,Pi);
                
                Pfp_rmv_i = V.mat\[Pfp_rmm(:,Pi);1];
                Pfp_rmv_i = round(Pfp_rmv_i(1:3));
                Pfp_rmv_i = min(max(Pfp_rmv_i,min_size),max_size);
            end
            Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) + prec*Pd_rmm(:,Pi);
        end
        
    case 'configMC_MCX'% after resizing, some S might be inside the volume
        Pfp_rmiv = round(job.Pfp_ancienne_rmiv);
        
        for Pi = 1:NP
            Pfp_rmiv_i = Pfp_rmiv(:,Pi);
            Pfp_rmm_i =  V.mat*[Pfp_rmiv_i;1];
            
            while Y(Pfp_rmiv_i(1),Pfp_rmiv_i(2),Pfp_rmiv_i(3))~=0
%                 if size(Pfp_rmiv_i,1)==3
%                     Pfp_rmm_i =  V.mat*[Pfp_rmiv_i;1];
%                 else
%                     Pfp_rmm_i =  V.mat*Pfp_rmiv_i;
%                 end
                Pfp_rmm_i = Pfp_rmm_i - prec*[Pd_rmm(:,Pi);0];
                Pfp_rmiv_i = round(V.mat\Pfp_rmm_i);
            end
            Pfp_rmm(:,Pi) = Pfp_rmiv_i(1:3,:);
        end
        
        case 'configMC_tMC'% after resizing, some S might be inside the volume
        Pfp_rmiv = round(job.Pfp_ancienne_rmiv);
        
        for Pi = 1:NP
            Pfp_rmiv_i = Pfp_rmiv(:,Pi);
            Pfp_rmm_i =  V.mat*[Pfp_rmiv_i;1];
            
            while Y(Pfp_rmiv_i(1),Pfp_rmiv_i(2),Pfp_rmiv_i(3))==0
%                 if size(Pfp_rmiv_i,1)==3
%                     Pfp_rmm_i =  V.mat*[Pfp_rmiv_i;1];
%                 else
%                     Pfp_rmm_i =  V.mat*Pfp_rmiv_i;
%                 end
                Pfp_rmm_i = Pfp_rmm_i + prec*[Pd_rmm(:,Pi);0];
                Pfp_rmiv_i = round(V.mat\Pfp_rmm_i);
            end
            Pfp_rmm(:,Pi) = Pfp_rmiv_i(1:3,:);
        end
end

out{1} = Pfp_rmm;
end
%previous version
%     while sum(Y(test_vx(1,1),test_vx(2,1),test_vx(3,1)))>0%==0
%         Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) - prec*Pd_rmm(:,Pi);
%         %     while sum(Y(test_vx(1,1),test_vx(2,1),test_vx(3,1)))>0%==0
%         %         Pfp_rmm(:,Pi) = Pp_c1_rmm(:,Pi) + count*prec*Pd_rmm(:,Pi);
%
%         test_MNI = Pfp_rmm(:,Pi);
%         test_vx = V.mat\[test_MNI;1];
%         test_vx = round(test_vx(1:3));
%         test_vx = min(max(test_vx,min_size),max_size);
%     end
%     Pfp_rmm(:,Pi) = Pfp_rmm(:,Pi) + 20*prec*Pd_rmm(:,Pi);

%     %%% pour respecter les differentes orientations on a besoin de
%     %%% travailler sur les voxels et non les mm
%     while sum(Y(test_vx(1,1),test_vx(2,1),test_vx(3,1)))>0%==0
%         test_rem = test_vx;
%         test_vx = test_vx - 10*prec*Pd_rmv(:,Pi);
%          test_vx = round(test_vx(1:3));
%         test_vx = min(max(test_vx,min_size),max_size);
%     end
%     Pfp2_rmv(:,Pi) =  test_rem;