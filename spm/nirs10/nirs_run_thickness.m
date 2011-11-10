function out = nirs_run_thickness(job)

prec=1;

% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    % Load NIRS.mat
    try
        NIRS = [];
        load(job.NIRSmat{iSubj,1});
        
        %segmented T1
        V = spm_vol(NIRS.Dt.ana.T1seg);
        Y = spm_read_vols(V);
        
        max_size = V.dim';
        min_size = [1;1;1];
        
        %fitted positions on skin and on cortex
        Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
        Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
        Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
        
        % SP sampling points
        SPp_rmm = Pp_rmm(job.SP);
        SPfp_rmm = Pfp_rmm(job.SP);
        SPp_c1_rmm = Pp_c1_rmm(job.SP);
        % - directions : pointing outside of the head
        SPd_rmm = -(SPp_rmm - SPp_c1_rmm)/norm(SPp_rmm - SPp_c1_rmm);
        
        % SP will contain thicknesses
        SP=zeros(5,size(job.SP));
        for Pi = 1:size(job.SP)
            SPi = zeros(3,5);
            
            SPi(:,1) = SPp_c1_rmm(:,Pi);%c1
            SPi1_v = V.mat\[SPi(:,1);1];
            SPi1_v = round(SPi1_v(1:3));
            %initialization : projection on WM
            while Y(SPi1_v(1),SPi1_v(2),SPi1_v(3))==1
                SPi1_m = SPi1_m + prec*SPd_rmm(:,Pi);
                
                SPi1_v = V.mat\[SPi1_m;1];
                SPi1_v = round(SPi1_v(1:3));
                SPi1_v = min(max(SPi1_v,min_size),max_size);
            end
            SPi(:,2) = SPi1_m - prec*SPd_rmm(:,Pi);
            SPi2_v = V.mat\[SPi(:,2);1];
            SPi2_v = round(SPi2_v(1:3));
            
            %current position
            Ki_v = min(max(SPi2_v,min_size),max_size);
            Ki_m = SPi(:,2);
            
            for j=[1 3 4 5] % juste sur les images sans masque
                while Y(Ki_v(1),Ki_v(2),Ki_v(3))==j
                    Ki_m = Ki_m - prec*SPd_rmm(:,Pi);
                    
                    Ki_v = V.mat\[Ki_m;1];
                    Ki_v = round(Ki_v(1:3));
                    Ki_v = min(max(Ki_v,min_size),max_size);
                end
                SPi(:,j) = Ki_m + prec*SPd_rmm(:,Pi);
                
                switch j % computing thicknesses
                    case 1
                        SP(1,Pi) = norm(SPi(:,2) - SPi(:,1));
                    case 2
                        SP(2,Pi) = 10;
                    case 3
                        SP(3,Pi) = norm(SPi(:,3) - SPi(:,1));
                    otherwise
                        SP(j,Pi) = norm(SPi(:,j-1) - SPi(:,j));
                end
            end
        end
        
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Estimate layer thickness failed for the ' int2str(iSubj) 'th subject.']);
    end
end
out.NIRSmat = job.NIRSmat;
end