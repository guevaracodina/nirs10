function out = do_fit(varargin)
prec=1;

Y = varargin{2};
Pfp_rmiv_i = varargin{3};
Pd_rmiv_i = varargin{4};
        
switch varargin{1}
    case 'drive_i'
        while Y(round(Pfp_rmiv_i(1)),round(Pfp_rmiv_i(2)),round(Pfp_rmiv_i(3)))==0
            Pfp_rmiv_i = Pfp_rmiv_i + prec*Pd_rmiv_i;
        end
        out{1} = Pfp_rmiv_i; % boundary sur la direction de Pi et passant par Pi

    case 'find_bound_i'
        while Y(round(Pfp_rmiv_i(1)),round(Pfp_rmiv_i(2)),round(Pfp_rmiv_i(3)))~=0
            Pfp_rmiv_i = Pfp_rmiv_i - prec*Pd_rmiv_i;
        end
        out{1} = (2*Pfp_rmiv_i + prec*Pd_rmiv_i)/2; % boundary sur la direction de Pi et passant par Pi
        out{2} = Pfp_rmiv_i;
        
    case 'find_bound_o'
        [nx ny nz] = size(Y);
        flag0 = 1;
        while Y(max(1,min(nx,round(Pfp_rmiv_i(1)))),max(1,min(ny,round(Pfp_rmiv_i(2)))),max(1,min(nz,round(Pfp_rmiv_i(3)))))~=5 && flag0 
            if nx == round(Pfp_rmiv_i(1)) || 1 == round(Pfp_rmiv_i(1)) || ...
               ny == round(Pfp_rmiv_i(2)) || 1 == round(Pfp_rmiv_i(2)) || ...
               nz == round(Pfp_rmiv_i(3)) || 1 == round(Pfp_rmiv_i(3)) 
                flag0 = 0;
            end
            Pfp_rmiv_i = Pfp_rmiv_i + prec*Pd_rmiv_i;
        end
        out{1} = (2*Pfp_rmiv_i - prec*Pd_rmiv_i)/2; % boundary sur la direction de Pi et passant par Pi
        out{2} = Pfp_rmiv_i- prec*Pd_rmiv_i; %%% on laisse les points � la surface. C'est � dire dans l'air mais avec la peau dans les voxels adjacents!
end