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
        
        while Y(round(Pfp_rmiv_i(1)),round(Pfp_rmiv_i(2)),round(Pfp_rmiv_i(3)))~=5
            Pfp_rmiv_i = Pfp_rmiv_i + prec*Pd_rmiv_i;
        end
        out{1} = (2*Pfp_rmiv_i - prec*Pd_rmiv_i)/2; % boundary sur la direction de Pi et passant par Pi
        out{2} = Pfp_rmiv_i;
        
    otherwise
        disp('Whaaat !!!!');
        out = [1 1 1];
end