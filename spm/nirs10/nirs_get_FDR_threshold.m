function [p_value flag first_k_value stop_k_value] = nirs_get_FDR_threshold(FDR_type,q_value,Nc,p)
flag = 0;
switch FDR_type
    case 1
        for i0 = Nc:-1:1
            if p(i0) <= i0*q_value/Nc
                flag = 1;
                break
            end
        end
    case 2
        qp = 0;
        for i0 = 1:Nc;
            qp = qp+ 1/i0;
        end
        for i0 = Nc:-1:1
            if p(i0) <= i0*q_value/(Nc*qp)
                flag = 1;
                break
            end
        end
end
first_k_value = Nc;
stop_k_value = i0;
p_value = p(i0); %t_value = spm_invTcdf(p_value,erdf);