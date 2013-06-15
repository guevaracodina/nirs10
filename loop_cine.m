function C = loop_cine(Q,W,Z)
try
    s1 = W.s1;
    s2 = W.s2;
    nD = size(Q.ibeta,1); %number of delays
       
    Q.ivar(isnan(Q.ivar)) = max(Q.ivar(~isnan(Q.ivar)));
    Q.ibeta(isnan(Q.ibeta)) = 0;
    Q.ivar(Q.ivar < 0) = max(Q.ivar(:)); %-Q.ivar(Q.ivar < 0);
    stat = Q.ibeta./sqrt(Q.ivar);    
    if ~isreal(stat)
        stat = real(stat);
        disp('Some imaginary values removed in t-stat calculation');
    end
    
    stat_map = zeros(nD,s1,s2); %T-stat
    beta_map = zeros(nD,s1,s2);
    cov_map = zeros(nD,s1,s2);
    
    stat_tmp = zeros(s1,s2);
    beta_tmp = zeros(s1,s2);
    cov_tmp = zeros(s1,s2);
 
    for c1=1:nD
        stat_tmp(Q.index_mask) = stat(c1,:);
        beta_tmp(Q.index_mask) = Q.ibeta(c1,:);
        cov_tmp(Q.index_mask) = Q.ivar(c1,:);        
        stat_map(c1,:,:) = stat_tmp;
        beta_map(c1,:,:) = beta_tmp;
        cov_map(c1,:,:) = cov_tmp;
    end
    
    if Z.spatial_LPF %does not work properly -- do not use
        K.k1 = s1;
        K.k2 = s2;
        K.radius = Z.radius;
        K = spatial_LPF('set',K);
        stat_map = spatial_LPF('lpf',K,stat_map);
        beta_map = spatial_LPF('lpf',K,beta_map);
        cov_map = spatial_LPF('lpf',K,cov_map);
    end
    C.stat_map = stat_map;
    C.beta_map = beta_map; %Used at the group level
    C.cov_map = cov_map;
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end