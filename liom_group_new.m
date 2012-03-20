function G = liom_group_new(cbeta,s1,s2)
try
    [ns np] = size(cbeta);
    tmp_mask = zeros(size(cbeta));
    tmp_mask(cbeta == 0) = 1;
    %mask counts the number of missing subjects at each pixel
    mask = sum(tmp_mask,1);
    
    calc_by_pixel = 0;
    if ~calc_by_pixel
        index_group = find(mask == 0); %intersection of all the subjects
        %design matrix
        X = ones(1,ns);
        pX = pinv(X);
        b = pX'*cbeta; %simple average
        s = std(cbeta,0,1);
        
        
        res = cbeta - X' * b;
        res = reshape(res,ns,s1,s2);
        % LKC calculation
        L2 = calc_LKC(index_group,[s1 s2], res, 'group');
        r = sqrt(L2./pi);
        L1 = pi * r;
        L0 = 1;
    else
        %try to improve on LKC
        X = ones(1,ns);
        min_subj = 4;
        b = zeros(1,np);
        s = zeros(1,np);
        sbeta = cbeta;
        sbeta(logical(tmp_mask)) = NaN;
        res = cbeta;
        eidfi = zeros(1,ns-min_subj+1);
        erdfi = zeros(1,ns-min_subj+1);
        first_idx = [];
        for i0=1:(ns-min_subj+1)
            %indices of pixels where exactly i0-1 subjects are missing
            %at least min_subj subjects must be present
            idx = find(mask == (i0-1));
            b(idx) = sum(cbeta(:,idx),1)/(ns-i0+1);
            s(idx) = std(sbeta(:,idx),1);
            res(:,idx) = cbeta(:,idx) - X' * b(1,idx);
            L2i = calc_LKC(idx,[s1 s2], reshape(res,ns,s1,s2), 'group');
            r = sqrt(L2i./pi);
            L1i = pi * r;
            L0i = 1;
            LKCi{i0} = [L0i L1i L2i];
            idxi{i0} = idx;
            eidfi(i0) = 1;
            erdfi(i0) = ns-i0;
            if isempty(first_idx)
                if ~isempty(idx)
                    first_idx = idx;
                end
                L2 = L2i;
                L1 = L1i;
                L0 = L0i;
            end
        end
        %L2b = []; for i0=1:(ns-min_subj+1) , L2b = [L2b LKCi{i0}(3)]; end; figure; plot(L2b)
        G.LKCi = LKCi;
        G.idxi = idxi;
        G.eidfi = eidfi;
        G.erdfi = erdfi;
        G.min_subj = min_subj;
    end
    beta_group = reshape(b,s1,s2);
    std_group = reshape(s,s1,s2);
    t = b./s;
    t(isnan(t)) = 0;
    
    G.ns = ns;
    G.erdf_group = ns-1;
    G.var_bs = [];
    G.beta_group = beta_group;
    G.std_group = std_group;
    G.tmap_group = reshape(t,s1,s2);
    G.LKC = [L0 L1 L2];
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end