function M = nirs_Mauchly(cbeta,p_value)
try
    [ns0 nS0 nC0 np] = size(cbeta);
    cbeta0 = reshape(cbeta,[ns0*nS0*nC0,np]);
    msk0 = zeros(size(cbeta0));
    msk0(logical(cbeta0)) = 1;
    msk = sum(msk0,1);
    M = zeros(1,np);
    for ip=1:np %loop over lots of points
        if msk(ip) %exclude points with null data
            tmp = reshape(cbeta(:,:,:,ip),[ns0,nS0*nC0]);
            M0 = Mausphercnst(tmp,p_value,0);
            M(ip) = M0.H; %Returns 1 if nonspherical
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem calculating Mauchly sphericity test');
end