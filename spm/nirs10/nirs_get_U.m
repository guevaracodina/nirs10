function U = nirs_get_U(TR,k,U)
ons = U.ons;
%Added division by TR; is that correct?
ons = ons(:)/TR; %what happens when there are more than 1 onset type?
dur = U.dur/TR;
dur = dur(:);
xP.name = 'none';
xP.h = 0;
T = 1;
dt = TR;
v = length(U);
for i = 1:v
    % peri-stimulus times {seconds}
    %----------------------------------------------------------------------
    pst   = [0:(k-1)]*T*dt - min(ons)*TR;
    for j = 1:length(ons)
        w      = [0:(k-1)]*T*dt - ons(j)*TR;
        v      = find(w >= 0);
        pst(v) = w(v);
    end
    
    
    % interaction with causes (u) - 1st = main effects
    %----------------------------------------------------------------------
   u     = ons.^0;
    
    % orthogonalize inputs
    %----------------------------------------------------------------------
    u      = spm_orth(u);
    
    % and scale so sum(u*dt) = number of events, if event-related
    %----------------------------------------------------------------------
    if ~any(dur)
        u  = u/dt;
    end
    
    % create stimulus functions (32 bin offset)
    %======================================================================
    ton       = round(ons*TR/dt) + 33;               % onsets
    tof       = round(dur*TR/dt) + ton + 1;          % offset
    sf        = sparse((k*T + 128),size(u,2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
        end
        if size(sf,1) > tof(j)
            sf(tof(j),:) = sf(tof(j),:) - u(j,:);
        end
    end
    sf        = cumsum(sf);                         % integrate
    sf        = sf(1:(k*T + 32),:);                 % stimulus
    
    % place in ouputs structure
    %----------------------------------------------------------------------
    U(i).dt   = dt;         % - time bin {seconds}
    U(i).u    = sf;         % - stimulus function matrix
    U(i).pst  = pst;        % - pst (seconds)
    U(i).P    = xP;         % - parameter struct
end

U.u = U.u(33:end); %?
