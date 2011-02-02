function [bias beta err] = detrend_wave_rev_v03(Yo,Xo,NT,wavelet_depth)
% ------------------------------------------------------------
% [bias beta err] = detrend_wave_rev_v01(Y,X,NT)
%  Y: (N by M) data matrix
%    N: # of time points
%    M: # of channels
%  X: Design matrix (Predictor for task-related signals)
%    NOTE: Please exclude the 'constant' column in X
%  NT: # of trial. Usually 4 is sufficient.
%    For a fast application, use 3.
%
% Revised version of
%  wavelet-detrending algorithm for long-signals
%   Especially for NIRS data
%
% Based on
%  "Wavelet-MDL detrending for near infrared spectroscopy"
%  to be appeared at Journal of biomedical optics, 2009.
%
% Revised in 2009.01.10. Kwang Eun Jang
% ------------------------------------------------------------

% ----- Basic size parameters ----- %
n_ch = size(Yo,2);    % # of channels
n_B = size(Xo,2);     % # of basis (predictor)
m_Yo = size(Yo,1);    % # of time point (Original)
m_W = wavelet_depth; %9;              % support size of wavelet filter
N_org = size(Yo,1);

% ----- Determine the maximum level of decomposition ----- %
%   Constraint:
%   the last level for which at least one coefficient is correct..
%
% Note: CDF 9/7 filter length = 9
%       level 1: Finest
%       level J: Coarsest
%
% Step: Symmetric extension-1
%       (For boundary problem. Just adjust the size parameter)
%      -> Find the # of wavelet coef. (initial)
%      -> Data elongation using symmetric extension-2
%         (for determine the maximum level of decomposition)

% Find the # of wavelet coef. (initial)
decJ_tmp = fix(log(m_Yo/(m_W-1))/log(2)); % initial

Jtmp = decJ_tmp+1;
tmp_1 = 0;            % # of detail coef. for j-th level
tmp_2 = m_Yo;         % # of residual coef. for j-th level

% Let's find the # of wavelet coef. in Jtmp
for i = Jtmp:-1:2
    tmp_1 = floor(tmp_2/2);
    tmp_2 = tmp_2-tmp_1;
end
ini_m = tmp_1; % This is the # of wavelet coef. in maximum
% level. (initial)

% So the # of data elongation should be
if ini_m == (m_W-1)
    disp('HERE');    % Happen to be..
    dm_left = 200;
    dm_right = 200;
else
    %tmp_1 = 2*m_W-ini_m;  %(M)
    tmp_1 = 2*m_W-ini_m-2;  %(M-1)
    tmp_2 = tmp_1*(2^decJ_tmp);    % this

    dm_left = ceil(tmp_2/2);       % # for decomposition
    dm_right = tmp_2 - dm_left;
end


% data elongation using half-point symmetric extension
Y = [Yo(dm_left:-1:1,:);
    Yo
    Yo(end:-1:end-dm_right+1,:)];

X = [Xo(dm_left:-1:1,:);
    Xo
    Xo(end:-1:end-dm_right+1,:)];


% # of data
N = size(Y,1)

% ----- New decomposition (maximum level) ----- %
decJ = fix(log(N/(m_W-1))/log(2));
%decJ = decJ_tmp+1;
J = decJ+1;          % # of decomposition
L = zeros(J,1);      % Container for # of wavelet coef.

tmp = N;
for i = J:-1:2
    L(i) = floor(tmp/2);
    tmp = tmp-L(i);
end
L(1) = tmp;

% Wavelet Decomposition
Y_wt = cell(1,n_ch);
for i = 1:n_ch
    Y_wt{i} = waveletcdf97(Y(:,i),decJ);
end

% Level-dependent power estimation (pre-calculation)
level_var = zeros(J,n_ch);

sIDX = 1;
for j = 1:n_ch
    level_var(1,j) = 1/0.6745 * median(abs(Y_wt{j}(1:L(1))));
end

for i = 2:J
    sIDX = sum(L(1:(i-1)))+1;
    eIDX = sIDX + L(i)-1;

    for j = 1:n_ch
        level_var(i,j) = 1/0.6745 * median(abs(Y_wt{j}(sIDX:eIDX)));
    end
end

for j = 1:n_ch
    level_var(1,j) = level_var(2,j);
end
%}

% save (weighting vector)
wgt_ce = cell(1,n_ch);
for i = 1:n_ch
    wgt_ce{i} = zeros(N,1);
    wgt_ce{i}(1:L(1)) = level_var(1,i);
    for j = 2:J
        sIDX = sum(L(1:(j-1)))+1;
        eIDX = sIDX + L(j) - 1;
        wgt_ce{i}(sIDX:eIDX) = level_var(j,i);
    end
end

% Finest level of trends  (In the paper, J0)
JF = max(J-NT+1,1);

% # of wavelet coef for CDF9/7: Nw = N
Nw = N;

% Wavelet transform of Basis
Wx = zeros(Nw,n_B);
for i = 1:n_B
    Wx(:,i) = waveletcdf97(X(:,i),decJ);
end

% ----- Prior distribution of coefficient ----- %
% using 'universal prior for integers'
P = modi_unvPrior(Nw,L);    % function of scale

% ----- Initial Detrending... ----- %
% Evaluate MDL value for each level and
%  evaluate within the level, too.
%  (small wavelet coef. will be delete first)
%
% Save the MDL value and the # of deletion within the scale

% delete # recording
dNum = zeros(NT,n_ch);
pXo = (Xo'*Xo)\Xo'; %PP
for ch = 1:n_ch
    disp(strcat('Detrending... ch ',num2str(ch)));

    % wavelet transform of data
    %wy = waveletcdf97(Y(:,ch),decJ);
    wy = Y_wt{ch};


    for J0 = JF:J
        iter = J0 - JF + 1;  % iteration flag
        disp(strcat('- Iteration... ',num2str(iter)));

        % ----- construct 'A' matrix & find a ML solution ----- %

        index_L = J-J0+1;               % scale @ 'L'vector
        n_Full = sum(L(1:index_L));        % # of total coef. in J0-scale

        n_A = n_Full+n_B;
        A = zeros(Nw,n_A);

        for i = 1:n_Full
            A(i,i) = 1;
        end;

        for i = n_Full+1:n_A
            A(:,i) = Wx(:,i-n_Full);
        end

        % Weighted A
        wA = zeros(size(A));
        for i = 1:n_A
            wA(:,i) = (1./wgt_ce{ch}).*A(:,i);
        end


        % find a ML solution
        sol = (A'*wA)\(wA'*wy); %PP % inv(A'*V*A)*A'*V*wy.


        % bias coef.
        c_Full = zeros(Nw,1);          % Full coef. in J0-scale
        c_Full(1:n_Full) = sol(1:end-n_B);

        % beta estimation
        beta_w = sol(end-n_B+1:end);


        % ----- To Measure the fidelity term in MDL ----- %

        %figure; plot(res(1:100)); pause;
        % time?
        %{
       bias_t = waveletcdf97(c_Full,-decJ);
       [us1 us2 res_t] = ...
	   regress(Y(:,ch)-bias_t, X);
       sig = 1/Nw * res_t'*res_t;
       term1 = Nw / 2 * log2(sig);
        %}
        
        bias_t = waveletcdf97(c_Full,-decJ);
        bias_t = bias_t(dm_left+1:end-dm_right);
        try
            [us1 us2 res_t] = ...
                regress(Yo(:,ch)-bias_t, Xo);
        catch
            dY = Yo(:,ch)-bias_t;
            res_t = dY - Xo*pXo*dY;
        end
        sig = 1/N_org * (res_t'*res_t);
        term1 = N_org / 2 * log2(sig);
        %}
        % Penalty term
        term2 = n_Full*(1/2 * log2(Nw)) + ...   % Standard MDL term
            n_Full*sum(-log2(P(1:n_Full))); % From the Prior

        % Saito's MDL
        %{
       %disp('Uniform');
       term2 = n_Full*(3/2 * log2(Nw));
        %}


        % ----- Let's measure MDL within the scale.----- %
        %  To find out it is possible to reduce # of coef.
        %  WithIN given J0-scale.
        MDLwn = zeros(L(index_L),1);
        MDLwn(1) = term1 + term2;

        % sorting along the magnitude of wavelet coef. of trend
        % small coef. will be removed if MDL value gets small
        n_Prev = sum(L(1:index_L-1));
        [useless1 ind] = sort(abs(c_Full(n_Prev+1:n_Full)),'ascend');

        for dcoef = 1:(L(index_L)-1)  % delete of coef.
            n0 = n_Full - dcoef;

            dIND = ind(1:dcoef);         % will be del
            dIND = sort(dIND,'ascend');

            colIND = (1:(n_Full+n_B))';  % position of remain coef
            for del = 1:dcoef
                delPos = dIND(del) - del+1;
                colIND = [colIND(1:n_Prev+delPos-1);colIND(n_Prev+delPos+1:end)];
            end
            Ac = A(:,colIND);            % remain A matrix
            wAc = wA(:,colIND);

            % find a ML solution, again.
            sol = (Ac'*wAc)\(wAc'*wy); %PP
            %sol = inv(Ac'*Ac)*Ac'*wy;

            % bias coef.
            biasCoef = zeros(Nw,1);
            biasCoef(colIND(1:end-n_B)) = sol(1:n_Full-dcoef);

            % beta estimation
            beta_w = sol(end-n_B+1:end);

            % time
            bias_t = waveletcdf97(biasCoef,-decJ);

            %{
	 [us1 us2 res_t] = ...
	     regress(Y(:,ch)-bias_t, X);
	 sig = 1/Nw * res_t'*res_t;
	 term1 = Nw / 2 * log2(sig);
            %}
            
            bias_t = bias_t(dm_left+1:end-dm_right);
            try
                [us1 us2 res_t] = ...
                    regress(Yo(:,ch)-bias_t, Xo);
            catch
                dY = Yo(:,ch)-bias_t;
                res_t = dY - Xo*pXo*dY;
            end
            sig = 1/N_org * (res_t'*res_t);
            term1 = N_org / 2 * log2(sig);
            %}

            %disp('Universal');
            term2 = n0*(1/2 * log2(Nw))+...        % standard
                n0*sum(-log2(P(1:n0)));        % Prior

            % Saito's MDL
            %{
	 %disp('Uniform');
	 term2 = n0*(3/2 * log2(Nw))
            %}
            MDLwn(1+dcoef) = term1 + term2;
        end

        %}
        % Where is minimum MDL within J0-scale?
        [useless1 minwn] = min(MDLwn);
        MDLS(iter,ch) = MDLwn(minwn);  % save MDL value
        dNum(iter,ch) = minwn-1;       % save # of deletion
    end
end

% ----- Actual Detrending... ----- %
% Based on the earlier step, find the trend signal
bias = zeros(m_Yo,n_ch);
beta = zeros(n_B,n_ch);
err = zeros(1,n_ch);

for ch = 1:n_ch
    [tmpVal tmpPos] = min(MDLS(:,ch));
    minPos = mod(tmpPos-1,NT)+1;            % min scale

    wy = waveletcdf97(Y(:,ch),decJ);        % wavelet{Y}
    index_L = J-(minPos+JF-1)+1;               % iteration number to

    % If the number of trial is too small,
    if index_L == NT
        disp(strcat('----- warning!! ---- @ ',num2str(ch)));
    end

    n_Full = sum(L(1:index_L));        % # of coefficients that describe the trend

    % construct 'A' matrix
    n_A = n_Full+n_B;
    A = zeros(Nw,n_A);
    for i = 1:n_Full
        A(i,i) = 1;
    end;

    for i = n_Full+1:n_A;
        A(:,i) = Wx(:,i-n_Full);
    end

    % Weighted A
    wA = zeros(size(A));
    for i = 1:n_A
        wA(:,i) = (1./wgt_ce{ch}).*A(:,i);
    end

    % find a ML solution
    sol = (A'*wA)\(wA'*wy); %PP
    %sol = inv(A'*A)*A'*wy;


    c_Full = zeros(Nw,1);
    c_Full(1:n_Full) = sol(1:end-n_B);

    if dNum(minPos,ch) == 0
        % no need to delete
        biasTmp = waveletcdf97(c_Full,-decJ);
    else
        % sorting
        n_Prev = sum(L(1:index_L-1));
        [useless ind] = sort(abs(c_Full(n_Prev+1:n_Full)),'ascend');

        % # of deletion
        dcoef = dNum(minPos,ch);
        dIND = ind(1:dcoef);        % will be del
        dIND = sort(dIND,'ascend');
        colIND = [1:(n_Full+n_B)]';  % for solve
        for del = 1:dcoef
            delPos = dIND(del) - del+1;
            colIND = [colIND(1:n_Prev+delPos-1);colIND(n_Prev+delPos+1:end)];
        end
        Ac = A(:,colIND);
        wAc = wA(:,colIND);

        % find a solution
        sol = (Ac'*wAc)\(wAc'*wy); %PP
        %sol = inv(Ac'*Ac)*Ac'*wy;


        % bias coef.
        biasCoef = zeros(Nw,1);
        biasCoef(colIND(1:end-n_B)) = sol(1:n_Full-dcoef);
        biasTmp = waveletcdf97(biasCoef,-decJ);
    end
    bias(:,ch) = biasTmp(dm_left+1:end-dm_right);
    try
        [beta(:,ch) useless errTmp] = regress(Yo(:,ch)-bias(:,ch),Xo);
    catch
        dY = Yo(:,ch)-bias(:,ch);
        beta(:,ch) = pXo*dY;
        errTmp = dY - Xo*beta(:,ch);
    end
    err(ch) = 1/m_Yo * sum(errTmp(:).^2);
end




