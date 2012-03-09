%nirs_spm_Volterra
SPM = [];
%SPM's session
s = 1;
fs = 20;
volt = 2;
ns = 1000;
norm_bf = 4.7506; % = 1/max(X(:,1))
SPM.nscan(s) = ns;
xBF.T = 1;
xBF.dt = 1/(fs*xBF.T); % - time bin length {seconds}
xBF.name = 'hrf'; %description of basis functions specified 
%xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);
bf  = xBF.bf;

SPM.xBF = xBF;
SPM.xBF.T = xBF.T; %get more precision on onsets position      
SPM.xBF.T0 = 1; %shouldn't need it - no offset
SPM.xBF.UNITS = 'secs';
SPM.xBF.Volterra = volt;
SPM.xY.RT = 1/fs; %xBF.dt;

U.name = {'one'};
U.dt = xBF.dt;
U.ons = [0 1];
U.dur = 0;
U.P.name = 'none'; 
U.P.h    = 0;
SPM.Sess(s).U = U;

%copied code from spm_get_ons(SPM,s):
% create stimulus functions (32 bin offset)
%======================================================================
U = spm_get_ons(SPM,s);
V = SPM.xBF.Volterra;
%convolve stimuli U with basis functions
[X,Xn,Fc] = nirs_spm_Volterra(U,bf,V); 
try
    X = X((0:(ns - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
end
SPM.xX.X = X;

try
    tSPM = SPM;
    tSPM.Sess(s).U.ons = U.ons(1); %first stim
    U1 = spm_get_ons(tSPM,s);
    [Xt,Xn,Fc] = nirs_spm_Volterra(U1,bf,V); 
    Xt = Xt((0:(ns - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
end
try
    uSPM = SPM;
    uSPM.Sess(s).U.ons = U.ons(2); %second stim
    U2 = spm_get_ons(uSPM,s);
    [Xu,Xn,Fc] = nirs_spm_Volterra(U2,bf,V); 
    Xu = Xu((0:(ns - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
    Xtu = Xt+Xu;
end

% and orthogonalise (within trial type)
%--------------------------------------
for i = 1:length(Fc)
    X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
end 
bf_norm = 1/max(X(:,1));
SPM.xsDes = struct('Basis_functions', SPM.xBF.name, 'Sampling_period_sec', ...
    num2str(SPM.xY.RT), 'Total_number_of_samples', num2str(SPM.nscan));
%spm_DesRep('DesMtx',SPM.xX,[],SPM.xsDes)
%figure; plot(linspace(0,50,1000),X)
X1 = X(:,1);
X2 = X(:,2);
corr(X1,X2)
Y = spm_orth(X);
Y1 = Y(:,1);
Y2 = Y(:,2);
corr(Y1,Y2)
%figure; plot(linspace(0,50,1000),Y)
%rescale
rX = X./(ones(size(X,1),1)*max(X,[],1));
%figure; plot(linspace(0,50,1000),rX)
rY = Y./(ones(size(Y,1),1)*max(Y,[],1));
%figure; plot(linspace(0,50,1000),rY)
rhoX = corr(X)
rhoY = corr(Y)
Z1 = X(:,1)*norm_bf;
Z2 = X(:,2)*norm_bf^2;

try
    %figure; plot(linspace(0,50,1000),Xtu,'r',linspace(0,50,1000),X,'k')
    figure; plot(linspace(0,50,1000),Xtu(:,1)-Xtu(:,2),'r',...
        linspace(0,50,1000),X(:,1)-X(:,2),'k',linspace(0,50,1000),X(:,2)-Xtu(:,2),'g')
end
try 
    %figure; plot(linspace(0,50,1000),Y)
end