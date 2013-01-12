function varargout = approx_trRV(varargin)
%Difficulty with calculation of TrRV and especially TrRVRV is that the
%nScan values that contribute to the traces are very badly distributed,
%with a few very large values. Strategy is to calculate the values that
%contribute to TrRV exactly, and since the indices where they occur will
%be the same as for trRVRV, approximate trRVRV by those values.
target_sampling_rate = 0;
if nargin == 6
    SPM = varargin{6};
    if isfield(SPM.job,'target_sampling_rate')
        if isfield(SPM.job.target_sampling_rate,'specified_sampling_rate')
            target_sampling_rate = 1;
            sampling_rate = SPM.job.target_sampling_rate.specified_sampling_rate.sampling_rate;
        end
    end
end
% if nargin == 5 || nargin == 6
%     exact = varargin{5};
% else
%     exact = 0;
% end
%Now always enforce exact calculation
exact = 1;

X = varargin{1};
pKX = varargin{2};
S = varargin{3};
if nargin == 4
    %F-test 
    c = varargin{4};
end

if nargout == 2
    trRVRVon = 1;
else
    trRVRVon = 0;
end
var1 = S * S';
if exist('c','var')
    sample_size = 100;
    X0 = X*(eye(size(c,1))-c*c');
    xKXs0 = spm_sp('Set',X0);
    pKX0 = spm_sp('x-', xKXs0);   
    var2 = X*(pKX *var1) - X0 * (pKX0 * var1); % RV * e(kk)   
else
    sample_size = 100;
    var2 = var1 - X * (pKX * var1); % RV * e(kk) 
end
RVd = diag(var2);
trRV = sum(RVd); %exact
varargout{1} = trRV;

if target_sampling_rate
    %downsample the variables X, pKX and S
    fs = SPM.fs;
    decimate_factor = round(fs/sampling_rate);
    X = X(round(decimate_factor/2):decimate_factor:end,:);
    pKX = pKX(:,round(decimate_factor/2):decimate_factor:end);
    S = S(round(decimate_factor/2):decimate_factor:end,:);
    S = S(:,round(decimate_factor/2):decimate_factor:end);
    %recalculate trRV
    var1 = S * S';
    if exist('c','var')
        sample_size = 100;
        X0 = X*(eye(size(c,1))-c*c');
        xKXs0 = spm_sp('Set',X0);
        pKX0 = spm_sp('x-', xKXs0);
        var2 = X*(pKX *var1) - X0 * (pKX0 * var1); % RV * e(kk)
    else
        sample_size = 100;
        var2 = var1 - X * (pKX * var1); % RV * e(kk)
    end
    RVd = diag(var2);
    trRV = sum(RVd); %exact
end
%computation of var1 and var2 is fairly quick - actually, no: var1 takes
%about 3.5 minutes. Why is it slower this time than in the earlier test?
nScan = size(X,1);

if trRVRVon
    if exact
        var3 = var1 * var2;
        if exist('c','var')
            var4 = X * (pKX * var3) - X0 * (pKX0 * var3);
        else
            var4 = var3 - X * (pKX * var3);
        end
        RVRVe = diag(var4);
        trRVRV = sum(RVRVe);
    else
        %approximation
        %var3 = var1*var2 is the long calculation - instead...
        %find extreme values (k2) of var2    
        RVs = std(RVd);
        RVm = median(RVd);
        k2 = abs(RVd-RVm)>RVs; %exceeding one standard deviation 
        %largest contributions to var3
        var3e = var1 * var2(:,k2); %e for extreme
        if exist('c','var')
            var4e = X * (pKX * var3e) - X0 * (pKX0 * var3e);
        else
            var4e = var3e - X * (pKX * var3e);
        end
        RVRVe = diag(var4e(k2,:));
        %find typical value of RVRV contribution:
        var3t = var1 * var2(:,1:sample_size); %t for typical
        if exist('c','var')
            var4t = X * (pKX * var3t) - X0 * (pKX0 * var3t);
        else
            var4t = var3t - X * (pKX * var3t);
        end
        RVRVt = median(diag(var4t(1:sample_size,:))); %median is insensitive to extreme values
        %trRVRV approximated by sum of extreme values and contribution from typical
        %values
        trRVRV = sum(RVRVe) + RVRVt*(nScan-length(find(k2)));

        %In one example, "exact" calculation of trRV, trRVRV gave:
        % 341.25, 238.84
        %while truncation to 100 random samples gave
        % 334.19 , 202.19
        %and trunction to 1000 random samples gave
        % 339.53, 206.81
        %so the error is as much as 15%
    end
end

if trRVRVon
    if target_sampling_rate
        varargout{2} =  1.0280 * varargout{1}^2 * trRVRV/(trRV)^2; %numerical factor is an estimated correction 
    else
        varargout{2} = trRVRV;
    end
end