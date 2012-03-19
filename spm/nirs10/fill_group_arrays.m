function [cbeta ccov_beta new_version] = fill_group_arrays(TOPO,big_TOPO,v1,c1,h1,Z,xCon,ns,shb)
Sess = TOPO.Sess;
Cp = TOPO.Cp;
new_version = 0;
if shb
    sign_hb = 1;
else
    sign_hb = -1;
end
%Why more than 8 cases???
%Either T or F stat
%Either new (LKC) or old (pre-LKC) version
%Group of subjects or group of sessions
%and for group of subjects, one needs to distinguish if the data come from
%individual sessions at the subject level, or from group of sessions
if xCon(c1).STAT == 'T'
    fc = 0; %used only for FFX || nS==1
    %fill in cbeta and ccov_beta
    for f1=1:ns
        if Z.FFX || Z.nS==1
            %select sessions which had the contrast
            if any(f1==Sess{c1})
                fc = fc+1;
                if isfield(TOPO.v{v1}.s{f1}.hb{h1},'stat_map')
                    tmp = sign_hb*squeeze(TOPO.v{v1}.s{f1}.hb{h1}.beta_map(Cp{c1,f1},:,:));
                    if f1 == 1 %set array
                        cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(fc,:) = tmp(:);
                    if ~Z.simple_sum
                        tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.beta_map(Cp{c1,f1},:,:))./squeeze(TOPO.v{v1}.s{f1}.hb{h1}.stat_map(Cp{c1,f1},:,:));
                    end
                    ccov_beta(fc,:) = tmp(:).^2;
                    new_version = 1;
                else
                    %now use Cp{c1,f1} to access the required c_interp_beta instead of c1
                    tmp = sign_hb*squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_beta(Cp{c1,f1},:,:));
                    if f1 == 1
                        cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(fc,:) = tmp(:);
                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_cov_interp_beta(Cp{c1,f1},:,:));
                    ccov_beta(fc,:) = tmp(:);
                end
            end
        else
            if ~isfield(big_TOPO{f1}.v{v1},'s')
                if isfield(big_TOPO{f1}.v{v1}.g{1}.hb{h1},'stat_map')
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.beta_map(c1,:,:));
                    if f1 == 1
                        cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(f1,:) = tmp(:);
                    if ~Z.simple_sum
                        tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.beta_map(c1,:,:))./squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.stat_map(c1,:,:));
                    end
                    ccov_beta(f1,:) = tmp(:).^2;
                    new_version = 1;
                else
                    %group analysis of a group of sessions analysis
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_interp_beta(c1,:,:));
                    if f1 == 1
                        cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(f1,:) = tmp(:);
                    tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                    ccov_beta(f1,:) = tmp(:);
                end
            else
                %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                is1 = Z.group_session_to_average;
                %do each session separately
                if isfield(big_TOPO{f1}.v{v1}.s{is1}.hb{h1},'stat_map')
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.beta_map(c1,:,:));
                    if f1 == 1
                        cbeta = zeros(1,length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(f1,:) = tmp(:);
                    if ~Z.simple_sum
                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.beta_map(c1,:,:))./squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.stat_map(c1,:,:));
                    end
                    ccov_beta(f1,:) = tmp(:).^2;
                    new_version = 1;
                else
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_beta(c1,:,:));
                    if f1 == 1
                        cbeta = zeros(1,length(tmp(:)));
                        ccov_beta = cbeta;
                    end
                    cbeta(f1,:) = tmp(:);
                    tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                    ccov_beta(f1,:) = tmp(:);
                end
            end
        end
    end
else %quick fix for F stats
    fc = 0; %used only for FFX || nS==1
    ccov_beta = [];  %not used
    for f1=1:ns
        if Z.FFX || Z.nS==1
            %select sessions which had the contrast
            if any(f1==Sess{c1})
                fc = fc+1;
                tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_F(Cp{c1,f1},:,:));
                if f1 == 1
                    cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                end
                cbeta(fc,:) = tmp(:);
            end
        else
            if ~isfield(big_TOPO{f1}.v{v1},'s')
                %group analysis of a group of sessions analysis
                tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_interp_F(c1,:,:));
                if f1 == 1
                    cbeta = zeros(length(Sess{c1}),length(tmp(:)));
                end
                cbeta(f1,:) = tmp(:);
            else
                %do each session separately
                is1 = Z.group_session_to_average;
                %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_F(c1,:,:));
                if f1 == 1
                    cbeta = zeros(1,length(tmp(:)));
                end
                cbeta(f1,:) = tmp(:);
                %end
            end
        end
    end
end