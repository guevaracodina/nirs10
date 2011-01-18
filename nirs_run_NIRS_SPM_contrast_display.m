function out = nirs_run_NIRS_SPM_contrast_display(job)
% 'Contrast estimates and 90% C.I.'
view = job.view;
preproc_file =job.NIRS_SPM_Coregistration_Channels{1,1};
data_file = job.data_file{1,1};
map_file = job.map_file{1,1};
cinterp_SPM_nirs = [];
preproc_info = [];
SPM_nirs = [];
try
    load(preproc_file);
    preproc_info;
    %give chpos.rchn and chpos.cchn (2D projection coordinates as row
    %and column)
    chpos = preproc_info.rend_ch_pos{view};
catch
    disp('Could not load coregistration file');
end
try
    load(map_file);
    try
        %in case this map comes from a group estimation    
        tstat = SPM_nirs.nirs.tstat;
    catch
        try 
            cinterp_SPM_nirs;
            %need to calculate tstat: (without correction for tube formula
            %-- see nirs_run_NIRS_SPM_contrast for code (which could be
            %saved in cinterp_ files if desired
            tstat = cinterp_SPM_nirs.cbeta ./sqrt(cinterp_SPM_nirs.ccovbeta);
            tstat = reshape(tstat,cinterp_SPM_nirs.s1,cinterp_SPM_nirs.s2);
        catch
            disp('Content of map file is unrecognized.');
        end
    end
catch
    disp('Could not load statistical map file');
end
try
    load(data_file);
    SPM_nirs;
catch
    disp('Could not load data file');
end
%regressors to display
reg_num = job.reg_num;


%
%view = 4;
%load();

%Specify channel to look at
%Find the channel with the minimum distance to the map's maximum

[r1 c1] = find(tstat == max(max(tstat)));
d1 = zeros(size(chpos.rchn,1),1);
for i1=1:size(chpos.rchn,1)
    d1(i1) = sum(([r1; c1]-[chpos.rchn(i1); chpos.cchn(i1)]).^2).^0.5;
end
[vd1 id1] = min(d1);
disp(['Channel is: ' int2str(id1) ' and distance to max is: ' num2str(vd1)]); 
%id1 is now the index of the channel closest to the maximum of the stat map

beta = SPM_nirs.nirs.beta(:,id1);
%PP use channel-wise normalization: ResMS = ResSS/trRV
Bcov = (SPM_nirs.nirs.ResSS(id1)/SPM_nirs.xX.trRV)*SPM_nirs.xX.Bcov;
con =  zeros(size(beta,1),size(reg_num,2)); %SPM.xCon(Ic).c;
for j1=1:size(reg_num,2)
    con(reg_num(j1),j1)=1;
end
CI    = 1.6449;                 % = spm_invNcdf(1 - 0.05);
%-Colour specifications and index;
%--------------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 .5 .5];

TITLE = 'Contrast estimates';

% compute contrast of parameter estimates and 90% C.I.
%------------------------------------------------------------------
cbeta = con'*beta;
SE    = sqrt(diag(con'*Bcov*con));
CI    = CI*SE;

contrast.contrast      = cbeta;
contrast.standarderror = SE;
contrast.interval      = 2*CI;
assignin('base','contrast',contrast)

% bar chart
%------------------------------------------------------------------
figure; %(Fgraph)
%subplot(2,1,2)
cla
hold on

% estimates
%------------------------------------------------------------------
h     = bar(cbeta);
set(h,'FaceColor',Col(2,:))

% standard error
%------------------------------------------------------------------
for j = 1:length(cbeta)
    line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
        'LineWidth',6,'Color',Col(3,:))
end

title(TITLE,'FontSize',12)
xlabel('contrast')
ylabel(['contrast estimate at channel ',int2str(id1)])
set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
hold off

% set Y to empty so outputs are assigned
%------------------------------------------------------------------
%Y = [];

% all fitted effects or selected effects

%%%%%%%%%%%%%%%%%%%
% Volterra kernels display
%%%%%%%%%%%%%%%%%
output_Volterra = 0;
if output_Volterra
    dt = SPM_nirs.xBF.dt; 
    u =  3; %LFCP x LFCP  %PP length(Sess(s).Fc);
    % Parameter estimates and basis functions
        %------------------------------------------------------------------
        bf    = SPM_nirs.xBF.bf/dt; %PP
        pst   = ([1:size(bf,1)] - 1)*dt;

        % second order kernel
        %------------------------------------------------------------------
        if u > length(SPM_nirs.Sess.U) %PP

            % Parameter estimates and kernel
            %--------------------------------------------------------------
            B     = beta(SPM_nirs.Sess.Fc(u).i); %PP ?
            i     = 1;
            Y     = 0;
            for p = 1:size(bf,2)
                for q = 1:size(bf,2)
                    Y = Y + B(i)*bf(:,p)*bf(:,q)';
                    i = i + 1;
                end
            end

            % plot
            %--------------------------------------------------------------
            figure; %PP (Fgraph)
            %subplot(2,2,3)
            imagesc(pst,pst,Y)
            axis xy
            axis image

            title('2nd order Kernel','FontSize',12);
            xlabel('peristimulus time {secs}')
            ylabel('peristimulus time {secs}')

            %subplot(2,2,4)
            figure;
            plot(pst,Y)
            axis square
            grid on

            title(SPM_nirs.Sess.Fc(u).name,'FontSize',12);
            xlabel('peristimulus time {secs}')


            % first  order kernel
            %--------------------------------------------------------------
        else
            B     = beta(Sess(s).Fc(u).i(1:size(bf,2)));
            Y     = bf*B;

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)
            subplot(2,1,2)
            plot(pst,Y)
            grid on
            axis square

            title({'1st order Volterra Kernel' Sess(s).Fc(u).name},...
                'FontSize',12);
            xlabel('peristimulus time {secs}')
            ylabel(['impulse response',XYZstr])
        end
end
%------------------------------------------------------------------
out = [];
end