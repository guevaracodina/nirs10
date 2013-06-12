function [f1 f2 f3 f4]=nirs_SCKS_display(SCKSall)
try
    SCKS = SCKSall{end,end};
    %Need to fill pU: v and r
    if SCKS.SCKS_OK
        M=SCKS;
        M(1).inSec=1;
        M(1).dt = SCKS.dt;
        nD=SCKS.M(1,1).E.nD;
        v=SCKS.pU.v{1};
        
        x = zeros(size(SCKS.qU.x{1}));
        N = size(x,2);
        %N = length(SCKS.pU.r); %just for fig1, fig3
        time=1:N;
        if M(1).inSec==1
            time=time*M(1).dt;
            if isfield(M(1),'offsetTime')&& ~isempty(M(1).offsetTime), time=time+M(1).offsetTime; end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 1: States
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f1 = spm_figure('Create','Graphics','SCKS State estimates');
        set(f1,'RendererMode','auto','Renderer','painter','toolbar','figure');
        clf(f1);
        for p = 1:2
            subplot(2,1,p),
            hax = gca;
            si    = spm_invNcdf(1 - 0.05);
            
            if p == 1,
                xxfig = SCKS.qU.x{2}(:,1:nD:end);
                Sfig  = SCKS.qU.S{2}(:,1:nD:end);
                tit   = 'Forward Estimate';
            else
                xxfig = SCKS.qU.x{1}(:,1:nD:end);
                Sfig  = SCKS.qU.S{1}(:,1:nD:end);
                tit   = 'Backward Estimate';
            end
            
            s = abs(Sfig);
            
            % conditional covariances
            %------------------------------------------------------------------
            j       = [1:size(xxfig,1)];
            ss      = si*s(j,:);
            s(j,:)  = [];
            [ill indss] = sort(mean(ss,2),'descend');
            pf = plot(time,xxfig,'linewidth',1.5);
            %     set(hax,'xlim',time([1 end]),'nextplot','add')
            hold on
            for ic = 1:size(xxfig,1)
                col0 = get(pf(indss(ic)),'color');
                col = (ones(1,3)-col0)*0.65 + col0;
                
                plot([time ],[xxfig(indss(ic),:)+ss(indss(ic),:) ;xxfig(indss(ic),:)-ss(indss(ic),:)],...
                    'r',...
                    'Color',col);
                %         fill([time fliplr(time)],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
                %             'r',...
                %             'FaceColor',col,...
                %             'EdgeColor',col);
                
                hold on
                COL{ic} = col0;
            end
            for ic = 1:size(xxfig,1)
                plot(time,xxfig(indss(ic),:),'color',COL{ic},'linewidth',0.75);
            end
            h1 = plot(time,(x(:,1:nD:end))','Color',[0 0 0],'linewidth',1);
            
            title(tit);
            grid(hax,'on')
            axis(hax,'tight')
            set(hax,'box','on','Layer','top');
            set(hax,'tickdir','out')
            legend([h1(1) pf(1)],{'True','x1'});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 2: Inputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f2 = spm_figure('Create','Graphics','SCKS Input estimate');
        set(f2,'RendererMode','auto','toolbar','figure');
        clf(f2);
        for p = 1:2
            subplot(2,1,p),
            hax = gca;
            si    = spm_invNcdf(1 - 0.05);
            if p == 1,
                xxfig = SCKS.qU.v{2};
                %         Sfig  = SCKS.qU.S{2};
                Sfig  = SCKS.qU.C{2}; %correction dub
                tit   = 'Forward Estimate';
            else
                xxfig = SCKS.qU.v{1};
                %         Sfig  = SCKS.qU.S{1};%correction dub
                Sfig  = SCKS.qU.C{1};%correction dub
                tit   = 'Backward Estimate';
            end
            xxfig=xxfig(:,1:nD:end);
            time2=1:(length(xxfig));
            if M(1).inSec==1
                time2=time2*M(1).dt;
                if isfield(M(1),'offsetTime')&& ~isempty(M(1).offsetTime), time2=time2+M(1).offsetTime; end
            end
            
            s = abs(Sfig);
            s=s(:,1:nD:end);
            % conditional covariances
            %------------------------------------------------------------------
            j       = [1:size(xxfig,1)];
            ss      = si*s(j,:);
            s(j,:)  = [];
            
            pf = plot(time2,xxfig,'linewidth',1.5);
            %     set(hax,'xlim',time2([1 end]),'nextplot','add')
            hold on
            for ic = 1:size(xxfig,1)
                col0 = get(pf(ic),'color');
                col = (ones(1,3)-col0)*0.65 + col0;
                plot(time2,[(xxfig(ic,:)+ss(ic,:)) ;((xxfig(ic,:)-ss(ic,:)))],...
                    'r',...
                    'Color',col);%,...
                %         fill([time2 fliplr(time2)],[(xxfig(ic,:) + ss(ic,:)) fliplr((xxfig(ic,:) - ss(ic,:)))],...
                %             'r',...
                %             'FaceColor',col,...
                %             'EdgeColor',col);%,...
                %'FaceAlpha',0.25);
                hold on
                plot(time2,xxfig(ic,:),'color',col0,'linewidth',0.75);
            end
            if ~isempty(v)
                v1 = zeros(length(time2),1);
                v = full(v);
                v = [v; zeros(length(time2)-length(v),1)];
                h1 = plot(time2,v(1:length(time2))','Color','r','linewidth',1); %PP red more visible
            end
            hold off
            title(tit); %drawnow;
            grid(hax,'on')
            axis(hax,'tight')
            set(hax,'box','on','Layer','top','tickdir','out');
            legend([h1(1) pf(1)],{'True','Input'});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 3: Response
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f3 = spm_figure('Create','Graphics','SCKS response estimate');
        set(f3,'RendererMode','auto','toolbar','figure');
        clf(f3);
        pf = plot(time,SCKS.qU.r{1}(:,1:nD:end),'o','markersize',5); %estimate
        set(gca,'xlim',[1,N*M(1).dt],'nextplot','add')
        pf2 = plot(time,SCKS.qU.z{1}(:,1:nD:end),'.'); %residual
        if ~isempty(SCKS.pU.r)
            pf3 = plot(time,SCKS.pU.r.y(:,:),'-'); %true
        end
        %doColor(M(1).YName,pf),doColor(M(1).YName,pf2),doColor(M(1).YName,pf3)
        
        
        title(tit); %drawnow;
        grid(gca,'on')
        axis(gca,'tight')
        set(gca,'box','on','tickdir','out');
        legend([pf(1) pf2(1) pf3(1)],{'Estimate','Residuals','Ori'});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4: Residuals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f4 = spm_figure('Create','Graphics','SCKS Input residual');
        set(f4,'RendererMode','auto','toolbar','figure');
        try
            subplot(3,1,1)
            pf2 = plot(time,SCKS.qU.z{1}(1,1:nD:end),'-'); %residual
            hold on
            h1 = plot(time2,.3*v(1:length(time2))','Color','r','linewidth',1); %PP red more visible
            axis tight
            subplot(3,1,2)
            pf2 = plot(time,SCKS.qU.z{1}(2,1:nD:end),'-'); %residual
            axis tight
            subplot(3,1,3)
            pf2 = plot(time,SCKS.qU.z{1}(3,1:nD:end),'-'); %residual
            axis tight
        end
        if SCKS.DO.save_figures
            nirs_SCKS_print_figures(SCKS,f1,['StateEst']);
            nirs_SCKS_print_figures(SCKS,f2,['InputEst']);
            nirs_SCKS_print_figures(SCKS,f3,['RespEst']);
            nirs_SCKS_print_figures(SCKS,f4,['Resid']);
        end
    else
        disp('SCKS estimation did not work.');
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    disp('Problem displaying SCKS results');
end
% function doColor(name,a)
%
% axis tight
% color='kBRGM';
% mode={'Flow','HbR','HbO','HbT','CMRO'};
%
% if length(a)<=4
%     for i2=1:min(a,length(name))
%         try
%             set(a(i2),'color',color(strmatch(name{i2},mode)))
%         end
%     end
% else
%     set(a(7:end),'linestyle',':')
% end