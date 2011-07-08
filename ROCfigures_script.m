function ROCfigures_script
Z = load('ROC.mat');
[nJob nSubj] = size(Z.TPu1);
if isfield(Z,'TPu1u') || isfield(Z,'TPu1Ou')
    compute_LU = 1;
else
    compute_LU = 0;
end
if isfield(Z,'TPu1O')
    compute_OR = 1; 
else
    compute_OR = 0;
end
if isfield(Z,'TPu1A')
    all_channels = 1; 
else
    all_channels = 0;
end
line_choice = 1;
figure_choice =10;
%1+3 or 4 plots
plot4 = 0; %Boolean
plotUL = 1;
%0: generic testing - 1 data set, 1st Volterra
%14: using all channels, just ROC plot
%1: figure for 1st Volterra
%12: V1 figure for different number of spikes
%13: V2 figure for different number of spikes
%10: figure for 1st Volterra - canonical HRF  - Also for HbO + HbR with
%plotUL = 1
%2: figure for 1st of 2 Volterra estimation
%3: figure for 2nd of 2 Volterra estimation
%21: WMDL 1st Volterra
%31: WMDL 2nd Volterra
%4: Compare negative 2nd to positive 2nd Volterra: figure for 2st of 2 Volterra estimation
%5: Compare negative 2nd to positive 2nd Volterra: figure for 1nd of 2 Volterra estimation
%6: Compare Gamma to Canonical HRF - 1st Volterra
%61: Compare Neg and Pos
%62: Compare V and noV
%7: Compare Gamma to Canonical HRF - 2nd Volterra
%8: Canonical HRF - figure for 1st of 2 Volterra estimation
%9: Canonical HRF - figure for 2nd of 2 Volterra estimation

%ROC plots
%Compute ROC area
%Compute area under ROC curve
for Jidx=1:nJob
    for Idx=1:nSubj
        %1st Volterra
        meanROC=(Z.TPu1{Jidx,Idx}(1:end-1)+Z.TPu1{Jidx,Idx}(2:end))/2;
        ROCarea1{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu1{Jidx,Idx}));
        %2nd Volterra
        meanROC=(Z.TPu2{Jidx,Idx}(1:end-1)+Z.TPu2{Jidx,Idx}(2:end))/2;
        ROCarea2{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu2{Jidx,Idx}));
        if compute_OR
            %HbO
            meanROC=(Z.TPu1O{Jidx,Idx}(1:end-1)+Z.TPu1O{Jidx,Idx}(2:end))/2;
            ROCareaO{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu1O{Jidx,Idx}));
            %HbR
            meanROC=(Z.TPu1R{Jidx,Idx}(1:end-1)+Z.TPu1R{Jidx,Idx}(2:end))/2;
            ROCareaR{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu1R{Jidx,Idx}));
            %HbO
            meanROC=(Z.TPu2O{Jidx,Idx}(1:end-1)+Z.TPu2O{Jidx,Idx}(2:end))/2;
            ROCareaO2{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu2O{Jidx,Idx}));
            %HbR
            meanROC=(Z.TPu2R{Jidx,Idx}(1:end-1)+Z.TPu2R{Jidx,Idx}(2:end))/2;
            ROCareaR2{Jidx,Idx} = -sum(meanROC.*diff(Z.FPu2R{Jidx,Idx}));
        end
        if all_channels
            Nch = size(Z.FPu1A{Jidx,Idx},1);
            ROCarea1A = zeros(1,Nch);
            ROCarea2A = zeros(1,Nch);
            
            for i=1:Nch
                meanROC=(Z.TPu1A{Jidx,Idx}(i,1:end-1)+Z.TPu1A{Jidx,Idx}(i,2:end))/2;
                ROCarea1A(i) =  -sum(meanROC.*diff(Z.FPu1A{Jidx,Idx}(i,:)));
                meanROC=(Z.TPu2A{Jidx,Idx}(i,1:end-1)+Z.TPu2A{Jidx,Idx}(i,2:end))/2;
                ROCarea2A(i) =  -sum(meanROC.*diff(Z.FPu2A{Jidx,Idx}(i,:)));
            end
        end
    end
end

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultTextFontSize',16)




% % set(gcf,'PaperUnits','centimeters')
% % %This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
% % xSize = 9; ySize = 7;
% % %These are my size variables, width of 8 and a height of 12, will be used a lot later.
% % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % %Additional coordinates to center the figure on A4-paper
% % set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% % %This command sets the position and size of the figure on the paper to the desired values.
% % set(gcf,'Position',[X Y xSize*50 ySize*50])

switch figure_choice
    case 0
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

             
        legend([h1(1)],['Generic, area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%']); hold on
        %set(gca,'xlim',[0,1.01]); hold on
        %set(gca,'ylim',[0,1.01]); hold on
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'Generic'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 14
        figure; hold on
        %if plot4
        %    subplot(2,2,1)
            title('ROC curves'); hold on
        %end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx = 1;
        %h = [];
        for nch=1:Nch
            if nch <= Nch/2
                h1 = plot(Z.FPu1A{Jidx,Idx}(nch,:),Z.TPu1A{Jidx,Idx}(nch,:),':r'); hold on
            else
                h2 = plot(Z.FPu1A{Jidx,Idx}(nch,:),Z.TPu1A{Jidx,Idx}(nch,:),'-b'); hold on          
            end
        end
        legend([h1 h2],'HbO','HbR')
        %if plot4
        %    subplot(2,2,2);
        %else
            figure; 
        %    subplot(1,3,1);
        %end
    
        title('ROC Area'); hold on
          
        Idx=1;
        Jidx=1;
%          tmp = [ROCarea1A' ROCarea2A']; 
%         BPL = {'V1','V2'};
        tmp = [ROCarea1A(1:Nch/2)' ROCarea1A(1+Nch/2:Nch)' ...
               ROCarea2A(1:Nch/2)' ROCarea2A(1+Nch/2:Nch)']; 
        BPL = {'O1','R1','O2','R2'};
        boxplot(tmp,'notch','on','labels',BPL);
        
    case 1
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
            
        legend([h1(1) h2(1) h3(1)],['(0.75,0), area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
            ['(0.50,0), area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
            ['(0.25,0), area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%']); hold on
        %set(gca,'xlim',[0,1.01]); hold on
        %set(gca,'ylim',[0,1.01]); hold on
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'(0.75,0)','(0.5,0)','(0.25,0)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 11

        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
             
%         Jidx = 4;
%         h4 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':g',...
%                  Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'sg'); hold on
%         Jidx = 5;
%         h5 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':y',...
%                  Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'dy'); hold on
%             
%         legend([h1(1) h2(1) h3(1) h4(1) h5(1)],...
%             ['720, area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
%             ['360, area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
%             ['180, area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%'],...
%             ['90, area=' num2str(ROCarea1{4,1}*100,'%3.1f') '%'],...
%             ['45, area=' num2str(ROCarea1{5,1}*100,'%3.1f') '%']); hold on
        
        legend([h1(1) h2(1) h3(1)],...
            ['360, V1 area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['180, V1 area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['90, V1 area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']); hold on
        %set(gca,'xlim',[0,1.01]); hold on
        %set(gca,'ylim',[0,1.01]); hold on
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        %BPL = {'720','360','180','90','45'};
        BPL = {'360','180','90'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 12

        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
             
%         Jidx = 4;
%         h4 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':g',...
%                  Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'sg'); hold on
%         Jidx = 5;
%         h5 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':y',...
%                  Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'dy'); hold on
%             
%         legend([h1(1) h2(1) h3(1) h4(1) h5(1)],...
%             ['720, area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
%             ['360, area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
%             ['180, area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%'],...
%             ['90, area=' num2str(ROCarea1{4,1}*100,'%3.1f') '%'],...
%             ['45, area=' num2str(ROCarea1{5,1}*100,'%3.1f') '%']); hold on
        
        legend([h1(1) h2(1) h3(1)],...
            ['360, V1 area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['180, V1 area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['90, V1 area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']); hold on
        %set(gca,'xlim',[0,1.01]); hold on
        %set(gca,'ylim',[0,1.01]); hold on
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        %BPL = {'720','360','180','90','45'};
        BPL = {'360','180','90'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 10
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                    Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
        else
            if ~compute_OR
                Jidx= 1;
                h1 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-r'); hold on
                herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-r'); hold on
                Jidx= 2;
                h2 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'--k'); hold on
                herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'--k'); hold on
                Jidx= 3;
                h3 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-.b'); hold on
                herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-.b'); hold on
            else
                Jidx= 1;
                h1 = errorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.TPu1Ol{Jidx,Idx},Z.TPu1Ou{Jidx,Idx},':r'); hold on
                herrorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.FPu1Ol{Jidx,Idx},Z.FPu1Ou{Jidx,Idx},':r'); hold on
                Jidx= 2;
                h2 = errorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.TPu1Ol{Jidx,Idx},Z.TPu1Ou{Jidx,Idx},':r'); hold on
                herrorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.FPu1Ol{Jidx,Idx},Z.FPu1Ou{Jidx,Idx},':r'); hold on
                Jidx= 3;
                h3 = errorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.TPu1Ol{Jidx,Idx},Z.TPu1Ou{Jidx,Idx},':r'); hold on
                herrorbar(Z.FPu1O{Jidx,Idx},Z.TPu1O{Jidx,Idx},Z.FPu1Ol{Jidx,Idx},Z.FPu1Ou{Jidx,Idx},':r'); hold on
                Jidx= 1;
                h4 = errorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.TPu1Rl{Jidx,Idx},Z.TPu1Ru{Jidx,Idx},'-b'); hold on
                herrorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.FPu1Rl{Jidx,Idx},Z.FPu1Ru{Jidx,Idx},'-b'); hold on
                Jidx= 2;
                h5 = errorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.TPu1Rl{Jidx,Idx},Z.TPu1Ru{Jidx,Idx},'-b'); hold on
                herrorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.FPu1Rl{Jidx,Idx},Z.FPu1Ru{Jidx,Idx},'-b'); hold on
                Jidx= 3;
                h6 = errorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.TPu1Rl{Jidx,Idx},Z.TPu1Ru{Jidx,Idx},'-b'); hold on
                herrorbar(Z.FPu1R{Jidx,Idx},Z.TPu1R{Jidx,Idx},Z.FPu1Rl{Jidx,Idx},Z.FPu1Ru{Jidx,Idx},'-b'); hold on   
            end
        end
        if ~compute_OR
            legend([h1(1) h2(1) h3(1)],['100%, area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
                ['75%, area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
                ['50%, area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%']); hold on
        else
%             legend([h1(1) h4(1) h2(1) h5(1) h3(1) h6(1)],['HbO (3, -3), area=' num2str(ROCareaO{1,1}*100,'%3.1f') '%'],...
%                 ['HbR (3, -3), area=' num2str(ROCareaR{1,1}*100,'%3.1f') '%'],...
%                 ['HbO (2, -2) area=' num2str(ROCareaO{2,1}*100,'%3.1f') '%'],...
%                 ['HbR (2, -2), area=' num2str(ROCareaR{2,1}*100,'%3.1f') '%'],...
%                 ['HbO (1, -1), area=' num2str(ROCareaO{3,1}*100,'%3.1f') '%'],...
%                 ['HbR (1, -1) area=' num2str(ROCareaR{3,1}*100,'%3.1f') '%']); hold on
legend([h1(1) h4(1) h2(1) h5(1) h3(1) h6(1)],...
                ['HbO (3, -3), V1: ' num2str(ROCareaO{1,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaO2{1,1}*100,'%3.1f') '%'],...
                ['HbR (3, -3), V1: ' num2str(ROCareaR{1,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaR2{1,1}*100,'%3.1f') '%'],...
                ['HbO (2, -2), V1: ' num2str(ROCareaO{2,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaO2{2,1}*100,'%3.1f') '%'],...
                ['HbR (2, -2), V1: ' num2str(ROCareaR{2,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaR2{2,1}*100,'%3.1f') '%'],...
                ['HbO (1, -1), V1: ' num2str(ROCareaO{3,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaO2{3,1}*100,'%3.1f') '%'],...
                ['HbR (1, -1), V1: ' num2str(ROCareaR{3,1}*100,'%3.1f') '%, ' 'V2: ' num2str(ROCareaR2{3,1}*100,'%3.1f') '%']); hold on
            title('ROC curves for HbO and HbR'); hold on
        end
        set(gca,'xlim',[0,1.01]); hold on
        set(gca,'ylim',[0,1.01]); hold on
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'(1,0)','(0.75,0)','(0.5,0)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 2
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
        
        legend([h1(1) h2(1) h3(1)],['(3,-3), area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
            ['(2,-2), area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
            ['(1,-1), area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'(3,-3)','(2,-2)','(1,-1)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 21
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
        
        legend([h1(1) h2(1) h3(1)], ['B,   V1 area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['W,   V1 area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['B+W,  V1 area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'B','W','B+W'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 3
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on    
        
        legend([h1(1) h2(1) h3(1)],['(3,-3), area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['(2,-2), area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['(1,-1), area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'(3,-3)','(2,-2)','(1,-1)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        set(gca,'ylim',[-4,2]);
        hold off
    case 31
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on    
        
        legend([h1(1) h2(1) h3(1)],['B, area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['W, area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['B+W, area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'B','W','B+W'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        set(gca,'ylim',[-4,2]);
        hold off     
     
     case 4
        %for negative vs positive 2nd Volterra amplitudes
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

        Jidx = 2;
        h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                 Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on

        
        legend([h1(1) h2(1) h3(1)],['(4,-4), area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
            ['(4,0), area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
            ['(4,4), area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'(4,-4)','(4,0)','(4,4)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
        
        %%%%%%%%%%%
    case 5
        %linespec = fn_linespec(2);
        %4 plots combined
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        Jidx= 1;
        h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

        %Jidx = 2;
        %h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
        %        Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

        Jidx = 3;
        h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on    
        
        legend([h1(1) h3(1)],['(4,-4), area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...           
            ['(4,4), area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'(4,-4)','(4,0)','(4,4)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        set(gca,'ylim',[-5,5]);
        hold off
    case 6
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                    Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
            Jidx = 4;
            h4 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':g',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'sg'); hold on
        else
            Jidx= 1;
            h1 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-.b'); hold on
            Jidx= 4;
            h4 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},':g'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1) h4(1)],...
            ['CC, V1 area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['GG, V1 area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['CG, V1 area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%'],...
            ['GC, V1 area=' num2str(ROCarea1{4,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{4,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'CC','GG','CG','GC'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 61
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                    Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on
            Jidx = 4;
            h4 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},':g',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'sg'); hold on
        else
            Jidx= 1;
            h1 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-.b'); hold on
            Jidx= 4;
            h4 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},':g'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1) h4(1)],...
            ['(2,-2), V1 area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['(2,-1), V1 area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['(2,1),  V1 area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%'],...
            ['(2,2),  V1 area=' num2str(ROCarea1{4,1}*100,'%3.1f') '%' ', V2 area=' num2str(ROCarea2{4,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'(2,-2)','(2,-1)','(2,1)','(2,2)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
%         title('T-Statistic'); hold on
%         tmp = [];    
%         for Idx=1:nSubj    
%             for Jidx= 1:nJob
%                tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
%             end   
%         end
%         boxplot(tmp,'notch','on','labels',BPL);  
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        %set(gca,'ylim',[-4,2]);
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
     case 62
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
             Jidx= 1;
            h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                    Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
            Jidx = 4;
            h4 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':g',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'sg'); hold on
        else
            Jidx= 1;
            h1 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-.b'); hold on
            Jidx= 4;
            h4 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},':g'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1) h4(1)],...
            ['VV, V1: ' num2str(ROCarea1{1,1}*100,'%3.1f') '%' ', V2: ' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['VN, V1: ' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
            ['NV,  V1: ' num2str(ROCarea1{3,1}*100,'%3.1f') '%' ', V2: ' num2str(ROCarea2{3,1}*100,'%3.1f') '%'],...
            ['NN,  V1: ' num2str(ROCarea1{4,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('V1 Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'VV','VN','NV','NN'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);  
%        title('V2 Estimated Amplitude'); hold on
%         tmp = [];    
%         for Idx=1:nSubj    
%             for Jidx= 1:nJob
%                tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
%             end   
%         end
%         BPL = {'(2,-2)','(2,-1)','(2,1)','(2,2)'};
%         boxplot(tmp,'notch','on','labels',BPL);
%                 
%         title('Ratio of 2nd to 1st amplitudes'); hold on
%         tmp = [];    
%         for Idx=1:nSubj    
%             for Jidx= 1:nJob
%                tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
%             end   
%         end
%         boxplot(tmp,'notch','on','labels',BPL); 
        %set(gca,'ylim',[-4,2]);
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 7
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                    Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on
            Jidx = 4;
            h4 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},':g',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'sg'); hold on
        else
            Jidx= 1;
            h1 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'-.b'); hold on
            Jidx= 4;
            h4 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},':g'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1) h4(1)],['CC, area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['GG, area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['CG, area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%'],...
            ['GC, area=' num2str(ROCarea2{4,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'CC','GG','CG','GC'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        set(gca,'ylim',[-3,1]);
        hold off
        %%%
    case 8
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-r',...
                    Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'--k',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'-.b',...
                     Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'xb'); hold on
%             Jidx = 4;
%             h4 = plot(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},':g',...
%                      Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},'sg'); hold on
        else
            Jidx= 1;
            h1 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},'-.b'); hold on
%             Jidx= 4;
%             h4 = errorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.TPu1l{Jidx,Idx},Z.TPu1u{Jidx,Idx},':g'); hold on
%             herrorbar(Z.FPu1{Jidx,Idx},Z.TPu1{Jidx,Idx},Z.FPu1l{Jidx,Idx},Z.FPu1u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1)],['(1,-1), area=' num2str(ROCarea1{1,1}*100,'%3.1f') '%'],...
            ['(2,-2), area=' num2str(ROCarea1{2,1}*100,'%3.1f') '%'],...
            ['(3,-3), area=' num2str(ROCarea1{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,1,:),3))];          
            end   
        end
        BPL = {'(1,-1)','(2,-2)','(3,-3)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,1,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Signal-to-Noise Ratio'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp Z.T{Jidx,Idx}.SNR(:)];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        hold off
    case 9
        %linespec = fn_linespec(2);
        %4 plots combined
        %1st Volterra
        figure; hold on
        if plot4
            subplot(2,2,1)
            title('ROC curves'); hold on
        end 
        box on; hold on
        xlabel('1-Specificity');
        ylabel('Sensitivity'); hold on
        if ~plotUL
            Jidx= 1;
            h1 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-r',...
                    Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'or'); hold on

            Jidx = 2;
            h2 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'--k',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'+k'); hold on

            Jidx = 3;
            h3 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'-.b',...
                     Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'xb'); hold on
%             Jidx = 4;
%             h4 = plot(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},':g',...
%                      Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},'sg'); hold on
         else
            Jidx= 1;
            h1 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'-r'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'-r'); hold on
            Jidx= 2;
            h2 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'--k'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'--k'); hold on
            Jidx= 3;
            h3 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},'-.b'); hold on
            herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},'-.b'); hold on
%             Jidx= 4;
%             h4 = errorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.TPu2l{Jidx,Idx},Z.TPu2u{Jidx,Idx},':g'); hold on
%             herrorbar(Z.FPu2{Jidx,Idx},Z.TPu2{Jidx,Idx},Z.FPu2l{Jidx,Idx},Z.FPu2u{Jidx,Idx},':g'); hold on
            set(gca,'xlim',[0,1.01]);
            set(gca,'ylim',[0,1.01]);
        end
      
        legend([h1(1) h2(1) h3(1)],['(1,-1), area=' num2str(ROCarea2{1,1}*100,'%3.1f') '%'],...
            ['(2,-2), area=' num2str(ROCarea2{2,1}*100,'%3.1f') '%'],...
            ['(3,-3), area=' num2str(ROCarea2{3,1}*100,'%3.1f') '%']);
        if plot4
            subplot(2,2,2);
        else
            figure; 
            subplot(1,3,1);
        end
        title('Estimated Amplitude'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,2,:),3))];          
            end   
        end
        BPL = {'(1,-1)','(2,-2)','(3,-3)'};
        boxplot(tmp,'notch','on','labels',BPL);
                
        if plot4
            subplot(2,2,3);
        else
            subplot(1,3,2);
        end
        title('T-Statistic'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.t(:,2,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL);    
        if plot4
            subplot(2,2,4);
        else
            subplot(1,3,3);
        end
        title('Ratio of 2nd to 1st amplitudes'); hold on
        tmp = [];    
        for Idx=1:nSubj    
            for Jidx= 1:nJob
               tmp = [tmp squeeze(mean(Z.T{Jidx,Idx}.RE(:,3,:),3))];          
            end   
        end
        boxplot(tmp,'notch','on','labels',BPL); 
        set(gca,'ylim',[-3,1]);
        hold off
end
end

function scatspec = fn_scatspec
    scatspec{1,1} = 'ob';
    scatspec{2,1} = '.k';
    scatspec{3,1} = '+r';
    scatspec{4,1} = 'sg';
end

function linespec = fn_linespec(line_choice)
switch line_choice
    case 1
        linespec{1,1} = '.';
        linespec{2,1} = '+';
        linespec{3,1} = 'o';
        linespec{4,1} = '^';
        linespec{5,1} = 's';
        linespec{6,1} = 'p';
        linespec{7,1} = 'd';
        linespec{8,1} = '*';
        linespec{9,1} = 'x';
        linespec{10,1} = 'h';
    case 2
        linespec{1,1} = '-b';
        linespec{2,1} = '--r';
        linespec{3,1} = ':k';
        linespec{4,1} = '-.g';
        linespec{5,1} = '-k';
        linespec{6,1} = '--g';
        linespec{7,1} = ':k';
        linespec{8,1} = '-.b';
    case 3        
        linespec{1,1} = '-';
        linespec{2,1} = ':';
        linespec{3,1} = '--';
        linespec{4,1} = '-.';
        linespec{5,1} = '-';
        linespec{6,1} = '--';
        linespec{7,1} = ':';
        linespec{8,1} = '-.';
    case 4
        linespec{1,1} = '.';
        linespec{2,1} = '+';
        linespec{3,1} = 'o';
        linespec{1,2} = '*';
        linespec{2,2} = 'x';
        linespec{3,2} = 's';
        linespec{1,3} = 'd';
        linespec{2,3} = '^';
        linespec{3,3} = 'v';
        linespec{1,4} = '<';
        linespec{2,4} = '>';
        linespec{3,4} = 'p';
        linespec{1,5} = 'h';
end 
end
% 
% %Tstat
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(squeeze(mean(T{Jidx,Idx}.t(:,1,:),1)),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% %Correlation
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.RE(:,3,:),3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,max(T{Jidx,Idx}.RE(:,3,:),[],3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.t(:,1,:),3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.t(:,2,:),3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,min(T{Jidx,Idx}.t(:,2,:),[],3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,median(T{Jidx,Idx}.t(:,2,:),3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,prctile(T{Jidx,Idx}.t(:,2,:),25,3),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,mean(T{Jidx,Idx}.SNR,2),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% %[b,bint,r,rint,stats] = regress(mean(T{Jidx,Idx}.SNR,2),T{Jidx,Idx}.CORR)
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         scatter(T{Jidx,Idx}.CORR,10.^(mean(T{Jidx,Idx}.SNR,2)/10),scatspec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% 
% %1st Volterra
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(FPb1{Jidx,Idx},TPb1{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %2nd Volterra
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(FPb2{Jidx,Idx},TPb2{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %uncorrected only
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},linespec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% 
% %1st and 2nd Volterra grouped with subplot - with error bars
% figure;
% subplot(2,1,1);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%          %add error bars
%         errorbar(FPu1{Jidx,Idx},TPu1{Jidx,Idx},TPu1l{Jidx,Idx},TPu1u{Jidx,Idx},linespec{Jidx,Idx}); hold on
%         herrorbar(FPu1{Jidx,Idx},TPu1{Jidx,Idx},FPu1l{Jidx,Idx},FPu1u{Jidx,Idx},linespec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% subplot(2,1,2);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%          errorbar(FPu2{Jidx,Idx},TPu2{Jidx,Idx},TPu2l{Jidx,Idx},TPu2u{Jidx,Idx},linespec{Jidx,Idx}); hold on
%          herrorbar(FPu2{Jidx,Idx},TPu2{Jidx,Idx},FPu2l{Jidx,Idx},FPu2u{Jidx,Idx},linespec{Jidx,Idx}); hold on
%     end   
% end
% hold off
% 
% %1st and 2nd Volterra grouped with subplot
% figure;
% subplot(2,1,1);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(FPb1{Jidx,Idx},TPb1{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% subplot(2,1,2);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(FPb2{Jidx,Idx},TPb2{Jidx,Idx},[linespec{Jidx,Idx} 'r'],...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% 
% %figure for SNR
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         m = median(T{Jidx,Idx}.SNR(:));
%         errorbar(median(T{Jidx,Idx}.a2),m,m-prctile(T{Jidx,Idx}.SNR(:),25),prctile(T{Jidx,Idx}.SNR(:),75)-m); hold on
%     end   
% end
% hold off
% %better as a boxplot?
% figure;
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp T{Jidx,Idx}.SNR(:)];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %COMPARE HbO and HbR
% %1st and 2nd Volterra grouped with subplot - with error bars
% figure;
% subplot(2,1,1);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%          %add error bars
%         errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%         herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% subplot(2,1,2);
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%          errorbar(FPu2O{Jidx,Idx},TPu2O{Jidx,Idx},TPu2Ol{Jidx,Idx},TPu2Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%          herrorbar(FPu2O{Jidx,Idx},TPu2O{Jidx,Idx},FPu2Ol{Jidx,Idx},FPu2Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%          errorbar(FPu2R{Jidx,Idx},TPu2R{Jidx,Idx},TPu2Rl{Jidx,Idx},TPu2Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%          herrorbar(FPu2R{Jidx,Idx},TPu2R{Jidx,Idx},FPu2Rl{Jidx,Idx},FPu2Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %Custom:
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%          %add error bars
%         errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%         herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% 
% %Custom:
% figure;
% for Idx=1:nSubj    
%     for Jidx=[2 4]; % 1:nJob
%          %add error bars
%         errorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},TPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         herrorbar(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},FPu1Ol{Jidx,Idx},FPu1Ou{Jidx,Idx},[linespec{Jidx,Idx} 'r']); hold on
%         errorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},TPu1Rl{Jidx,Idx},TPu1Ru{Jidx,Idx},[linespec{Jidx+2,Idx} 'b']); hold on
%         herrorbar(FPu1R{Jidx,Idx},TPu1R{Jidx,Idx},FPu1Rl{Jidx,Idx},FPu1Ru{Jidx,Idx},[linespec{Jidx+2,Idx} 'b']); hold on
%     end   
% end
% hold off
% 
% %SNR vs t-stat
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         tmp= T{Jidx,Idx}.t(:,1,1:half_chan_len);
%         scatter(T{Jidx,Idx}.SNR(:),tmp(:));
%         scatter(mean(T{Jidx,Idx}.SNR,2),squeeze(mean(tmp,3)));
%     end
% end
% addpath('J:\NIRS_nonlinear');
% %careful: order of series is the opposite!!!
% XMatrix1 = [FPu1R{2,1};FPu1O{2,1};FPu1R{1,1};FPu1O{1,1};]';
% YMatrix1 = [TPu1R{2,1};TPu1O{2,1};TPu1R{1,1};TPu1O{1,1}]';
% LMatrix1 = [TPu1Rl{2,1};TPu1Ol{2,1};TPu1Rl{1,1};TPu1Ol{1,1}]';
% UMatrix1 = [TPu1Ru{2,1};TPu1Ou{2,1};TPu1Ru{1,1};TPu1Ou{1,1}]';        
% Figure3_createfigure(XMatrix1, YMatrix1, LMatrix1, UMatrix1);
% %Figure3_createfigure(FPu1O{Jidx,Idx},TPu1O{Jidx,Idx},TPu1Ol{Jidx,Idx},
% %TPu1Ou{Jidx,Idx}
% 
% %1st Volterra
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(mean(FPb1{Jidx,Idx}),mean(TPb1{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
%              mean(FPu1{Jidx,Idx}),mean(TPu1{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %2nd Volterra
% figure;
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(mean(FPb2{Jidx,Idx}),mean(TPb2{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
%              mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %1st and 2nd Volterra grouped as subplots
% figure;
% subplot(1,2,1);
% %1st Volterra
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(mean(FPb1{Jidx,Idx}),mean(TPb1{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
%              mean(FPu1{Jidx,Idx}),mean(TPu1{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% subplot(1,2,2);
% %2nd Volterra
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(mean(FPb2{Jidx,Idx}),mean(TPb2{Jidx,Idx}),[linespec{Jidx,Idx} 'r'],...
%              mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% %Uncorrected only, no colors
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%         plot(mean(FPu2{Jidx,Idx}),mean(TPu2{Jidx,Idx}),[linespec{Jidx,Idx} 'b']); hold on
%     end   
% end
% hold off
% 
% %Relative error
% %1st Volterra
% figure;
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
%     end   
% end
% boxplot(tmp); 
% %2nd Volterra
% figure;
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %Relative error grouped as subplots
% %1st Volterra
% figure;
% subplot(1,2,1);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %2nd Volterra
% subplot(1,2,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %Relative error for 1st Volterra and ratio of 2nd to 1st Volterra for
% %second plot
% figure;
% subplot(1,2,1);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %ratio of 2nd to 1st estimated Volterra less simulated
% subplot(1,2,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,3,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %absolute value
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.RE(:,3,:),3)))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %t-stat
% %1st Volterra
% figure;
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %2nd Volterra
% figure;
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %absolute value
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.t(:,2,:),3)))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %1st and 2nd Volterra grouped as subplots
% %1st Volterra
% figure;
% subplot(1,2,1);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %2nd Volterra
% subplot(1,2,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% 
% %max(reshape(cell2mat(TPu),quarter_chan_len,[]))
% %br = squeeze(beta(:,2,(quarter_chan_len+1):half_chan_len))./squeeze(beta(:,1,(quarter_chan_len+1):half_chan_len));
% %
% 
% %figures for box-whisker plot
% %3 plots combined
% %1st Volterra
% figure;
% subplot(1,3,1)
% %for Idx=1:nSubj    
%     Jidx= 1;
%         plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-b',...
%             FPu1{Jidx,Idx},TPu1{Jidx,Idx},'ob'); hold on
% 
%       Jidx = 2;
%        plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'--r',...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
%             Jidx = 3;
%        plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.g',...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+g'); hold on
% 
%          Jidx = 4;
%        plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},':k',...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+k'); hold on
%              Jidx = 5;
%        plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},'-.r',...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+r'); hold on
% 
%          Jidx = 6;
%        plot(FPu1{Jidx,Idx},TPu1{Jidx,Idx},':b',...
%              FPu1{Jidx,Idx},TPu1{Jidx,Idx},'+b'); hold on
% %   end   
% %end
% hold off   
% subplot(1,3,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% subplot(1,3,3);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,1,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% 
% %4 plots combined
% %2nd Volterra
% figure;
% subplot(2,2,1)
% %for Idx=1:nSubj    
%     Jidx= 1;
%         plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-b',...
%             FPu2{Jidx,Idx},TPu2{Jidx,Idx},'ob'); hold on
% 
%       Jidx = 2;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--r',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
% 
%          Jidx = 3;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.g',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+g'); hold on
% 
%          Jidx = 4;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--k',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+k'); hold on
%           Jidx = 5;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.r',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
% 
%          Jidx = 6;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':b',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+b'); hold on
% 
% 
% %   end   
% %end
% hold off   
% subplot(2,2,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% subplot(2,2,3);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(abs(mean(T{Jidx,Idx}.t(:,2,:),3)))];          
%     end   
% end
%  hold off   
%  boxplot(tmp,'notch','on'); 
% subplot(2,2,4);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,3,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% %figures for box-whisker plot
% %3 plots combined
% %2nd Volterra
% figure;
% subplot(1,3,1)
% %for Idx=1:nSubj    
%     Jidx= 1;
%         plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-b',...
%             FPu2{Jidx,Idx},TPu2{Jidx,Idx},'ob'); hold on
% 
%       Jidx = 2;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'--r',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
%             Jidx = 3;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.g',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+g'); hold on
% 
%          Jidx = 4;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':k',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+k'); hold on
%             Jidx = 5;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},'-.r',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+r'); hold on
% 
%          Jidx = 6;
%        plot(FPu2{Jidx,Idx},TPu2{Jidx,Idx},':b',...
%              FPu2{Jidx,Idx},TPu2{Jidx,Idx},'+b'); hold on
% %   end   
% %end
% hold off   
% subplot(1,3,2);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.RE(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 
% subplot(1,3,3);
% tmp = [];    
% for Idx=1:nSubj    
%     for Jidx= 1:nJob
%        tmp = [tmp squeeze(mean(T{Jidx,Idx}.t(:,2,:),3))];          
%     end   
% end
% boxplot(tmp,'notch','on'); 