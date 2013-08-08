function nirs_HDM_display(HDM)
try
    if isfield(HDM,'HDM_OK')
        M=HDM;
        Model_Choice = HDM.O.PhysioModel_Choice;
        U = M.U;
        m = size(M.U.u,2);
        Ep = HDM.Ep;
        Cp = HDM.Cp;
        pE = HDM.pE;
        %subj_id = HDM.subj_id;
        plot_algebraic_CMRO2 = HDM.DO.plot_algebraic_CMRO2;
        
        Finter = spm_figure('GetWin','Interactive');
        header = get(Finter,'Name');
        set(Finter,'Name','Hemodynamic modelling')
        
        % Figure 1:
        % display input parameters
        %--------------------------------------------------------------------------
        t       = [1:M.N]*M.dt;
        Fhdm    = spm_figure;
        set(Fhdm,'name','Hemodynamic Modeling')
        
        subplot(2,2,1)
        [P C] = private_getPC(Model_Choice,Ep,Cp);
        
        [dummy, j] = max(abs(P));
        spm_barh(P,C)
        axis square
        title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
        switch Model_Choice
            case {0,1,3,4} %Buxton-Friston
                set(gca,'Ytick',[1:m],'YTickLabel',U.name,'FontSize',8)
                str = {};
                for i = 1:m
                    str{end + 1} = U.name{i};
                    str{end + 1} = sprintf('mean = %0.2f',P(i));
                    str{end + 1} = '';
                end
                set(gca,'Ytick',[1:m*3]/3 + 1/2,'YTickLabel',str)
            case 2
                m2 = 2*m;
                for m3=1:m
                    temp_labels{m3} = ['CMRO2' U.name{m3}];
                end
                labels2 = [U.name temp_labels];
                set(gca,'Ytick',[1:m2],'YTickLabel',labels2,'FontSize',8)
                str = {};
                for i = 1:m2
                    str{end + 1} = [labels2{i} '*100']; %U.name{i};
                    str{end + 1} = sprintf('mean = %0.2f',100*P(i));
                    str{end + 1} = '';
                end
                set(gca,'Ytick',[1:m2*3]/3 + 1/2,'YTickLabel',str)
            otherwise
        end
        xlabel('relative efficacy per event/sec')
        
        % display hemodynamic parameters
        %---------------------------------------------------------------------------
        subplot(2,2,3)
        switch Model_Choice
            case {0,4} %Buxton-Friston
                P     = Ep(1:5);
                pE    = pE(1:5);
                C     = diag(Cp(1:5,1:5));
                spm_barh(P,C,pE)
                title({ 'hemodynamic parameters'},'FontSize',10)
                set(gca,'Ytick',[1:18]/3 + 1/2)
                set(gca,'YTickLabel',{  'SIGNAL decay',...
                    sprintf('%0.2f per sec',P(1)),'',...
                    'FEEDBACK',...
                    sprintf('%0.2f per sec',P(2)),'',...
                    'TRANSIT TIME',...
                    sprintf('%0.2f seconds',P(3)),'',...
                    'EXPONENT',...
                    sprintf('%0.2f',P(4)),'',...
                    'EXTRACTION',...
                    sprintf('%0.0f %s',P(5)*100,'%'),''},'FontSize',8)
            case 1 %Zheng-Mayhew
                P     = Ep(1:7);            %MODIFIER
                pE    = pE(1:7);            %MODIFIER
                C     = diag(Cp(1:7,1:7));  %MODIFIER
                spm_barh(P,C,pE)
                title({ 'hemodynamic parameters'},'FontSize',10)
                set(gca,'Ytick',[1:18+6]/3 + 1/2)   %MODIFIER? was 18
                set(gca,'YTickLabel',{  'SIGNAL decay',...
                    sprintf('%0.2f per sec',P(1)),'',...
                    'FEEDBACK',...
                    sprintf('%0.2f per sec',P(2)),'',...
                    'TRANSIT TIME',...
                    sprintf('%0.2f seconds',P(3)),'',...
                    'EXPONENT',...
                    sprintf('%0.2f',P(4)),'',...
                    'EXTRACTION',...
                    sprintf('%0.0f %s',P(5)*100,'%'),'',...
                    'VASCULAR TONE',...                                %MODIFIER
                    sprintf('%0.2f *10 seconds',P(6)),'',...              %MODIFIER
                    'GAIN PARAMETER',...                               %MODIFIER
                    sprintf('%0.3f *10 seconds',P(7)),''},'FontSize',8)   %MODIFIER
            case 2 %Huppert1
                P     = Ep(1:11);
                P(7) = P(7)/10;
                pE    = pE(1:11);
                pE(7) = pE(7)/10;
                C     = diag(Cp(1:11,1:11));
                C(7) = C(7)/10;
                spm_barh(P,C,pE)
                title({ 'hemodynamic parameters'},'FontSize',10)
                set(gca,'Ytick',[1:18+18]/3 + 1/2)
                set(gca,'YTickLabel',{  'Flow decay (ksr)',...
                    sprintf('%0.2f per sec',P(1)),'',...
                    'Flow feedback (kr)',...
                    sprintf('%0.2f per sec',P(2)),'',...
                     'CMRO2 decay (ksm)',...
                    sprintf('%0.2f per sec',P(3)),'',...
                    'CMRO2 feedback (km)',...
                    sprintf('%0.2f per sec',P(4)),'',...
                    'Pressure (P0r)',...
                    sprintf('%0.2f seconds',P(5)),'',...
                    'Veinous volume (Vw0)',...
                    sprintf('%0.4f',P(6)),'',...
                    'beta/10',...
                    sprintf('%0.2f',P(7)/10),'',...
                    'delta',...
                    sprintf('%0.2f',P(8)),'',...
                    'HGB',...
                    sprintf('%0.2f per sec',P(9)),'',...
                    'Ra0r',...
                    sprintf('%0.2f per sec',P(10)),'',...
                    'M0',...
                    sprintf('%0.2f per sec',P(11)),''},'FontSize',8)
        end
        
        % get display state kernels (i.e. state dynamics)
        %==========================================================================
        H1 = M.H1;
        j = 1; %PP ?
        subplot(3,2,2)
        if plot_algebraic_CMRO2
            tmp_H1 = exp(H1(:,:,j));
            %Algebraic relation for m = CMRO2, in arbitrary units
            %m = f * HbR /HbT; assuming gamma_R and gamma_T = 1; f: flow
            tmp_ma = tmp_H1(:,2) .* tmp_H1(:,4) ./ tmp_H1(:,3);
            tmp_H1 = [tmp_H1 tmp_ma];
            plot(t,tmp_H1)
            axis square
            title({['1st order kernels for ' U.name{j}];...
                'state variables'},'FontSize',9)
            ylabel('normalized values')
            switch Model_Choice
                case {0,4} %Buxton-Friston
                    legend('s','f','v','q','ma',0);
                case 1 %Zheng-Mayhew
                    legend('s','f','v','q','w','ma',0); %MODIFIER
                case 2 %Huppert1
                    %legend('s','f','v','q','s2','m','Ct','Cv','ma',0)
                    legend('s','f','s2','f2','v','t','p','ma',0)
                case 3 %Buxton-Friston part 1
                    legend('s','f',0);
                otherwise
                    legend('s','f','v','q','ma',0);
            end
            grid on
        else
            plot(t,exp(H1(:,:,j)))
            axis square
            title({['1st order kernels for ' U.name{j}];...
                'state variables'},'FontSize',9)
            ylabel('normalized values')
            switch Model_Choice
                case {0,4} %Buxton-Friston
                    legend('s','f','v','q',0);
                case 1 %Zheng-Mayhew
                    legend('s','f','v','q','w', 0); %MODIFIER
                case 2 %Huppert1
                    %legend('s','f','v','q','s2','m','Ct','Cv',0)
                    legend('s','f','s2','f2','v','t','p',0)
                case 3 %Buxton-Friston part 1
                    legend('s','f',0);
                otherwise
                    legend('s','f','v','q',0);
            end
            grid on
        end
        
        % display output kernels (i.e. hemodynamic response)
        %--------------------------------------------------------------------------
        K1 = M.K1;
        K2 = M.K2;
        subplot(3,2,4)
        plot(t,K1(:,:,j))
        axis square        
        switch Model_Choice
            case {0,1,3,4}
                title({ '1st order kernel';'output: HbT, HbR'},'FontSize',9)
            case 2
                title({ '1st order kernel';'output: HbT, HbR, HbO'},'FontSize',9)
        end        
        ylabel('normalized measure responses')        
        grid on        
        subplot(3,2,6)
        axis square
%         imagesc(t,t,K2(:,:,1,j,j))
%         title({ '2nd order kernel';'output: flow'},'FontSize',9)
        imagesc(t,t,K2(:,:,2,j,j))
        title({ '2nd order kernel';'output: HbR'},'FontSize',9)        
        xlabel({'time {seconds} for'; U.name{j}})
        grid on
                        
        if HDM.DO.save_figures
            nirs_HDM_print_figures(HDM,Fhdm,['HDM ']);
        end
    else
        disp('HDM estimation did not work.');
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    disp('Problem displaying HDM results');
end

function [P C] = private_getPC(Model_Choice,Ep,Cp)
switch Model_Choice
    case 0 %Buxton-Friston
        P     = Ep(6:end);
        C     = diag(Cp(6:end,6:end));
    case 1 %Zheng-Mayhew
        P     = Ep(8:end);                  %MODIFIER
        C     = diag(Cp(8:end,8:end));      %MODIFIER
    case 2 %Huppert1
        P     = Ep(12:end);                 %MODIFIER
        C     = diag(Cp(12:end,12:end));    %MODIFIER
end