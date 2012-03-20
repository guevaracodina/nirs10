function H = nirs_copy_figure(H,DF,CF,c1,hb,Inv,tstr,fwe,write_neg_pos)
if fwe
    fh0P = H.Pt;
    fh0N = H.Nt;
    fh0C = H.Ct;
else
    fh0P = H.Pu;
    fh0N = H.Nu;
    fh0C = H.Cu;
end
%Begin with no colorbars
load split
cbar1 = 0; cbar2 = 0;
if ~isempty(DF) %&& DF.GroupColorbars %Separated colorbars too far from working -- just skip
    try
        GInv = CF.GInv;
        nC = CF.nC;
        if isfield(DF,'fh1')
            fh1 = DF.fh1;
            ax1 = DF.ax1;
            if isfield(DF,'hc1')
                cbar1 = 1;
                %hc1 = DF.hc1;
                hc1_min = DF.hc1_min;
                hc1_max = DF.hc1_max;
            end
            split1 = DF.split1;
            %split1 = split;
        end
        tick_number = DF.tick_number-2;
        fontsize_choice = DF.fontsize_choice-5;
        if GInv
            if isfield(DF,'fh2')
                fh2 = DF.fh2;
                ax2 = DF.ax2;
                if DF.GroupColorbars
                    split2 = DF.split2;
                else
                    split1 = DF.split1;
                end
                if isfield(DF,'hc2')
                    cbar2 = 1;
                    %hc2 = DF.hc2;
                    if isfield(DF,'hc2_min')
                        hc2_min = DF.hc2_min;
                        hc2_max = DF.hc2_max;
                    end
                end
            end
        end
        
        %Assemble figures onto 2 grids (one for positive responses and one for negatives responses)
        %of subplots of size number of contrasts by 3 (for HbR, HbO and HbT)
        
        %GInv: Boolean to specify if inverted (i.e. negative) responses are
        %generated at all. Note that a decrease in HbR is considered a positive
        %response
        %Inv: Boolean to specify whether an increase or decrease chart is the
        %current one looked at - if Inv = 1, it is an increase (perhaps
        %confusingly), and if Inv = 0, it is a decrease
        
        %Assumptions:
        %1st contrast is the 1st Volterra kernel, 2nd contrast is the 2nd (T-contrasts)
        %3rd contrast is an F contrast for 1st Volterra, 4th contrast is an
        %F-contrast for 2nd Volterra
        if write_neg_pos
            if exist('fh1','var')
                if GInv %Both positive and negative responses generated
                    if Inv %Increase if 1
                        if mod(c1,2) == 1 %1st Volterra
                            fhP = fh0P; %Positive response: (i.e. HbR-, HbO+, HbT+), depending on Inv = 1 or 0
                            fhN = fh0N;
                        else %2nd Volterra - looking for inverted response only for T contrast
                            fhP = fh0N; %Looking for inhibitive or saturation response
                            fhN = fh0P;
                        end
                    else %Decrease if 0
                        if mod(c1,2) == 1 %1st Volterra
                            fhN = fh0P; %
                            fhP = fh0N;
                        else
                            fhN = fh0N;
                            fhP = fh0P;
                        end
                        
                    end
                    switch hb
                        case 'HbO'
                            figure(fhP);
                            axNew=subplot(nC,3,3*(c1-1)+2);
                        case 'HbT'
                            figure(fhP);
                            axNew=subplot(nC,3,3*(c1-1)+3);
                        case 'HbR'
                            figure(fhN);
                            axNew=subplot(nC,3,3*(c1-1)+1);
                    end
                else %Only positive responses (i.e. HbR-, HbO+, HbT+)
                    if Inv %This must be HbO+, HbT+
                        fhP = fh0P;
                        fhN = fh0N;
                    else %This must be HbR-
                        fhP = fh0N;
                        fhN = fh0P;
                    end
                    switch hb
                        case 'HbO' %2nd column
                            figure(fhP);
                            axNew=subplot(nC,3,3*(c1-1)+2);
                        case 'HbT' %3rd column
                            figure(fhP);
                            axNew=subplot(nC,3,3*(c1-1)+3);
                        case 'HbR' %^first column
                            figure(fhN);
                            axNew=subplot(nC,3,3*(c1-1)+1);
                    end
                end
                set(fh1,'CurrentAxes',ax1);
                copyobj(allchild(ax1),axNew);
                colormap(split)
                %colormap(split1)
                axis(axNew, 'off')
                axis(axNew, 'image')
                axis(axNew, 'ij') %otherwise bottom and up are inverted
                %copyobj(allchild(hc1),hcNew);
                if cbar1
                    hcNew = colorbar;
                
                
                %make the subplots larger
                p = get(axNew, 'pos');
                p(3) = p(3)+0.05;
                p(4) = p(4)+0.05;
                set(axNew, 'pos', p);
                freezeColors(axNew);
                %caxis manual
                set(hcNew,'YLim', [hc1_min hc1_max]);
                y_tick = linspace(hc1_min, hc1_max, tick_number)';
                set(hcNew, 'YTick', y_tick);
                set(hcNew, 'FontSize', fontsize_choice);
                %Customize here number of decimals
                set(hcNew,'YTickLabel',sprintf('%.1f |',get(hcNew,'YTick')'));
                end
                %cbfreeze(hcNew);
                
                if ~GInv || (tstr == 'T' && GInv)
                    close(fh1);
                    %otherwise, don't close the figure as it is needed below (for F
                    %contrasts only
                end
            end
        end
        %Repeat for combined positive + negative figure:
        if exist('fh2','var')
            if GInv && Inv
                switch hb
                    case 'HbO'
                        figure(fh0C); %actual increases or decreases, no sign reversals
                        axNew2=subplot(nC,3,3*(c1-1)+2);
                    case 'HbT'
                        figure(fh0C);
                        axNew2=subplot(nC,3,3*(c1-1)+3);
                    case 'HbR'
                        figure(fh0C);
                        axNew2=subplot(nC,3,3*(c1-1)+1);
                end
                set(fh2,'CurrentAxes',ax2);
                copyobj(allchild(ax2),axNew2);
                if DF.GroupColorbars
                    colormap(axNew2,split2);
                else
                    colormap(axNew2,split1);
                end
                
                axis(axNew2, 'off')
                axis(axNew2, 'image')
                axis(axNew2, 'ij') %otherwise bottom and up are inverted
                %                 if cbar1
                %                     clear hcNew
                % %                     if ~DF.GroupColorbars
                % %                         hcNew = colorbar('EastOutside');
                % %                     else
                %                         hcNew = colorbar;
                % %                     end
                %                 end
                %                 if cbar2 && ~cbar1
                % %                     if ~DF.GroupColorbars
                % %                         hcNew2 = colorbar('WestOutside');
                % %                     else
                %                         hcNew2 = colorbar;
                % %                     end
                %                 end
                if cbar2 
                hcNew2 = colorbar;
                end
                %make the subplots larger
                if DF.GroupColorbars
                    p = get(axNew2, 'pos');
                    p(3) = p(3)+0.05;
                    p(4) = p(4)+0.05;
                    set(axNew2, 'pos', p);
                end
                drawnow
                freezeColors(axNew2)
                %caxis manual
                %                 if cbar1
                %                     if ~DF.GroupColorbars
                %                         set(hcNew,'YLim', [hc1_min hc1_max]);
                %                         y_tick = linspace(hc1_min, hc1_max, tick_number)';
                %                         set(hcNew, 'YTick', y_tick);
                %                         set(hcNew, 'FontSize', fontsize_choice);
                %                         %Customize here number of decimals
                %                         set(hcNew,'YTickLabel',sprintf('%.1f |',get(hcNew,'YTick')'));
                %                         %hold on
                %                         %caxis([hc1_min hc1_max]);
                %                         %v1 = caxis(hcNew);
                %                         cbfreeze(hcNew);
                %                     end
                %                 end
                
                %                 if cbar2
                %if DF.GroupColorbars
                %                     set(hcNew2,'YLim', [hc2_min hc2_max]);
                %                     y_tick = linspace(hc2_min, hc2_max, tick_number)';
                if ~exist('hc2_min','var')
                    hc2_min = hc1_min;
                end
                if ~exist('hc1_max','var')
                    hc1_max = hc2_max;
                end
                if cbar2
                set(hcNew2,'YLim', [hc2_min hc1_max]);
                y_tick = linspace(hc2_min, hc1_max, tick_number)';
                set(hcNew2, 'YTick', y_tick);
                set(hcNew2, 'FontSize', fontsize_choice);
                %Customize here number of decimals
                set(hcNew2,'YTickLabel',sprintf('%.1f |',get(hcNew2,'YTick')'));
                %v2 = caxis(hcNew2);
                %hold on
                %caxis([hc2_min hc2_max+(hc2_max-hc2_min)])
                cbfreeze(hcNew2);
                end
                %                 end
                %
                close(fh2); %For F contrasts, this will be the same as fh1
            end
        end
        
        %Repeat for combined positive + negative figure:
        if ~exist('fh2','var') && strcmp(tstr,'F')
            %if GInv && Inv
            figure(fh0C); %actual increases or decreases, no sign reversals
            switch hb
                case 'HbO'
                    axNew=subplot(nC,3,3*(c1-1)+2);
                case 'HbT'
                    axNew=subplot(nC,3,3*(c1-1)+3);
                case 'HbR'
                    axNew=subplot(nC,3,3*(c1-1)+1);
            end
            set(fh1,'CurrentAxes',ax1);
            copyobj(allchild(ax1),axNew);
            if DF.GroupColorbars
                colormap(axNew,split1);
            else
                colormap(axNew,split1);
            end
            
            axis(axNew, 'off')
            axis(axNew, 'image')
            axis(axNew, 'ij') %otherwise bottom and up are inverted
            %                 if cbar1
            %                     clear hcNew
            % %                     if ~DF.GroupColorbars
            % %                         hcNew = colorbar('EastOutside');
            % %                     else
            %                         hcNew = colorbar;
            % %                     end
            %                 end
            %                 if cbar2 && ~cbar1
            % %                     if ~DF.GroupColorbars
            % %                         hcNew2 = colorbar('WestOutside');
            % %                     else
            %                         hcNew2 = colorbar;
            % %                     end
            %                 end
            hcNew = colorbar;
            %make the subplots larger
            if DF.GroupColorbars
                p = get(axNew, 'pos');
                p(3) = p(3)+0.05;
                p(4) = p(4)+0.05;
                set(axNew, 'pos', p);
            end
            drawnow
            freezeColors(axNew)
            %caxis manual
            %                 if cbar1
            %                     if ~DF.GroupColorbars
            %                         set(hcNew,'YLim', [hc1_min hc1_max]);
            %                         y_tick = linspace(hc1_min, hc1_max, tick_number)';
            %                         set(hcNew, 'YTick', y_tick);
            %                         set(hcNew, 'FontSize', fontsize_choice);
            %                         %Customize here number of decimals
            %                         set(hcNew,'YTickLabel',sprintf('%.1f |',get(hcNew,'YTick')'));
            %                         %hold on
            %                         %caxis([hc1_min hc1_max]);
            %                         %v1 = caxis(hcNew);
            %                         cbfreeze(hcNew);
            %                     end
            %                 end
            
            %                 if cbar2
            %if DF.GroupColorbars
            %                     set(hcNew2,'YLim', [hc2_min hc2_max]);
            %                     y_tick = linspace(hc2_min, hc2_max, tick_number)';
            %                 if ~exist('hc2_min','var')
            %                     hc2_min = hc1_min;
            %                 end
            %                 if ~exist('hc1_max','var')
            %                     hc1_max = hc2_max;
            %                 end
            set(hcNew,'YLim', [hc1_min hc1_max]);
            y_tick = linspace(hc1_min, hc1_max, tick_number)';
            set(hcNew, 'YTick', y_tick);
            set(hcNew, 'FontSize', fontsize_choice);
            %Customize here number of decimals
            set(hcNew,'YTickLabel',sprintf('%.1f |',get(hcNew,'YTick')'));
            %v2 = caxis(hcNew2);
            %hold on
            %caxis([hc2_min hc2_max+(hc2_max-hc2_min)])
            %cbfreeze(hcNew);
            %                 end
            %
            close(fh1); %For F contrasts, this will be the same as fh1
        end
        %end
        if fwe
            H.Pt = fh0P;
            H.Nt = fh0N;
            H.Ct = fh0C;
        else
            H.Pu = fh0P;
            H.Nu = fh0N;
            H.Cu = fh0C;
        end
        try %just to see if that solves the memory leak
            close(fh1);
            close(fh2);
        end
        try %try to remove all dangling objects
            
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
    end
end
end

