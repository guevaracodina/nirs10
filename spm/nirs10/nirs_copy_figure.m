function [fh0P,fh0N,fh0C] = nirs_copy_figure(fh0P,fh0N,fh0C,DF,CF,c1,hb,Inv,tstr)
try
GInv = CF.GInv;
nC = CF.nC;
%split = CF.split;
fh1 = DF.fh1;
ax1 = DF.ax1;
hc1 = DF.hc1;
split1 = DF.split1;
hc1_min = DF.hc1_min;
hc1_max = DF.hc1_max;
tick_number = DF.tick_number-2;
fontsize_choice = DF.fontsize_choice-5;
if GInv
    fh2 = DF.fh2;
    ax2 = DF.ax2;
    hc2 = DF.hc2;
    split2 = DF.split2;
    hc2_min = DF.hc2_min;
    hc2_max = DF.hc2_max;
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
    colormap(split1)
    hcNew = colorbar;
    axis(axNew, 'off')
    axis(axNew, 'image')
    axis(axNew, 'ij') %otherwise bottom and up are inverted
    %copyobj(allchild(hc1),hcNew);
    set(hcNew,'YLim', [hc1_min hc1_max]);
    y_tick = linspace(hc1_min, hc1_max, tick_number)';
    set(hcNew, 'YTick', y_tick);
    set(hcNew, 'FontSize', fontsize_choice);
    %Customize here number of decimals
    set(hcNew,'YTickLabel',sprintf('%.1f |',get(hcNew,'YTick')'));  
    
    %make the subplots larger
    p = get(axNew, 'pos');
    p(3) = p(3)+0.05;
    p(4) = p(4)+0.05;
    set(axNew, 'pos', p);
    
    %cbfreeze(hcNew);
    if ~GInv || (tstr == 'T' && GInv) 
        close(fh1);
        %otherwise, don't close the figure as it is needed below (for F
        %contrasts only
    end

%Repeat for combined positive + negative figure:
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
    colormap(axNew2,split2)
    hcNew2 = colorbar;
    axis(axNew2, 'off')
    axis(axNew2, 'image')
    axis(axNew2, 'ij') %otherwise bottom and up are inverted
    %copyobj(allchild(hc2),hcNew2); %doesn't work?
    set(hcNew2,'YLim', [hc2_min hc2_max]);
    y_tick = linspace(hc2_min, hc2_max, tick_number)';
    set(hcNew2, 'YTick', y_tick);
    set(hcNew2, 'FontSize', fontsize_choice);
    %Customize here number of decimals
    set(hcNew2,'YTickLabel',sprintf('%.1f |',get(hcNew2,'YTick')'));  

    %make the subplots larger
    p = get(axNew2, 'pos');
    p(3) = p(3)+0.05;
    p(4) = p(4)+0.05;
    set(axNew2, 'pos', p);
    drawnow
    freezeColors(axNew2);
    cbfreeze(hcNew2);
    close(fh2); %For F contrasts, this will be the same as fh1
end
catch exception
    disp(exception.identifier);
end
end

