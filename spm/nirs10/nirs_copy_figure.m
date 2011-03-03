function [fh0P,fh0N] = nirs_copy_figure(fh0P,fh0N,c1,nC,hb,GInv,fh1,ax1,hc1,Inv,split)
if GInv 
    if Inv
        if mod(c1,2) == 1
            fhP = fh0P;
            fhN = fh0N;
        else
            fhP = fh0N;
            fhN = fh0P;
        end
    else
        if mod(c1,2) == 1
            fhN = fh0P;
            fhP = fh0N;
        else
            fhN = fh0N;
            fhP = fh0P;
        end
        
    end
    switch hb
        case 'HbO'
            figure(fhP);
            ax2=subplot(nC,3,3*(c1-1)+2);
        case 'HbT'
            figure(fhP);
            ax2=subplot(nC,3,3*(c1-1)+3);
        case 'HbR'
            figure(fhN);
            ax2=subplot(nC,3,3*(c1-1)+1);

    end
else 
    if Inv 
        fhP = fh0P;
        fhN = fh0N;
    else
        fhP = fh0N;
        fhN = fh0P;
    end
    switch hb
        case 'HbO'
            figure(fhP);
            ax2=subplot(nC,3,3*(c1-1)+2);
        case 'HbT'
            figure(fhP);            
            ax2=subplot(nC,3,3*(c1-1)+3);
        case 'HbR'
            figure(fhN);
            ax2=subplot(nC,3,3*(c1-1)+1);
    end
end
try
    set(fh1,'CurrentAxes',ax1);
    copyobj(allchild(ax1),ax2);
    colormap(split)
    hc2 = colorbar;
    axis(ax2, 'off')
    axis(ax2, 'image')
    axis(ax2, 'ij') %otherwise bottom and up are inverted
    copyobj(allchild(hc1),hc2);
    %make the subplots larger
    p = get(ax2, 'pos');
    p(3) = p(3)+0.05;
    p(4) = p(4)+0.05;
    set(ax2, 'pos', p);
end
try close(fh1); end
end

