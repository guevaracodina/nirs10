function hc = nirs_set_colorbar(hc,hc_min,hc_max,tick_number,fontsize_choice)
if hc_min==hc_max %quick fix
    if hc_min < 0
        hc_min=1.001*hc_max;
    else
        hc_min=0.999*hc_max; 
    end
end
set(hc, 'YLim', [hc_min hc_max]);
y_tick = linspace(hc_min, hc_max, tick_number)';
set(hc, 'YTick', y_tick);
set(hc, 'FontSize', fontsize_choice);
%Customize here number of decimals
set(hc,'YTickLabel',sprintf('%.1f |',get(hc,'YTick')'));